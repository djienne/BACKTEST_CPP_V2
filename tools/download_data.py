#!/usr/bin/env python3
"""
Download Binance OHLCV and funding-rate history for the backtester.

Pulls from Binance's public archive (data.binance.vision) rather than the REST API:
the archive serves whole months as ZIPs, so there is no pagination and no rate limit,
and every file has a SHA256 checksum alongside it.

    python3 tools/download_data.py --check          # what's missing? downloads nothing
    python3 tools/download_data.py                  # fetch everything the config needs
    python3 tools/download_data.py --coins BTC --timeframes 5m --markets spot
    python3 tools/download_data.py --start 2021-01 --jobs 8 --force

Standard library only -- no pip install. Python 3.8+.

Output goes exactly where the C++ loaders look (see tools.cpp):

    spot     data/data/binance/<tf>/<COIN>-USDT.csv
    futures  data/data/futures/<COIN>_USDT-<tf>-futures.json
    funding  data/data/futures/<COIN>_USDT-8h-funding_rate.json

Downloaded ZIPs are cached under .data_cache/ so re-runs only fetch new months; delete
that directory to force a clean re-download.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import sys
import time
import urllib.error
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = REPO_ROOT / "backtest_config.json"
CACHE_DIR = REPO_ROOT / ".data_cache"

ARCHIVE = "https://data.binance.vision/data"
SPOT_API = "https://api.binance.com/api/v3/klines"
FUTURES_API = "https://fapi.binance.com/fapi/v1/klines"

# Binance kline CSVs carry 12 columns; we keep the first six.
#   open time | open | high | low | close | volume | close time | ...
KLINE_OPEN_TIME = 0
KLINE_VOLUME = 5

USER_AGENT = "BACKTEST_CPP_V2-downloader/1.0"
MAX_RETRIES = 4

# The four files the tracked regression fixtures are computed from. Binance restates
# historical spot klines: re-downloading BTC 1h reproduces 50,960 of 50,976 rows exactly
# but 16 differ, 5 of them in `close` -- the field every strategy trades on. Overwriting
# these would therefore shift regression_expected/ without anyone asking for it, so
# --force skips them unless --force-fixtures is also given.
#
# (Futures klines and funding rates matched the tracked copies exactly, so the
# restatement seems confined to spot.)
FIXTURE_PINNED = {
    "binance/1h/BTC-USDT.csv",
    "binance/1h/ETH-USDT.csv",
    "futures/BTC_USDT-1h-futures.json",
    "futures/BTC_USDT-8h-funding_rate.json",
}


# --------------------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------------------


class DownloadError(RuntimeError):
    pass


def log(msg: str) -> None:
    print(msg, flush=True)


def month_range(start: str, end_exclusive: datetime) -> List[Tuple[int, int]]:
    """Every (year, month) from `start` ("YYYY-MM") up to but excluding end_exclusive."""
    year, month = (int(p) for p in start.split("-"))
    out: List[Tuple[int, int]] = []
    while (year, month) < (end_exclusive.year, end_exclusive.month):
        out.append((year, month))
        month += 1
        if month == 13:
            year, month = year + 1, 1
    return out


def days_in_month(year: int, month: int) -> int:
    nxt = datetime(year + (month == 12), (month % 12) + 1, 1, tzinfo=timezone.utc)
    return (nxt - timedelta(days=1)).day


def normalize_ms(raw: int) -> int:
    """
    Binance switched the archive timestamp unit from milliseconds to **microseconds**
    with the 2025-01 files. Months either side of that boundary land in the same output
    file, so mixing them unconverted sorts the series wrongly and puts dates in the year
    58595 -- which is exactly how this surfaced.

    There is no unit marker in the CSV, so it is inferred from magnitude. In milliseconds
    the years 2000-2100 span roughly 0.95e12 .. 4.1e12; in microseconds they start above
    9.4e14. A 1e14 threshold sits an order of magnitude clear of both.
    """
    return raw // 1000 if raw > 100_000_000_000_000 else raw


def trim_number(text: str) -> str:
    """
    "4261.48000000" -> "4261.48",  "47.18100900" -> "47.181009",  "0.00000000" -> "0".

    Done on the string rather than via float() so no precision is lost, and so the output
    matches the existing tracked CSVs, which are stored in this trimmed form.
    """
    if "." not in text:
        return text
    trimmed = text.rstrip("0").rstrip(".")
    return trimmed if trimmed not in ("", "-") else "0"


def http_get(url: str, timeout: int = 60) -> bytes:
    """GET with retries. Raises FileNotFoundError on 404 so callers can skip missing months."""
    last: Optional[Exception] = None
    for attempt in range(MAX_RETRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise FileNotFoundError(url) from exc
            last = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
        # Back off before retrying; the archive is a CDN and occasionally 5xx's.
        time.sleep(1.5 * (attempt + 1))
    raise DownloadError(f"{url}: {last}")


def write_atomic(path: Path, payload: str) -> None:
    """
    Write via a temp file then rename, so an interrupted run can never leave a truncated
    dataset behind -- the C++ loaders would happily read a half-written file.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(payload, encoding="utf-8", newline="\n")
    os.replace(tmp, path)


# --------------------------------------------------------------------------------------
# Archive access
# --------------------------------------------------------------------------------------


def archive_url(market: str, symbol: str, kind: str, timeframe: str, year: int, month: int,
                day: Optional[int] = None) -> str:
    period = "daily" if day else "monthly"
    stamp = f"{year:04d}-{month:02d}" + (f"-{day:02d}" if day else "")
    if kind == "fundingRate":
        base = f"{ARCHIVE}/futures/um/{period}/fundingRate/{symbol}"
        return f"{base}/{symbol}-fundingRate-{stamp}.zip"
    root = "spot" if market == "spot" else "futures/um"
    base = f"{ARCHIVE}/{root}/{period}/klines/{symbol}/{timeframe}"
    return f"{base}/{symbol}-{timeframe}-{stamp}.zip"


def fetch_zip(url: str, cache_key: Path, verify: bool = True) -> Optional[bytes]:
    """Return the ZIP bytes, from cache when possible. None if the archive has no such file."""
    if cache_key.exists():
        return cache_key.read_bytes()

    try:
        blob = http_get(url)
    except FileNotFoundError:
        return None

    if verify:
        try:
            checksum_blob = http_get(url + ".CHECKSUM").decode()
            # sha256sum -c format: "<hex>  <filename>"
            expected = checksum_blob.split()[0].lower()
            actual = hashlib.sha256(blob).hexdigest()
            if actual != expected:
                raise DownloadError(f"{url}: checksum mismatch ({actual} != {expected})")
        except FileNotFoundError:
            pass  # older archives lack a checksum sibling

    cache_key.parent.mkdir(parents=True, exist_ok=True)
    cache_key.write_bytes(blob)
    return blob


def rows_from_zip(blob: bytes, expected_columns: int) -> List[List[str]]:
    """
    Extract CSV rows from a Binance archive ZIP.

    Archives changed at some point to include a header row, so rather than assuming
    either way, a first cell that is not an integer is treated as a header and dropped.
    """
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        name = zf.namelist()[0]
        text = zf.read(name).decode("utf-8", errors="replace")

    rows: List[List[str]] = []
    for row in csv.reader(io.StringIO(text)):
        if len(row) < expected_columns:
            continue
        try:
            int(row[0])
        except ValueError:
            continue  # header line
        rows.append(row)
    return rows


# --------------------------------------------------------------------------------------
# Series assembly
# --------------------------------------------------------------------------------------


@dataclass
class Series:
    rows: List[List[str]]          # [open_time_ms, o, h, l, c, v] as strings
    missing_months: List[str]
    gaps: int


def collect_klines(market: str, symbol: str, timeframe: str, months: Sequence[Tuple[int, int]],
                   verify: bool, quiet: bool) -> Series:
    by_time: Dict[int, List[str]] = {}
    missing: List[str] = []

    for (year, month) in months:
        url = archive_url(market, symbol, "klines", timeframe, year, month)
        cache = CACHE_DIR / market / symbol / timeframe / f"{year:04d}-{month:02d}.zip"
        blob = fetch_zip(url, cache, verify)

        if blob is None:
            # The monthly file is only published once the month is complete; fall back to
            # the daily archives, which is what makes "today" reachable.
            day_rows: List[List[str]] = []
            for day in range(1, days_in_month(year, month) + 1):
                d_url = archive_url(market, symbol, "klines", timeframe, year, month, day)
                d_cache = CACHE_DIR / market / symbol / timeframe / f"{year:04d}-{month:02d}-{day:02d}.zip"
                d_blob = fetch_zip(d_url, d_cache, verify)
                if d_blob is not None:
                    day_rows.extend(rows_from_zip(d_blob, 6))
            if not day_rows:
                missing.append(f"{year:04d}-{month:02d}")
                continue
            rows = day_rows
        else:
            rows = rows_from_zip(blob, 6)

        for row in rows:
            open_time = normalize_ms(int(row[KLINE_OPEN_TIME]))
            by_time[open_time] = [
                str(open_time),
                *[trim_number(row[i]) for i in range(1, KLINE_VOLUME + 1)],
            ]

    ordered = [by_time[t] for t in sorted(by_time)]
    gaps = count_gaps([int(r[0]) for r in ordered], timeframe)
    if gaps and not quiet:
        log(f"    note: {gaps} gap(s) in {symbol} {timeframe} -- missing candles, not an error")
    return Series(ordered, missing, gaps)


def collect_funding(symbol: str, months: Sequence[Tuple[int, int]], verify: bool) -> Series:
    by_time: Dict[int, List[str]] = {}
    missing: List[str] = []

    for (year, month) in months:
        url = archive_url("futures", symbol, "fundingRate", "", year, month)
        cache = CACHE_DIR / "funding" / symbol / f"{year:04d}-{month:02d}.zip"
        blob = fetch_zip(url, cache, verify)
        if blob is None:
            missing.append(f"{year:04d}-{month:02d}")
            continue
        # fundingRate CSV: calc_time, funding_interval_hours, last_funding_rate
        for row in rows_from_zip(blob, 3):
            calc_time = normalize_ms(int(row[0]))
            by_time[calc_time] = [str(calc_time), row[-1]]

    ordered = [by_time[t] for t in sorted(by_time)]
    return Series(ordered, missing, 0)


TIMEFRAME_MS = {
    "1m": 60_000, "3m": 180_000, "5m": 300_000, "15m": 900_000, "30m": 1_800_000,
    "1h": 3_600_000, "2h": 7_200_000, "4h": 14_400_000, "6h": 21_600_000,
    "8h": 28_800_000, "12h": 43_200_000, "1d": 86_400_000,
}


def count_gaps(times: Sequence[int], timeframe: str) -> int:
    step = TIMEFRAME_MS.get(timeframe)
    if not step or len(times) < 2:
        return 0
    return sum(1 for a, b in zip(times, times[1:]) if b - a != step)


# --------------------------------------------------------------------------------------
# Writers -- these must match tools.cpp exactly
# --------------------------------------------------------------------------------------


def write_spot_csv(path: Path, series: Series) -> None:
    lines = ["date,open,high,low,close,volume"]
    lines.extend(",".join(row) for row in series.rows)
    write_atomic(path, "\n".join(lines) + "\n")


def write_futures_json(path: Path, series: Series) -> None:
    payload = [
        [int(r[0]), float(r[1]), float(r[2]), float(r[3]), float(r[4]), float(r[5])]
        for r in series.rows
    ]
    write_atomic(path, json.dumps(payload, separators=(",", ":")))


def write_funding_json(path: Path, series: Series) -> None:
    # Six elements in freqtrade's shape: the loader reads [0] and [1], and the trailing
    # zeros keep the file interchangeable with freqtrade's own output.
    payload = [[int(r[0]), float(r[1]), 0.0, 0.0, 0.0, 0.0] for r in series.rows]
    write_atomic(path, json.dumps(payload, separators=(",", ":")))


# --------------------------------------------------------------------------------------
# Config
# --------------------------------------------------------------------------------------


@dataclass
class Job:
    market: str        # "spot" | "futures" | "funding"
    coin: str
    timeframe: str

    @property
    def symbol(self) -> str:
        return f"{self.coin}USDT"

    def output(self, data_dir: Path) -> Path:
        if self.market == "spot":
            return data_dir / "binance" / self.timeframe / f"{self.coin}-USDT.csv"
        if self.market == "futures":
            return data_dir / "futures" / f"{self.coin}_USDT-{self.timeframe}-futures.json"
        return data_dir / "futures" / f"{self.coin}_USDT-8h-funding_rate.json"

    def label(self) -> str:
        return f"{self.market:<7} {self.coin:<5} {self.timeframe}"


def load_config(path: Path) -> dict:
    if not path.exists():
        raise SystemExit(f"config not found: {path}\nExpected {DEFAULT_CONFIG.name} at the repo root.")
    with path.open(encoding="utf-8") as fh:
        return json.load(fh)


def jobs_from_config(cfg: dict) -> List[Job]:
    """
    One job per (market, coin, timeframe) any strategy needs. Derived from the same file
    the strategies read, so the download set cannot drift from what they open.
    """
    coins: List[str] = cfg["coins"]
    seen = set()
    jobs: List[Job] = []

    def add(market: str, timeframe: str) -> None:
        for coin in coins:
            key = (market, coin, timeframe)
            if key not in seen:
                seen.add(key)
                jobs.append(Job(market, coin, timeframe))

    for name, entry in cfg.get("strategies", {}).items():
        market = entry.get("market", "spot")
        for tf_key in ("timeframe", "htf"):
            tf = entry.get(tf_key)
            if tf:
                add(market, tf)
        if market == "futures":
            add("funding", "8h")

    return jobs


# --------------------------------------------------------------------------------------
# Commands
# --------------------------------------------------------------------------------------


def cmd_check(jobs: Sequence[Job], data_dir: Path, cfg: dict) -> int:
    log(f"Coverage against {data_dir}\n")
    missing = 0
    for job in sorted(jobs, key=lambda j: (j.market, j.coin, j.timeframe)):
        out = job.output(data_dir)
        if out.exists():
            size = out.stat().st_size
            log(f"  present  {job.label()}  {size / 1e6:8.1f} MB  {out.relative_to(REPO_ROOT)}")
        else:
            missing += 1
            log(f"  MISSING  {job.label()}  {out.relative_to(REPO_ROOT)}")

    log("")
    strategies = cfg.get("strategies", {})
    coins: List[str] = cfg["coins"]
    blocked = []
    for name, entry in sorted(strategies.items()):
        market = entry.get("market", "spot")
        needed: List[Job] = []
        for tf_key in ("timeframe", "htf"):
            if entry.get(tf_key):
                needed += [Job(market, c, entry[tf_key]) for c in coins]
        if market == "futures":
            needed += [Job("funding", c, "8h") for c in coins]
        absent = [j for j in needed if not j.output(data_dir).exists()]
        if absent:
            blocked.append((name, absent))

    if blocked:
        log(f"{len(blocked)} of {len(strategies)} strategies cannot run:")
        for name, absent in blocked:
            what = ", ".join(sorted({f"{j.market} {j.timeframe}" for j in absent}))
            log(f"  {name:<45} needs {what}")
    else:
        log(f"All {len(strategies)} strategies have their data.")

    return 1 if missing else 0


def is_fixture_pinned(job: Job, data_dir: Path) -> bool:
    try:
        return job.output(data_dir).relative_to(data_dir).as_posix() in FIXTURE_PINNED
    except ValueError:
        return False


def run_job(job: Job, data_dir: Path, months: Sequence[Tuple[int, int]], verify: bool,
            force: bool, quiet: bool, force_fixtures: bool = False) -> str:
    out = job.output(data_dir)
    if out.exists() and not force:
        return f"  skip     {job.label()}  (exists; --force to rebuild)"
    if out.exists() and force and is_fixture_pinned(job, data_dir) and not force_fixtures:
        return (f"  PINNED   {job.label()}  regression fixtures are computed from this file; "
                f"Binance restates history, so rebuilding it would shift them "
                f"(--force-fixtures to override)")

    if job.market == "funding":
        series = collect_funding(job.symbol, months, verify)
        if not series.rows:
            return f"  EMPTY    {job.label()}  no funding data in range"
        write_funding_json(out, series)
    else:
        series = collect_klines(job.market, job.symbol, job.timeframe, months, verify, quiet)
        if not series.rows:
            return f"  EMPTY    {job.label()}  archive has nothing in range"
        if job.market == "spot":
            write_spot_csv(out, series)
        else:
            write_futures_json(out, series)

    first = datetime.fromtimestamp(int(series.rows[0][0]) / 1000, timezone.utc).date()
    last = datetime.fromtimestamp(int(series.rows[-1][0]) / 1000, timezone.utc).date()
    return f"  ok       {job.label()}  {len(series.rows):>8,} rows  {first} .. {last}"


def cmd_download(jobs: Sequence[Job], data_dir: Path, start: str, verify: bool, force: bool,
                 jobs_n: int, quiet: bool, force_fixtures: bool = False) -> int:
    now = datetime.now(timezone.utc)
    # Include the current month; collect_klines falls back to daily files for it.
    months = month_range(start, now + timedelta(days=31))
    log(f"Downloading {len(jobs)} series from {start}, {jobs_n} parallel workers")
    log(f"Cache: {CACHE_DIR.relative_to(REPO_ROOT)}  (delete to force a clean re-download)\n")

    failures = 0
    with ThreadPoolExecutor(max_workers=jobs_n) as pool:
        futures = {pool.submit(run_job, j, data_dir, months, verify, force, quiet, force_fixtures): j
                   for j in jobs}
        for fut in as_completed(futures):
            job = futures[fut]
            try:
                log(fut.result())
            except Exception as exc:  # noqa: BLE001 - report and keep going
                failures += 1
                log(f"  FAILED   {job.label()}  {exc}")

    log("")
    log("Done." if not failures else f"Done with {failures} failure(s).")
    return 1 if failures else 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Download Binance OHLCV/funding history for the backtester.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Standard library only")[0].strip(),
    )
    ap.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    ap.add_argument("--check", action="store_true", help="report coverage and exit")
    ap.add_argument("--coins", help="comma-separated override, e.g. BTC,ETH")
    ap.add_argument("--timeframes", help="comma-separated override, e.g. 5m,1h")
    ap.add_argument("--markets", help="comma-separated subset of spot,futures,funding")
    ap.add_argument("--start", help='first month, "YYYY-MM" (default: config history_start)')
    ap.add_argument("--jobs", type=int, default=4,
                    help="parallel downloads (default 4; network-bound, raise if your link allows)")
    ap.add_argument("--force", action="store_true", help="rebuild outputs that already exist")
    ap.add_argument("--force-fixtures", action="store_true",
                    help="also rebuild the four files the regression fixtures depend on "
                         "(this will shift regression_expected/ -- rebaseline deliberately)")
    ap.add_argument("--no-verify", action="store_true", help="skip SHA256 checksum verification")
    ap.add_argument("--quiet", action="store_true", help="suppress per-series gap notes")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    if args.coins:
        cfg["coins"] = [c.strip().upper() for c in args.coins.split(",") if c.strip()]

    data_dir = (REPO_ROOT / cfg.get("data_dir", "./data/data")).resolve()
    jobs = jobs_from_config(cfg)

    if args.markets:
        wanted = {m.strip() for m in args.markets.split(",")}
        jobs = [j for j in jobs if j.market in wanted]
    if args.timeframes:
        wanted = {t.strip() for t in args.timeframes.split(",")}
        jobs = [j for j in jobs if j.timeframe in wanted or j.market == "funding"]

    if not jobs:
        log("Nothing to do -- filters excluded every series.")
        return 0

    if args.check:
        return cmd_check(jobs, data_dir, cfg)

    start = args.start or cfg.get("history_start", "2017-08")
    return cmd_download(jobs, data_dir, start, not args.no_verify, args.force,
                        max(1, args.jobs), args.quiet, args.force_fixtures)


if __name__ == "__main__":
    sys.exit(main())
