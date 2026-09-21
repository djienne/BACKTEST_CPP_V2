#!/usr/bin/env python3
"""Validate and repair Binance research data before a run.

Archives are the first source; public REST endpoints repair remaining holes.
Only completed UTC candles are used. An unresolved hole is an error, never a
synthetic price. Test fixtures are separate from this mutable research cache.
"""
from __future__ import annotations
import argparse
import calendar
import csv
import hashlib
import io
import json
import math
import os
import sys
import shutil
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CONFIG = ROOT / "backtest_config.json"
ARCHIVE = "https://data.binance.vision/data"
API = {"spot": "https://api.binance.com/api/v3", "futures": "https://fapi.binance.com/fapi/v1"}
TIMEFRAME_MS = {f"{n}m": n * 60_000 for n in (1,3,5,15,30)}
TIMEFRAME_MS.update({f"{n}h": n * 3_600_000 for n in (1,2,4,6,8,12)})
TIMEFRAME_MS["1d"] = 86_400_000
CACHE_DIR = ROOT / ".data_cache"

class DownloadError(RuntimeError):
    pass

def log(message):
    print(message, flush=True)

def utc_ms(value):
    """Parse an ISO date/time as UTC; reject non-UTC offsets."""
    if len(value) == 7:
        value += "-01"
    dt = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if dt.tzinfo is not None and dt.utcoffset().total_seconds() != 0:
        raise DownloadError("Dates must be UTC")
    return int(dt.replace(tzinfo=timezone.utc).timestamp() * 1000)

def normalize_ms(raw):
    return int(raw) // 1000 if int(raw) > 100_000_000_000_000 else int(raw)

def http_get(url, timeout=60):
    """Bounded retries, including server-directed backoff. No credentials required."""
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "BACKTEST_CPP_V2/2"}), timeout=timeout) as response:
                return response.read()
        except urllib.error.HTTPError as error:
            if error.code == 404:
                raise FileNotFoundError(url) from error
            if error.code in (400,401,403,451):
                raise DownloadError(f"{url}: HTTP {error.code}") from error
            delay = min(60, int(error.headers.get("Retry-After", 2 ** attempt)))
            last = error
        except (urllib.error.URLError, TimeoutError, ConnectionError) as error:
            delay, last = 2 ** attempt, error
        if attempt < 3:
            time.sleep(delay)
    raise DownloadError(f"{url}: {last}")

def write_atomic(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".{os.getpid()}.tmp")
    try:
        tmp.write_text(payload, encoding="utf-8", newline="\n")
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)

def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def months(start, end):
    dt = datetime.fromtimestamp(start / 1000, timezone.utc)
    year, month = dt.year, dt.month
    while utc_ms(f"{year:04}-{month:02}") < end:
        yield year, month
        year, month = (year + 1, 1) if month == 12 else (year, month + 1)

def fetch_zip(url, key):
    key = Path(key)
    checksum = key.with_suffix(".sha256")
    if key.exists() and checksum.exists() and digest(key) == checksum.read_text().strip():
        return key.read_bytes()
    try:
        blob = http_get(url)
        expected = http_get(url + ".CHECKSUM").decode().split()[0].lower()
    except FileNotFoundError:
        return None
    if hashlib.sha256(blob).hexdigest() != expected:
        raise DownloadError(f"Checksum mismatch: {url}")
    # ZIP writes are atomic too; a cancelled transfer cannot poison the cache.
    key.parent.mkdir(parents=True, exist_ok=True)
    tmp = key.with_name(key.name + f".{os.getpid()}.tmp")
    try:
        tmp.write_bytes(blob)
        os.replace(tmp, key)
    finally:
        tmp.unlink(missing_ok=True)
    write_atomic(checksum, expected)
    return blob

def zip_rows(blob):
    with zipfile.ZipFile(io.BytesIO(blob)) as archive:
        names = [n for n in archive.namelist() if n.endswith(".csv")]
        if len(names) != 1:
            raise DownloadError("Expected one CSV in archive")
        rows = list(csv.reader(io.StringIO(archive.read(names[0]).decode("utf-8-sig"))))
    if rows and not rows[0][0].isdigit():
        rows = rows[1:]
    return rows

def validate_row(row, step):
    if len(row) < 6:
        raise DownloadError("Truncated OHLCV row")
    ts = normalize_ms(row[0])
    v = [float(x) for x in row[1:6]]
    if ts % step or not all(math.isfinite(x) for x in v):
        raise DownloadError(f"Invalid candle timestamp/value: {ts}")
    o,h,l,c,volume = v
    if min(o,h,l,c) <= 0 or volume < 0 or l > min(o,c) or h < max(o,c) or h < l:
        raise DownloadError(f"Impossible OHLCV candle: {ts}")
    return [ts, *v]

def merge_row(rows, row):
    """Accept repeated identical candles; reject competing values for one timestamp."""
    if row[0] in rows and rows[row[0]] != row:
        raise DownloadError(f"Conflicting duplicate candle at {row[0]}")
    rows[row[0]] = row


def read_rows(path, step):
    """Reject conflicting duplicates; identical duplicates can be safely collapsed."""
    path = Path(path)
    if not path.exists():
        return {}
    if path.suffix == ".csv":
        with path.open() as stream:
            rows = list(csv.reader(stream))
        if not rows or rows.pop(0) != ["date","open","high","low","close","volume"]:
            raise DownloadError(f"Invalid CSV header: {path}")
    else:
        rows = json.loads(path.read_text())
    out = {}
    for raw in rows:
        row = validate_row(raw, step)
        merge_row(out, row)
    return out

def missing_ranges(rows, start, end, step):
    """Return half-open contiguous holes without allocating an expected-timestamp set."""
    holes, cursor = [], start
    for ts in sorted(t for t in rows if start <= t < end):
        if ts > cursor:
            holes.append((cursor, ts))
        cursor = ts + step
    if cursor < end:
        holes.append((cursor, end))
    return holes

def kline_path(data, market, coin, tf):
    return (data / "binance" / tf / f"{coin}-USDT.csv" if market == "spot"
            else data / "futures" / f"{coin}_USDT-{tf}-futures.json")

def archive_candles(market, coin, tf, year, month, end):
    symbol = coin + "USDT"
    root = "spot" if market == "spot" else "futures/um"
    stamp = f"{year:04}-{month:02}"
    leaf = f"{symbol}-{tf}-{stamp}"
    url = f"{ARCHIVE}/{root}/monthly/klines/{symbol}/{tf}/{leaf}.zip"
    blob = fetch_zip(url, CACHE_DIR / market / coin / tf / (leaf + ".zip"))
    if blob is not None:
        return zip_rows(blob)
    # Daily archives cover the unpublished monthly tail. Older missing months
    # are repaired with REST rather than issuing 31 usually pointless requests.
    now = datetime.now(timezone.utc)
    if (year, month) != (now.year, now.month):
        return []
    rows = []
    for day in range(1, calendar.monthrange(year, month)[1] + 1):
        if utc_ms(f"{stamp}-{day:02}") >= end:
            break
        leaf = f"{symbol}-{tf}-{stamp}-{day:02}"
        url = f"{ARCHIVE}/{root}/daily/klines/{symbol}/{tf}/{leaf}.zip"
        blob = fetch_zip(url, CACHE_DIR / market / coin / tf / (leaf + ".zip"))
        if blob is not None:
            rows.extend(zip_rows(blob))
    return rows

def rest_candles(market, coin, tf, start, end):
    cursor = start
    while cursor < end:
        query = urllib.parse.urlencode(dict(symbol=coin+"USDT", interval=tf, startTime=cursor, endTime=end-1, limit=1000))
        page = json.loads(http_get(f"{API[market]}/klines?{query}"))
        if not isinstance(page, list):
            raise DownloadError(f"Invalid kline response for {coin}")
        if not page:
            break
        for row in page:
            if cursor <= normalize_ms(row[0]) < end:
                yield row
        following = normalize_ms(page[-1][0]) + TIMEFRAME_MS[tf]
        if following <= cursor:
            raise DownloadError("Kline pagination made no progress")
        cursor = following

def save_rows(path, rows):
    values = [rows[t] for t in sorted(rows)]
    if path.suffix == ".csv":
        stream = io.StringIO()
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["date","open","high","low","close","volume"])
        writer.writerows(values)
        write_atomic(path, stream.getvalue())
    else:
        write_atomic(path, json.dumps(values, separators=(",",":"), allow_nan=False))

def ensure_candles(data, market, coin, tf, start, end, offline=False, allow_prefix=False):
    step = TIMEFRAME_MS[tf]
    path = kline_path(data, market, coin, tf)
    try:
        rows = read_rows(path, step)
    except (ValueError, TypeError, DownloadError) as error:
        if offline:
            raise
        backup = path.with_name(path.name + f".invalid-{time.time_ns()}")
        shutil.copyfile(path, backup)
        log(f"Repairing invalid data in {path}: {error}; original saved as {backup}")
        rows = {}
    holes = missing_ranges(rows, start, end, step)
    if holes and offline:
        raise DownloadError(f"{coin} {tf}: missing {holes[0]} (offline)")
    if holes:
        log(f"Repairing {coin} {market} {tf}: {len(holes)} missing interval(s)")
        requested_months = sorted({m for a,b in holes for m in months(a,b)})
        for year, month in requested_months:
            for raw in archive_candles(market, coin, tf, year, month, end):
                row = validate_row(raw, step)
                if start <= row[0] < end:
                    merge_row(rows, row)
        for a,b in missing_ranges(rows, start, end, step):
            for raw in rest_candles(market, coin, tf, a,b):
                row = validate_row(raw, step)
                merge_row(rows, row)
        save_rows(path, rows)
    actual_start = max(start, min(rows)) if allow_prefix and rows else start
    holes = missing_ranges(rows, actual_start, end, step)
    if holes:
        a,b = holes[0]
        raise DownloadError(f"{coin} {market} {tf}: unresolved gap [{a},{b}); archives and REST could not supply it")
    selected = [rows[t] for t in sorted(rows) if actual_start <= t < end]
    if not selected:
        raise DownloadError(f"No candles for {coin} {tf}")
    return path, selected

def archive_funding(coin, start, end):
    """Published monthly rates and historical intervals, checked against REST marks."""
    out = {}
    for year, month in months(start, end):
        leaf = f"{coin}USDT-fundingRate-{year:04}-{month:02}"
        url = f"{ARCHIVE}/futures/um/monthly/fundingRate/{coin}USDT/{leaf}.zip"
        blob = fetch_zip(url, CACHE_DIR / "funding_archive" / (leaf + ".zip"))
        if blob is None:
            continue
        for raw in zip_rows(blob):
            ts, hours, rate = normalize_ms(raw[0]), float(raw[1]), float(raw[2])
            if not 0 < hours <= 8 or not math.isfinite(rate):
                raise DownloadError(f"Invalid funding archive: {coin} {ts}")
            if start <= ts < end:
                out[ts // 1000] = dict(timestamp_ms=ts, rate=rate, interval_hours=hours)
    return out


def ensure_funding(data, coin, start, end, offline=False):
    """Cache a fully paginated authoritative ledger, including settlement marks.

    Coverage + content digest distinguish 'no event' from a deleted/missing record.
    We do not infer an eight-hour schedule from the legacy filename.
    """
    path = data / "futures" / f"{coin}_USDT-funding.json"
    cached, good = None, False
    if path.exists():
        try:
            cached = json.loads(path.read_text())
            events = cached["events"]
            encoded = json.dumps(events, sort_keys=True, separators=(",",":"), allow_nan=False).encode()
            good = cached.get("schema") == 2 and cached["events_sha256"] == hashlib.sha256(encoded).hexdigest()
            if good:
                validate_funding(events)
                if not isinstance(cached["start_ms"], int) or not isinstance(cached["end_ms"], int) or cached["start_ms"] >= cached["end_ms"]:
                    raise DownloadError("Invalid funding coverage bounds")
                if cached["start_ms"] <= start and cached["end_ms"] >= end:
                    return path, events
        except (KeyError, ValueError, TypeError, DownloadError):
            good = False
    if offline:
        raise DownloadError(f"Missing/unverified funding coverage: {coin} [{start},{end})")
    log(f"Downloading settlement ledger and marks: {coin}")
    # Extend a valid cache without dropping previously downloaded research coverage.
    if good:
        start = min(start, cached.get("start_ms", start))
        end = max(end, cached.get("end_ms", end))
    archived = archive_funding(coin, start, end)
    events, cursor = [], start
    while cursor < end:
        query = urllib.parse.urlencode(dict(symbol=coin+"USDT", startTime=cursor, endTime=end-1, limit=1000))
        page = json.loads(http_get(f"{API['futures']}/fundingRate?{query}"))
        write_atomic(CACHE_DIR / "funding_responses" / f"{coin}-{cursor}-{end}.json", json.dumps(page))
        if not isinstance(page, list):
            raise DownloadError(f"Invalid funding response: {coin}")
        if not page:
            break
        for item in page:
            ts = normalize_ms(item["fundingTime"])
            if not cursor <= ts < end:
                continue
            mark = item.get("markPrice")
            if not mark or float(mark) <= 0:
                query_mark = urllib.parse.urlencode(dict(symbol=coin+"USDT", interval="1m", startTime=ts, endTime=ts, limit=1))
                marks = json.loads(http_get(f"{API['futures']}/markPriceKlines?{query_mark}"))
                if not marks or normalize_ms(marks[0][0]) != ts:
                    raise DownloadError(f"Missing settlement mark: {coin} {ts}")
                mark = marks[0][1]
            events.append({"timestamp_ms":ts, "rate":float(item["fundingRate"]), "mark_price":float(mark)})
        following = normalize_ms(page[-1]["fundingTime"]) + 1
        if following <= cursor:
            raise DownloadError("Funding pagination made no progress")
        cursor = following
    # The archive uses nominal boundaries; REST preserves millisecond offsets.
    by_second = {e["timestamp_ms"] // 1000: e for e in events}
    for second, published in archived.items():
        if second not in by_second:
            ts = published["timestamp_ms"]
            query = urllib.parse.urlencode(dict(symbol=coin+"USDT", interval="1m", startTime=ts, endTime=ts, limit=1))
            marks = json.loads(http_get(f"{API['futures']}/markPriceKlines?{query}"))
            if not marks or normalize_ms(marks[0][0]) != ts:
                raise DownloadError(f"Cannot repair archived settlement mark: {coin} {ts}")
            by_second[second] = dict(timestamp_ms=ts, rate=published["rate"], mark_price=float(marks[0][1]))
        event = by_second[second]
        # Published archive rates have eight decimal places.
        if abs(event["rate"] - published["rate"]) > 5.1e-9:
            raise DownloadError(f"Funding archive/REST rate conflict: {coin} {second}")
        event["interval_hours"] = published["interval_hours"]
    events = sorted(by_second.values(), key=lambda e:e["timestamp_ms"])
    validate_funding(events)
    # Binance crypto perpetual intervals are at most eight hours. A larger hole
    # needs source investigation, not an inferred zero rate.
    times = [start//1000] + [e["timestamp_ms"]//1000 for e in events] + [end//1000]
    if any(b-a > 8*3600 for a,b in zip(times,times[1:])):
        raise DownloadError(f"Funding source has an unresolved interval over eight hours: {coin}")
    encoded = json.dumps(events, sort_keys=True, separators=(",",":"), allow_nan=False).encode()
    payload = dict(schema=2, start_ms=start, end_ms=end, events=events, events_sha256=hashlib.sha256(encoded).hexdigest(),
                   source=[ARCHIVE+"/futures/um/monthly/fundingRate", API["futures"]+"/fundingRate"])
    write_atomic(path, json.dumps(payload, separators=(",",":"), allow_nan=False))
    return path, events

def validate_funding(events):
    previous = -1
    for e in events:
        ts, rate, mark = e["timestamp_ms"], e["rate"], e["mark_price"]
        if ts <= previous or not isinstance(ts,int) or not math.isfinite(rate) or not math.isfinite(mark) or mark <= 0:
            raise DownloadError(f"Invalid or duplicate funding settlement: time={ts}, previous={previous}, rate={rate}, mark={mark}")
        previous = ts

def ensure_strategy(config, name, offline=False, warmup_bars=0):
    cfg = json.loads(Path(config).read_text())
    s, run = cfg["strategies"][name], cfg.get("run", {})
    coins = cfg["coins"][:1] if name == "2EMA_crossover" else cfg["coins"]
    if not coins or len(coins) > 16 or len(set(coins)) != len(coins) or any(not c.isalnum() or c.upper()!=c for c in coins):
        raise DownloadError("Expected 1..16 unique uppercase symbols")
    data = (ROOT / cfg.get("data_dir", "data/research")).resolve()
    fixtures = (ROOT/"data/fixtures").resolve()
    if (data == fixtures or fixtures in data.parents) and not offline:
        raise DownloadError("Research downloads cannot modify fixtures")
    tf, market = s["timeframe"], s["market"]
    if market not in API or tf not in TIMEFRAME_MS or warmup_bars < 0:
        raise DownloadError("Invalid market, timeframe, or warmup")
    step = TIMEFRAME_MS[tf]
    explicit_start = bool(run.get("start"))
    start = utc_ms(run["start"]) if explicit_start else utc_ms(cfg.get("history_start", "2017-08"))
    end = utc_ms(run["end"]) if run.get("end") else int(time.time()*1000)//step*step
    if start % step or end % step or start >= end or end > int(time.time()*1000)//step*step:
        raise DownloadError("Run bounds must be ordered completed UTC candle boundaries")
    history = start - warmup_bars*step if explicit_start else start
    if market == "futures":
        history = max(history, utc_ms("2019-09"))
    signals, rows, identities = [], [], {}
    for coin in coins:
        path, series = ensure_candles(data,market,coin,tf,history,end,offline,not explicit_start)
        signals.append(str(path)); rows.append(series)
        identities[str(path)] = digest(path)
    history = max(r[0][0] for r in rows)
    first_trade = start if explicit_start else history + warmup_bars*step
    if first_trade >= end:
        raise DownloadError("Insufficient common history after indicator warmup")
    funding, events = [], []
    execution_tf = tf
    if market == "futures":
        for coin in coins:
            path, series = ensure_funding(data,coin,first_trade,end,offline)
            funding.append(str(path)); events.extend(series); identities[str(path)] = digest(path)
        supported = sorted(((v,k) for k,v in TIMEFRAME_MS.items() if step%v == 0), reverse=True)
        # Retain the reported millisecond offset; the executor applies offset
        # settlements after opening orders. Whole-second offsets still require
        # a compatible finer candle and cannot silently be rounded away.
        execution_tf = next((k for v,k in supported if all((e["timestamp_ms"]//1000)%(v//1000) == 0 for e in events)), "")
        if not execution_tf:
            raise DownloadError("Funding times cannot be represented by supported execution candles")
    execution = signals
    if execution_tf != tf:
        execution = []
        for coin in coins:
            path,_ = ensure_candles(data,market,coin,execution_tf,first_trade,end,offline)
            execution.append(str(path)); identities[str(path)] = digest(path)
    return dict(signal_files=signals, execution_files=execution, funding_files=funding, coins=coins,
                history_start_ms=history, start_ms=first_trade, end_ms=end,
                signal_seconds=step//1000, execution_seconds=TIMEFRAME_MS[execution_tf]//1000,
                data_sha256=identities, sources=[ARCHIVE, API[market]], automatic_repair=not offline)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--ensure", help="strategy name; otherwise prepare every configured strategy")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--warmup-bars", type=int, default=0)
    parser.add_argument("--check", action="store_true", help="validate locally without downloading")
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    names = [args.ensure] if args.ensure else list(cfg["strategies"])
    # Advisory lock protects atomic read/merge/write across simultaneous binaries.
    # The runtime is Linux (native or Docker); no background download service.
    import fcntl
    data = (ROOT/cfg.get("data_dir","data/research")).resolve()
    CACHE_DIR.mkdir(parents=True,exist_ok=True)
    lock_path = CACHE_DIR / ("download-" + hashlib.sha256(str(data).encode()).hexdigest() + ".lock")
    failures = 0
    with lock_path.open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        for name in names:
            try:
                result = ensure_strategy(args.config,name,args.offline or args.check,args.warmup_bars)
                if args.manifest:
                    write_atomic(args.manifest,json.dumps(result,allow_nan=False))
                log(f"{name}: complete {result['start_ms']} .. {result['end_ms']}")
            except (DownloadError, ValueError, KeyError, OSError) as error:
                failures += 1
                log(f"{name}: FAILED: {error}")
    return bool(failures)

if __name__ == "__main__":
    sys.exit(main())
