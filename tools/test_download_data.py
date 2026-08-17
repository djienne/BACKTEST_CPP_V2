#!/usr/bin/env python3
"""
Unit tests for the pure parts of download_data.py -- no network access, so CI can run
them on every push.

    python3 -m unittest discover -s tools -v
"""

import io
import json
import sys
import tempfile
import unittest
import zipfile
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import download_data as dd  # noqa: E402


class TestNormalizeMs(unittest.TestCase):
    """
    Binance switched the archive timestamp unit from milliseconds to microseconds with
    the 2025-01 files. Getting this wrong put dates in the year 58595 and produced 14,261
    spurious gaps in a single series, so it is pinned hard.
    """

    def test_milliseconds_pass_through(self):
        # 2017-08-17 04:00 UTC, the first row of the tracked BTC 1h CSV.
        self.assertEqual(dd.normalize_ms(1502942400000), 1502942400000)
        # 2026-08-16, near the end of the archive.
        self.assertEqual(dd.normalize_ms(1755300000000), 1755300000000)

    def test_microseconds_converted(self):
        # 2025-01-01 00:00 UTC as microseconds -> the same instant in milliseconds.
        self.assertEqual(dd.normalize_ms(1735689600000000), 1735689600000)
        self.assertEqual(dd.normalize_ms(1748736000000000), 1748736000000)

    def test_both_units_agree_on_the_same_instant(self):
        instant = datetime(2025, 6, 1, tzinfo=timezone.utc)
        ms = int(instant.timestamp() * 1000)
        us = ms * 1000
        self.assertEqual(dd.normalize_ms(ms), dd.normalize_ms(us))

    def test_threshold_has_headroom_either_side(self):
        # Year 2100 in ms is still far below the threshold...
        year_2100_ms = 4102444800000
        self.assertEqual(dd.normalize_ms(year_2100_ms), year_2100_ms)
        # ...and year 2000 in us is still far above it.
        year_2000_us = 946684800000000
        self.assertEqual(dd.normalize_ms(year_2000_us), 946684800000)


class TestTrimNumber(unittest.TestCase):
    """Output must match the trimmed form the tracked CSVs are stored in."""

    def test_trailing_zeros_removed(self):
        self.assertEqual(dd.trim_number("4261.48000000"), "4261.48")
        self.assertEqual(dd.trim_number("47.18100900"), "47.181009")

    def test_all_zero_decimals_collapse_to_int(self):
        self.assertEqual(dd.trim_number("0.00000000"), "0")
        self.assertEqual(dd.trim_number("10000.00000000"), "10000")

    def test_integers_untouched(self):
        self.assertEqual(dd.trim_number("1502942400000"), "1502942400000")

    def test_no_precision_lost(self):
        # Done on the string, so a value float() would round must survive intact.
        self.assertEqual(dd.trim_number("0.000000012345678901"), "0.000000012345678901")

    def test_negative(self):
        self.assertEqual(dd.trim_number("-0.00003766"), "-0.00003766")


class TestMonthRange(unittest.TestCase):
    def test_spans_year_boundary(self):
        months = dd.month_range("2024-11", datetime(2025, 2, 1, tzinfo=timezone.utc))
        self.assertEqual(months, [(2024, 11), (2024, 12), (2025, 1)])

    def test_end_is_exclusive(self):
        months = dd.month_range("2024-01", datetime(2024, 1, 15, tzinfo=timezone.utc))
        self.assertEqual(months, [])

    def test_single_month(self):
        months = dd.month_range("2024-01", datetime(2024, 2, 1, tzinfo=timezone.utc))
        self.assertEqual(months, [(2024, 1)])


class TestDaysInMonth(unittest.TestCase):
    def test_lengths(self):
        self.assertEqual(dd.days_in_month(2024, 1), 31)
        self.assertEqual(dd.days_in_month(2024, 2), 29)  # leap year
        self.assertEqual(dd.days_in_month(2023, 2), 28)
        self.assertEqual(dd.days_in_month(2024, 4), 30)
        self.assertEqual(dd.days_in_month(2024, 12), 31)


def make_zip(rows, header=None):
    buf = io.BytesIO()
    lines = ([",".join(header)] if header else []) + [",".join(r) for r in rows]
    with zipfile.ZipFile(buf, "w") as zf:
        zf.writestr("data.csv", "\n".join(lines) + "\n")
    return buf.getvalue()


class TestRowsFromZip(unittest.TestCase):
    """Newer archives carry a header row; older ones do not. Both must parse."""

    KLINE = [["1502942400000", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"],
             ["1502946000000", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]]

    def test_without_header(self):
        rows = dd.rows_from_zip(make_zip(self.KLINE), 6)
        self.assertEqual(len(rows), 2)

    def test_with_header(self):
        header = ["open_time", "open", "high", "low", "close", "volume",
                  "close_time", "qav", "trades", "tbbav", "tbqav", "ignore"]
        rows = dd.rows_from_zip(make_zip(self.KLINE, header), 6)
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0][0], "1502942400000")

    def test_short_rows_dropped(self):
        rows = dd.rows_from_zip(make_zip([["123", "1"]]), 6)
        self.assertEqual(rows, [])


class TestWriters(unittest.TestCase):
    """Output must be exactly what the C++ loaders in tools.cpp parse."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_spot_csv_shape(self):
        series = dd.Series([["1502942400000", "4261.48", "4313.62", "4261.32", "4308.83", "47.181009"]], [], 0)
        out = self.dir / "BTC-USDT.csv"
        dd.write_spot_csv(out, series)
        lines = out.read_text().splitlines()
        self.assertEqual(lines[0], "date,open,high,low,close,volume")
        self.assertEqual(lines[1], "1502942400000,4261.48,4313.62,4261.32,4308.83,47.181009")

    def test_futures_json_shape(self):
        series = dd.Series([["1567965300000", "10000", "10000", "10000", "10000", "0.002"]], [], 0)
        out = self.dir / "f.json"
        dd.write_futures_json(out, series)
        rows = json.loads(out.read_text())
        self.assertEqual(rows, [[1567965300000, 10000.0, 10000.0, 10000.0, 10000.0, 0.002]])
        self.assertIsInstance(rows[0][0], int)  # timestamp must not become a float

    def test_funding_json_is_six_elements(self):
        # freqtrade's shape: the loader reads [0] and [1], the zeros keep it compatible.
        series = dd.Series([["1660464000000", "-0.00000333"]], [], 0)
        out = self.dir / "fr.json"
        dd.write_funding_json(out, series)
        rows = json.loads(out.read_text())
        self.assertEqual(len(rows[0]), 6)
        self.assertEqual(rows[0][0], 1660464000000)
        self.assertAlmostEqual(rows[0][1], -0.00000333)
        self.assertEqual(rows[0][2:], [0.0, 0.0, 0.0, 0.0])

    def test_write_is_atomic(self):
        out = self.dir / "sub" / "x.csv"
        dd.write_atomic(out, "hello\n")
        self.assertEqual(out.read_text(), "hello\n")
        self.assertFalse(out.with_suffix(".csv.tmp").exists())


class TestGapDetection(unittest.TestCase):
    def test_contiguous_series_has_no_gaps(self):
        step = dd.TIMEFRAME_MS["1h"]
        times = [1502942400000 + i * step for i in range(10)]
        self.assertEqual(dd.count_gaps(times, "1h"), 0)

    def test_missing_candle_counted(self):
        step = dd.TIMEFRAME_MS["1h"]
        times = [1502942400000 + i * step for i in range(10)]
        del times[5]
        self.assertEqual(dd.count_gaps(times, "1h"), 1)

    def test_unknown_timeframe_is_not_an_error(self):
        self.assertEqual(dd.count_gaps([1, 2, 3], "7m"), 0)


class TestJobPaths(unittest.TestCase):
    """Paths must land exactly where the strategies look."""

    def test_output_paths(self):
        data = Path("/data")
        self.assertEqual(dd.Job("spot", "BTC", "1h").output(data),
                         data / "binance/1h/BTC-USDT.csv")
        self.assertEqual(dd.Job("futures", "ETH", "5m").output(data),
                         data / "futures/ETH_USDT-5m-futures.json")
        self.assertEqual(dd.Job("funding", "BTC", "8h").output(data),
                         data / "futures/BTC_USDT-8h-funding_rate.json")

    def test_symbol(self):
        self.assertEqual(dd.Job("spot", "BTC", "1h").symbol, "BTCUSDT")


class TestArchiveUrls(unittest.TestCase):
    def test_spot_monthly(self):
        url = dd.archive_url("spot", "BTCUSDT", "klines", "1h", 2024, 3)
        self.assertTrue(url.endswith("/data/spot/monthly/klines/BTCUSDT/1h/BTCUSDT-1h-2024-03.zip"))

    def test_futures_monthly(self):
        url = dd.archive_url("futures", "ETHUSDT", "klines", "5m", 2024, 12)
        self.assertTrue(url.endswith("/data/futures/um/monthly/klines/ETHUSDT/5m/ETHUSDT-5m-2024-12.zip"))

    def test_funding_monthly(self):
        url = dd.archive_url("futures", "BTCUSDT", "fundingRate", "", 2023, 1)
        self.assertTrue(url.endswith("/data/futures/um/monthly/fundingRate/BTCUSDT/BTCUSDT-fundingRate-2023-01.zip"))

    def test_daily_variant(self):
        url = dd.archive_url("spot", "BTCUSDT", "klines", "1h", 2026, 8, day=5)
        self.assertIn("/daily/", url)
        self.assertTrue(url.endswith("BTCUSDT-1h-2026-08-05.zip"))


class TestJobsFromConfig(unittest.TestCase):
    def test_derives_every_series_a_strategy_needs(self):
        cfg = {
            "coins": ["BTC", "ETH"],
            "strategies": {
                "A": {"market": "spot", "timeframe": "1h"},
                "B": {"market": "futures", "timeframe": "5m", "htf": "1h"},
            },
        }
        jobs = dd.jobs_from_config(cfg)
        got = {(j.market, j.coin, j.timeframe) for j in jobs}
        self.assertIn(("spot", "BTC", "1h"), got)
        self.assertIn(("futures", "ETH", "5m"), got)
        self.assertIn(("futures", "BTC", "1h"), got)   # the htf
        self.assertIn(("funding", "BTC", "8h"), got)   # futures implies funding
        self.assertNotIn(("funding", "BTC", "1h"), got)

    def test_no_duplicates(self):
        cfg = {
            "coins": ["BTC"],
            "strategies": {
                "A": {"market": "spot", "timeframe": "1h"},
                "B": {"market": "spot", "timeframe": "1h"},
            },
        }
        self.assertEqual(len(dd.jobs_from_config(cfg)), 1)


class TestFixtureGuard(unittest.TestCase):
    """
    Binance restates historical spot klines, so re-downloading the files the regression
    fixtures are computed from would silently shift them.
    """

    def test_pinned_files_recognised(self):
        data = Path("/data")
        self.assertTrue(dd.is_fixture_pinned(dd.Job("spot", "BTC", "1h"), data))
        self.assertTrue(dd.is_fixture_pinned(dd.Job("futures", "BTC", "1h"), data))
        self.assertTrue(dd.is_fixture_pinned(dd.Job("funding", "BTC", "8h"), data))

    def test_other_files_not_pinned(self):
        data = Path("/data")
        self.assertFalse(dd.is_fixture_pinned(dd.Job("spot", "BTC", "5m"), data))
        self.assertFalse(dd.is_fixture_pinned(dd.Job("spot", "ETH", "4h"), data))


class TestRealConfig(unittest.TestCase):
    """The shipped config must stay loadable and self-consistent."""

    def setUp(self):
        self.cfg = dd.load_config(dd.DEFAULT_CONFIG)

    def test_has_coins_and_strategies(self):
        self.assertTrue(self.cfg["coins"])
        self.assertTrue(self.cfg["strategies"])

    def test_no_delisted_symbols(self):
        # These have no USDT-M perpetual, so their data cannot be downloaded at all.
        for dead in ("MATIC", "XMR", "EOS", "POL"):
            self.assertNotIn(dead, self.cfg["coins"],
                             f"{dead} is not listed on Binance USDT-M futures")

    def test_every_strategy_has_a_market_and_timeframe(self):
        for name, entry in self.cfg["strategies"].items():
            self.assertIn(entry.get("market"), ("spot", "futures"), name)
            self.assertTrue(entry.get("timeframe"), name)

    def test_jobs_derive_cleanly(self):
        self.assertTrue(dd.jobs_from_config(self.cfg))


if __name__ == "__main__":
    unittest.main(verbosity=2)
