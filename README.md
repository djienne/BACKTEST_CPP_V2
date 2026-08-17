# Trading Strategy Backtester C++

A high-performance parameter-sweep backtester for trading strategies (EMA crossover,
SuperTrend, Bollinger breakout, TRIX, and others), built on TA-Lib.

**Notes:**
- Strategies with an `F_` prefix target futures markets; the rest are spot.
- Price and indicator arrays are `float` (32-bit), which roughly halves numeric array
  memory versus `double`. Wallet balances, fees and result metrics are `double` — those
  accumulate across thousands of operations, and they are what the sweep ranks on.
- All timestamps are interpreted as **UTC**, so results do not depend on the machine's
  timezone.

## Performance

End-to-end benchmark for the simple 2-EMA crossover on `BTC-USDT 1h`, all three measured
on the same machine:

| Implementation | Parameter sweep | Total time | Relative to C++ |
| --- | ---: | ---: | ---: |
| C++ | 174,936 pairs | 12.99 s | 1.00x |
| Python + Numba | 174,936 pairs | 39.80 s | 3.06x slower |
| Pure Python | 174,936 pairs | ~6323 s | 486x slower (estimated) |

Absolute times are strongly hardware- and filesystem-dependent — the same C++ binary
takes roughly 21 s under WSL on an NTFS mount. Treat the ratios as the claim, not the
seconds. The Python implementations are in `python/`.

Two deliberate trade-offs affect the C++ figure:

- Wallet and fee arithmetic moved from `float` to `double`, which costs roughly 10–15% on
  this particular strategy (measured before/after on one machine, with overlapping run-to-run
  ranges) and buys a wallet that no longer accumulates float rounding error into the very
  number the sweep ranks on. Winning parameters are unchanged.
- Strategies that recompute an indicator the sweep does not vary now cache it. That runs
  the other way, and much further: `BigWill` went from 33 to 163 backtests/second with
  peak memory falling from 2,194 MB to 224 MB.

## Getting started

### Required packages (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install -y build-essential libjsoncpp-dev nlohmann-json3-dev gcc g++ unzip
```

### Build TA-Lib

```bash
unzip -o talib.zip
cd talib
./configure --prefix="$(pwd)/talib_install/"
make
make install
cd ..
```

Build TA-Lib **serially** — its own `src/tools/gen_code` target races under `make -j` and
will fail intermittently. The backtester itself builds fine in parallel.

On Windows, do this inside WSL. `./install.sh` does this and the build, but not the data
download — that is the next step.

A docker-compose setup is bundled as `docker_easy.zip`; unzip it and see the README
inside.

### Get the data

The repo ships only the small slice the regression fixtures need. Everything else is
downloaded on demand from Binance's public archive:

```bash
python3 tools/download_data.py --check    # what's missing, and which strategies it blocks
python3 tools/download_data.py            # fetch it (stdlib only, no pip install)
```

The downloader reads the same `backtest_config.json` the strategies do, so it cannot
fetch a different universe than they open. Re-runs are incremental — completed months are
cached under `.data_cache/`. See [Data](#data) for formats and caveats.

### Build and run

```bash
make -j all              # all strategies + regression drivers + unit tests
./backtest_double_EMA_float.exe
```

Binaries carry an rpath to the in-tree TA-Lib, so no `LD_LIBRARY_PATH` is needed.

| Command | Purpose |
| --- | --- |
| `make` | default strategy only |
| `make <Name>` | one strategy, e.g. `make BigWill` |
| `make tests && ./tests.exe` | unit suite (fast, no data dependency) |
| `make BUILD=debug all` | `-O0 -g` |
| `make BUILD=asan tests` | AddressSanitizer + UBSan |
| `make BUILD=tsan all` | ThreadSanitizer |
| `make format` | clang-format in place |
| `make help` | list discovered strategies |

Strategies are discovered from `*.cpp` automatically — adding one needs no Makefile edit.

## Configuration

`backtest_config.json` is the single source of truth for the traded universe. Both the
strategies and the downloader read it, so the data that gets fetched cannot drift from the
data the strategies open:

```jsonc
{
  "data_dir": "./data/data",
  "coins": ["BTC", "ETH"],
  "history_start": "2017-08",
  "strategies": {
    "BBTREND":           { "market": "spot",    "timeframe": "2h" },
    "F_SuperReversal_mtf": { "market": "futures", "timeframe": "5m", "htf": "1h" }
  }
}
```

Keys are each binary's `STRAT_NAME` constant, **not** its file name. A strategy with no
entry aborts at startup rather than guessing a universe nobody downloaded. Override the
file location with `BACKTEST_CONFIG=/path/to/config.json`.

Changing `coins` changes every strategy at once. The count is a runtime value;
`backtest_config::MAX_PAIRS` (16) is only the compile-time ceiling that keeps the
per-pair arrays fixed-size, so raising the universe past 16 means bumping that constant
and rebuilding.

**Symbols must exist on both Binance spot and USDT-M futures.** Three that earlier
versions of this repo referenced no longer do, so no data can be downloaded for them:

| Symbol | Status |
| --- | --- |
| MATIC | Swapped to POL 2024-09; USDT-M perpetual settled and delisted. POL has no USDT-M perpetual either. |
| XMR | Delisted from Binance entirely 2024-02-20 — spot, margin and futures. |
| EOS | No `EOSUSDT` USDT-M perpetual. |

## Strategies

All of them trade the configured `coins` on the timeframe their config entry names.

| Strategy | Market | Timeframe | Indicators |
| --- | --- | --- | --- |
| `backtest_double_EMA_float` | spot | 1h | 2 EMAs (single asset) |
| `backtest_double_EMA_StochRSI_float_muti_pair` | spot | 1h | 2 EMAs + StochRSI |
| `BigWill` | spot | 1h | Williams %R, EMAs, Awesome Oscillator |
| `backtest_TRIX_multi_pair` | spot | 1h | TRIX, EMA, StochRSI |
| `BBTREND` | spot | 2h | Bollinger breakout + EMA |
| `SuperTrend_EMA_ATR` | spot | 4h | SuperTrend, EMA, ATR |
| `3EMA_SRSI_ATR` | spot | 5m | 3 EMAs, StochRSI, ATR |
| `SuperReversal_mtf` | spot | 5m/1h | SuperTrend + EMAs, multi-timeframe |
| `F_BigWill` | futures | 1h | as `BigWill`, long + short |
| `F_BBTREND` | futures | 2h | as `BBTREND`, long + short |
| `F_SuperReversal_mtf` | futures | 5m/1h | as `SuperReversal_mtf` |
| `F_3EMA_SRSI_ATR` | futures | 5m | as `3EMA_SRSI_ATR`, long + short |

`python3 tools/download_data.py --check` reports which of them currently have their data.

`STRATEGY_TEMPLATE.cpp` is a complete, working RSI/EMA example to copy when adding your
own — it runs in about two seconds.

### Futures model

Leverage is 1 and **there is no liquidation model**. A short is valued as its mirrored
long, `2 * entry - price`, which is exact at leverage 1 but has no floor: a position that
runs far enough against the book can drive the wallet negative and keep trading. Read
futures results as an un-liquidated upper bound.

## Indicators

All in `indicators.hh` / `indicators.cpp`. Output length always equals input length, with
the warmup region left-padded with zeros and its length reported — see
[CONTRIBUTING.md](CONTRIBUTING.md).

| Group | Indicators |
| --- | --- |
| Moving averages | SMA, EMA, WMA, DEMA, TEMA, KAMA, HMA |
| Momentum | RSI, MACD, Stochastic (%K/%D), StochRSI (K/D/raw), CCI, ROC, MOM, Williams %R, Aroon, Ultimate Oscillator, Awesome Oscillator |
| Trend strength | ADX, +DI, −DI, DMI, Parabolic SAR, TRIX, SuperTrend |
| Volatility | ATR, NATR, True Range, StdDev, Bollinger Bands, Keltner Channels, Donchian Channels |
| Volume | OBV, MFI, A/D, Chaikin A/D Oscillator, rolling VWAP, relative volume |
| Transforms | HL2, HLC3, OHLC4, Heikin-Ashi |
| Multi-timeframe | `RESAMPLE_TIMEFRAME`, `PROJECT_HTF_TO_LTF` |
| Rolling extremes | MIN, MAX |

**StochRSI caveat:** `TALIB_STOCHRSI_*` round the normalized RSI to 3 decimals before
smoothing. This is inherited behaviour, kept because tracked results and tuned parameter
sets depend on it, but it means they do not match TA-Lib's own `STOCHRSI` exactly.

## Architecture

- `backtest_config.json` — traded universe and per-strategy market/timeframe
- `config.*` — loads that file; shared with the downloader so the two cannot drift
- `tools.*` — data loading, timestamp alignment, funding lookup, ranges, reporting
- `indicators.*` — the indicator library and TA-Lib integration
- `indicator_cache.hh` — per-pair keyed store of precomputed series
- `trade_core.*` — position open/close, fees, funding, wallet snapshots, drawdown, result
  metrics
- `strategy_runner.hh` — the parameter-sweep loop, banner and result printing, score-file
  cadence
- `Klinef.hh` — the OHLCV container plus its `IndicatorCache`
- `tools/download_data.py` — fetches market data; `tools/smoke_all_strategies.sh` runs
  every strategy end to end

Strategy files own only their signal logic, their indicator precomputation, and their
parameter ranges.

### Exit fills are modelled at the bar close

Every strategy books its exit at `close[ii]`, including exits triggered by an intrabar
event. A stop-loss reached by the bar's low but closing higher is therefore filled at
that higher close rather than at the stop, which flatters stop-outs on bars that recover.
The same applies to take-profits. This is a repo-wide simplification worth knowing when
reading results; making fills realistic would change every strategy's numbers.

## Testing

```bash
./run_regression.sh              # build + unit suite + the three numeric fixtures
./run_regression.sh --rebaseline # regenerate fixtures after an intended result change
```

On Windows, `powershell -File .\run_regression.ps1` wraps the same script via WSL.

The gate runs two suites and three numeric harnesses:

- **`tests.exe`** — 279 checks over trade-core fee maths, result metrics, drawdown,
  ranges, the UTC calendar helpers, every indicator, and regressions for the specific
  bugs found in review (see the comments on each test for what it pins).
- **`tools/test_download_data.py`** — 37 checks on the downloader's pure logic. No
  network, so it is safe to gate every push on.
- `verification_regression` — CSV/JSON loading, TA-Lib wrappers, alignment, funding
- `strategy_regression` — one fixed multi-pair spot case through the shared trade core
- `backtest_double_EMA_float` — a full 175k-point sweep, pinned to its winning parameters

The three harnesses compare against `regression_expected/`, and run against the small
committed data slice, so the gate needs no network.

`tools/smoke_all_strategies.sh [seconds]` is the end-to-end companion: it launches every
strategy in turn and reports whether each loaded its data, swept, and traded. It is not
part of the gate — it needs the downloaded dataset and takes minutes.

CI runs the gate plus ASan/UBSan and TSan builds of the unit suite.

A fixture change means a real result change: explain it in the commit message.

## Data

`data/` holds only the slice the regression fixtures are computed from — four files,
about 7.5 MB. Everything else is downloaded:

```bash
python3 tools/download_data.py --check                       # coverage report
python3 tools/download_data.py                               # fetch what's missing
python3 tools/download_data.py --coins BTC --timeframes 5m   # a subset
python3 tools/download_data.py --start 2021-01 --jobs 8      # narrower history, faster
```

It pulls Binance's public monthly archives (`data.binance.vision`) rather than the REST
API, so there is no pagination or rate limiting, and it verifies the SHA256 checksum that
ships beside each file. Downloads are cached under `.data_cache/`; both that and the
downloaded datasets are gitignored.

Two things worth knowing:

- **Binance restates historical spot klines.** Re-downloading BTC 1h reproduces 50,960 of
  50,976 rows exactly, but 16 differ — 5 of them in `close`. Futures klines and funding
  rates matched exactly. The four files the fixtures depend on are therefore pinned:
  `--force` skips them unless you also pass `--force-fixtures`, which will move
  `regression_expected/` and needs a deliberate rebaseline.
- **The archive timestamp unit changed.** Files from 2025-01 onward use microseconds
  where earlier ones use milliseconds; the downloader normalizes both to milliseconds.

Fresh data goes in `data/data/binance/<timeframe>/<COIN>-USDT.csv` for spot, or
`data/data/futures/<COIN>_USDT-<timeframe>-futures.json` for futures.

Spot CSVs have a header row and six columns — open time in **milliseconds** (the column
is labelled `date`), then open, high, low, close, volume:

```
date,open,high,low,close,volume
1502942400000,4261.48,4313.62,4261.32,4308.83,47.181009
```

Futures JSON is an array of `[open_time_ms, open, high, low, close, volume]`, and funding
rates are `[timestamp_ms, rate]` in `<COIN>_USDT-8h-funding_rate.json`. Both are what
Freqtrade's `download-data` produces. No downloader is bundled with this repo.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for how to add an indicator or a strategy.

## Example output

```
--------------------------------------------------------------------------
BEST PARAMETER SET FOUND:
--------------------------------------------------------------------------
Strategy         : TEMPLATE_RSI_EMA
Parameters       :
  RSI: 7 ; EMA: 50
  buy<: 35.000000 ; sell>: 60.000000
Max Open Trades  : 2
Gain             : 109.336%
Porfolio         : 2093.36$ (started with 1000$)
Win rate         : 72.9075%
max DD           : -17.1226%
Gain/DDC         : 5.29213
Score            : 385.836
Calmar ratio     : 0.476041
Number of trades : 454
--------------------------------------------------------------------------
```

(`STRATEGY_TEMPLATE.exe`, BTC + ETH on 1h, 1200 backtests in under a second.)
