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

End-to-end benchmark for the simple 2-EMA crossover on `BTC-USDT 1h`:

| Implementation | Parameter sweep | Total time | Relative to C++ |
| --- | ---: | ---: | ---: |
| C++ | 174,936 pairs | 12.99 s | 1.00x |
| Python + Numba | 174,936 pairs | 39.80 s | 3.06x slower |
| Pure Python | 174,936 pairs | ~6323 s | 486x slower (estimated) |

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

On Windows, do this inside WSL. `install.sh` performs the same steps.

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

## Strategies

Whether a strategy runs depends on the bundled data, which is incomplete (see
[Data](#data)).

| Strategy | Market | Timeframe | Indicators | Runs on bundled data |
| --- | --- | --- | --- | :---: |
| `backtest_double_EMA_float` | spot | 1h | 2 EMAs | yes |
| `backtest_double_EMA_StochRSI_float_muti_pair` | spot | 1h | 2 EMAs + StochRSI | yes |
| `BigWill` | spot | 1h | Williams %R, EMAs, Awesome Oscillator | yes |
| `backtest_TRIX_multi_pair` | spot | 1h | TRIX, EMA, StochRSI | yes |
| `BBTREND` | spot | 2h | Bollinger breakout + EMA | yes |
| `SuperTrend_EMA_ATR` | spot | 4h | SuperTrend, EMA, ATR | yes |
| `3EMA_SRSI_ATR` | spot | 5m | 3 EMAs, StochRSI, ATR | no — no 5m spot data |
| `SuperReversal_mtf` | spot | 5m/1h | SuperTrend + EMAs, multi-timeframe | no — no 5m spot data |
| `F_BigWill` | futures | 1h | as `BigWill`, long + short | no — no MATIC/XLM futures |
| `F_BBTREND` | futures | 2h | as `BBTREND`, long + short | no — no 2h futures data |
| `F_SuperReversal_mtf` | futures | 5m/1h | as `SuperReversal_mtf` | no — no MATIC/XLM futures |
| `F_3EMA_SRSI_ATR` | futures | 5m | as `3EMA_SRSI_ATR`, long + short | no — no MATIC/XLM futures |

`STRATEGY_TEMPLATE.cpp` is a complete, working RSI/EMA example to copy when adding your
own — it runs against the bundled data in about two seconds.

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

- `tools.*` — data loading, timestamp alignment, funding lookup, ranges, reporting
- `indicators.*` — the indicator library and TA-Lib integration
- `indicator_cache.hh` — per-pair keyed store of precomputed series
- `trade_core.*` — position open/close, fees, funding, wallet snapshots, drawdown, result
  metrics
- `strategy_runner.hh` — the parameter-sweep loop, banner and result printing, score-file
  cadence
- `Klinef.hh` — the OHLCV container plus its `IndicatorCache`

Strategy files own only their signal logic, their indicator precomputation, and their
parameter ranges.

## Testing

```bash
./run_regression.sh              # build + unit suite + the three numeric fixtures
./run_regression.sh --rebaseline # regenerate fixtures after an intended result change
```

On Windows, `powershell -File .\run_regression.ps1` wraps the same script via WSL.

The gate runs the unit suite (`tests.cpp`, 230+ checks covering trade-core fee maths,
result metrics, drawdown, ranges, the calendar helpers, and every indicator), then
compares three harnesses against `regression_expected/`:

- `verification_regression` — CSV/JSON loading, TA-Lib wrappers, alignment, funding
- `strategy_regression` — one fixed multi-pair spot case through the shared trade core
- `backtest_double_EMA_float` — a full 175k-point sweep, pinned to its winning parameters

CI runs all of that plus ASan/UBSan and TSan builds of the unit suite.

A fixture change means a real result change: explain it in the commit message.

## Data

Provided in `data/` for testing (and outdated). **Coverage is incomplete:**

- Spot: 15m, 1h, 2h, 4h — **no 5m**
- Futures: 1h, 4h, 5m — **no 2h**, and **no MATIC or XLM** at any timeframe

That is why six of the twelve strategies cannot run as configured. Fresh data goes in
`data/data/binance/<timeframe>/<COIN>-USDT.csv` (spot) or
`data/data/futures/<COIN>_USDT-<timeframe>-futures.json` (futures); the archived
`data_downloader_freqtrade` project downloads it via Freqtrade and Docker.

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
  buy<: 30.000000 ; sell>: 60.000000
Max Open Trades  : 2
Gain             : 98.6643%
Porfolio         : 1986.64$ (started with 1000$)
Win rate         : 75.5556%
max DD           : -11.1129%
Gain/DDC         : 7.89175
Score            : 596.265
Calmar ratio     : 0.826031
Number of trades : 180
--------------------------------------------------------------------------
```
