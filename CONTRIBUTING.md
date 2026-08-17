# Contributing

Two walkthroughs: adding an indicator, and adding a strategy. Both end with "run the
regression gate", which is the only sign-off that matters here.

## Build and test loop

```bash
unzip -o talib.zip
cd talib && ./configure --prefix="$(pwd)/talib_install/" && make && make install && cd ..
```

Build TA-Lib **serially** — its own `src/tools/gen_code` target races under `make -j`.

Then:

```bash
make -j all          # every strategy, both regression drivers, the unit suite
./run_regression.sh  # build + both unit suites + the three numeric fixtures
```

The gate runs against the small committed data slice and needs no network. To actually
run a strategy you need the rest of the data:

```bash
python3 tools/download_data.py --check   # what's missing and which strategies it blocks
python3 tools/download_data.py           # fetch it
./tools/smoke_all_strategies.sh 60       # launch every strategy, check each one trades
```

Other build modes:

| Command | What it gives you |
| --- | --- |
| `make BUILD=debug all` | `-O0 -g` |
| `make BUILD=asan tests` | AddressSanitizer + UBSan |
| `make BUILD=tsan all` | ThreadSanitizer (the `F_*` strategies are multithreaded) |
| `make format` | clang-format, in place |
| `make help` | lists discovered strategies |

Objects live under `build/$(BUILD)/`, so the modes never clobber each other, and header
dependencies are tracked — editing a header rebuilds exactly what depends on it.

## Adding an indicator

Everything lives in `indicators.hh` / `indicators.cpp`.

**If TA-Lib already implements it**, it is a one-line wrapper. Pick the adapter matching
the input shape — `run_single_input`, `run_ohlc`, `run_ohlc_no_period`, `run_ohlcv` — and
forward to the `TA_S_*` function:

```cpp
std::vector<float> TALIB_CCI(const std::vector<float> &high, const std::vector<float> &low,
                             const std::vector<float> &close, const int period, size_t *warmup)
{
    return run_ohlc(high, low, close, period, TA_S_CCI, "TALIB_CCI", warmup);
}
```

Multi-output indicators return a small struct carrying a `warmup` field — see
`TALIB_MACD` or `TALIB_STOCH`.

**If TA-Lib does not implement it**, compose it from the wrappers and keep the same
conventions — `KELTNER_CHANNELS`, `DONCHIAN_CHANNELS`, `VWAP_ROLLING` and `TALIB_HMA` are
worked examples.

Three conventions every indicator follows:

1. **Output length equals input length.** TA-Lib returns a compacted array; the adapters
   scatter it back and left-pad the warmup with zeros so `series[i]` lines up with candle
   `i`.
2. **Warmup is reported, never hidden.** Single-output wrappers take an optional
   `size_t *warmup`; multi-output structs carry a `warmup` field. `series[i]` is a real
   value only for `i >= warmup` — reading earlier compares against zero padding, which is
   how a strategy silently trades on nothing. Derive the loop start with
   `strategy_runner::first_tradable_index({warm_a, warm_b}, start_indexes[0], lookback)`
   rather than a hard-coded constant, so raising an indicator period cannot quietly
   outgrow it.
3. **Empty input gives empty output**, rather than aborting.

Then add a test in `tests.cpp`. Assert output length, warmup position, **and at least one
value with a known analytic answer** — a shape-only test passes on a wrapper plumbed to
the wrong TA-Lib function. Existing examples: SMA and WMA on a linear ramp (closed form),
OBV on a hand-checked up/down sequence, MACD's `histogram == macd - signal` identity,
bounded `[0, 100]` ranges for the oscillators.

## Adding a strategy

Copy `STRATEGY_TEMPLATE.cpp`. `make` discovers `*.cpp` in the repo root automatically, so
there is no list to update — `make MyStrategy` and `make MyStrategy.exe` both work
immediately, and `run_regression.sh` compiles it because it runs `make all`.

**Add a `backtest_config.json` entry first**, keyed by your `STRAT_NAME` constant (not the
file name). Without one the binary aborts at startup, deliberately: a strategy that
guesses its own market and timeframe is how this repo ended up shipping binaries that
read files nobody downloads. Then:

```bash
python3 tools/download_data.py --check   # confirms your entry's data exists
python3 tools/download_data.py           # fetches it if not
```

### Pair count is a runtime value

`backtest_config::MAX_PAIRS` (16) is a compile-time ceiling used only to size the
fixed-length per-pair arrays — `PortfolioState<MAX_PAIRS>`, `std::array<float,
MAX_PAIRS>`. The number actually traded is the runtime `NB_PAIRS`, filled from the config
in `main()`. So:

- size arrays and template arguments with `MAX_PAIRS`
- bound every loop with `NB_PAIRS`
- construct the portfolio as `PortfolioState<MAX_PAIRS> portfolio(initial, NB_PAIRS)`

They are fixed-size rather than `std::vector` on purpose: `PROCESS` runs millions of times
per sweep and a heap allocation there would dominate.

Anything derived from the pair count must be computed **after** the config load, not at
declaration time — filling a `max_open_trades` sweep before `NB_PAIRS` was set left it
empty and made `F_BigWill` run zero backtests while still reporting success.

The shape of a strategy is three functions:

- **`CALCULATE_INDICATORS`** — compute everything the sweep does *not* vary, once, into
  each pair's `IndicatorCache`.
- **`PROCESS`** — one backtest for one parameter set. Hoist your indicator references out
  of the bar loop (see below), then walk bars and call the `trade_core` helpers.
- **`main`** — load data, build the parameter list, hand it to `strategy_runner::sweep`.

### Indicator storage

Each `KLINEf` owns an `IndicatorCache`, keyed by name plus parameters:

```cpp
PAIRS[ic].indicators.put(IndicatorCache::key("EMA", period), TALIB_EMA(PAIRS[ic].close, period));
```

`get()` aborts on a missing key rather than returning an empty series, so a typo fails
loudly instead of making every comparison quietly false.

**Hoist lookups out of the bar loop.** `get()` hashes a string:

```cpp
// once per backtest
const std::vector<float> &ema = PAIRS[ic].indicators.get(IndicatorCache::key("EMA", n));
// then in the loop
if (ema[ii] > ...)
```

not `PAIRS[ic].indicators.get(...)[ii]` inside the loop.

### Cache anything the sweep re-derives

If an indicator depends on only *some* of the swept parameters, compute it lazily in
`PROCESS` and keep it in the cache:

```cpp
const std::string ao_key = IndicatorCache::key("AO", fast, slow);
if (!PAIRS[ic].indicators.has(ao_key))
{
    PAIRS[ic].indicators.put(ao_key, TALIB_AO(PAIRS[ic].high, PAIRS[ic].low, fast, slow));
}
```

This is not a micro-optimisation. `BigWill` recomputed AO on every parameter set even
though AO depends only on `(fast, slow)` — about 1,050 distinct values out of five million
combinations.

### Trade mechanics

Use `trade_core` rather than hand-rolling fills: `open_spot_long` / `close_spot_long`,
`open_futures_long` / `close_futures_long` / `open_futures_short` / `close_futures_short`,
`apply_funding_fee`, and `record_spot_snapshot` / `record_futures_snapshot` for the equity
curve. Close positions **before** opening new ones on each bar.

Futures caveat: leverage is 1 and there is no liquidation model. A short is valued as its
mirrored long, so a position that runs far enough against the book can drive the wallet
negative and keep trading. Read futures results as an un-liquidated upper bound.

### Multi-timeframe

Use `RESAMPLE_TIMEFRAME` and `PROJECT_HTF_TO_LTF`. Do not hand-roll the slicing: the
resampler finds the first true higher-timeframe boundary (the bundled 5m data starts at
16:55, so slicing from index 0 builds candles spanning 16:55 to 17:50) and the projection
applies the shift that keeps a higher-timeframe value invisible until its candle closes.

## Before opening a PR

```bash
./run_regression.sh
make BUILD=asan tests && ./tests.exe
```

If a numeric fixture changes, that is a real result change — say why in the commit
message, and regenerate deliberately with `./run_regression.sh --rebaseline`. A fixture
edit with no explanation is the one thing that will get a change sent back.

## Data

Only the slice the regression fixtures are computed from is committed (four files, ~7.5
MB). Everything else comes from `tools/download_data.py`, which reads the same
`backtest_config.json` the strategies do.

```bash
python3 tools/download_data.py --check                    # coverage + what it blocks
python3 tools/download_data.py                            # fetch missing series
python3 -m unittest discover -s tools -p 'test_*.py'      # downloader tests (no network)
```

Two traps worth knowing about, both covered by tests:

- **Binance restates historical spot klines.** Re-downloading BTC 1h reproduces all but
  16 of 50,976 rows, and 5 of those differ in `close`. The four fixture files are pinned;
  `--force` skips them unless `--force-fixtures` is given, because rebuilding them moves
  `regression_expected/`.
- **Archive timestamps switched from milliseconds to microseconds in 2025-01.** Mixing
  the two sorts a series wrongly and produced dates in the year 58595; `normalize_ms()`
  handles it.

Symbols must exist on both Binance spot and USDT-M futures. MATIC (now POL, which has no
perpetual), XMR and EOS do not, so no data can be fetched for them.
