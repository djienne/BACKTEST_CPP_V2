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
./run_regression.sh  # build + unit tests + the three numeric fixtures
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
   how a strategy silently trades on nothing.
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

The bundled dataset does not cover every strategy. Spot has 15m/1h/2h/4h (no 5m); futures
has 1h/4h/5m (no 2h) and lacks MATIC and XLM entirely. Six of the twelve strategies
therefore cannot run as configured — see the table in `README.md`. Fresh data goes in
`data/data/binance/<timeframe>/<COIN>-USDT.csv` (spot) or
`data/data/futures/<COIN>_USDT-<timeframe>-futures.json` (futures).
