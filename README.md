# Trading Strategy Backtester

C++17 research backtests for Binance spot and USDT-margined perpetual futures.
Thirteen strategy programs share one data loader, indicator cache, search runner,
and OHLC execution model. Results describe this model; they are not live-trading
performance estimates.

## Build and check

Docker Compose is the supported build and verification environment on Windows
and Linux. From the repository directory:

```sh
docker compose build
docker compose run --rm -T backtest bash install.sh
docker compose run --rm -T backtest bash run_regression.sh
docker compose -f compose.yaml -f compose.sanitizers.yaml run --rm -T backtest setarch x86_64 -R bash run_regression.sh asan
docker compose -f compose.yaml -f compose.sanitizers.yaml run --rm -T backtest setarch x86_64 -R bash run_regression.sh tsan
```

On Windows, `./run_regression.ps1 -Mode release` runs the same gate.
The image includes the compiler, Python, NumPy, Numba, JSON headers, and
clang-format. The bundled `talib.zip` is compiled serially by the installer.
Release, debug, ASan/UBSan, and TSan executables have separate directories under
`build/`. The `.exe` suffix is historical: these are Linux executables.

The regression gate is offline and compares hand-computed accounts, analytic
indicator cases, causal properties, and an independent Python EMA ledger.
It does not download market data or overwrite expected outputs. Sanitizer
failures are failures; a host/runtime incompatibility must be reported explicitly.
On x86-64, the sanitizer commands above disable address randomization for each process
only. The optional Compose override permits that operation in the disposable
container; no host-wide setting is changed. On another architecture, use its
machine name in place of x86_64.

## Run a study

Edit `backtest_config.json`, then run the desired binary, for example:

```sh
docker compose run --rm -T backtest build/release/BigWill.exe
```

Every direct strategy run checks data first and tries to download missing or
invalid candles and missing funding coverage automatically. Explicit dates are
recommended for a repeatable study:

```json
"run": {
  "start": "2024-01-01",
  "end": "2024-04-01",
  "max_trials": 1000,
  "seed": 42,
  "workers": 2,
  "holdout_fraction": 0.2
}
```

Dates are UTC and describe a half-open interval: start included, end excluded.
They must align with the configured signal timeframe; the end cannot include
an unfinished candle. Indicator prehistory is fetched before the start.
Empty dates request the configured history through the latest complete candle,
with trading starting after common prehistory and warmup. Genuine exchange
outages in that history can prevent a study from running.

`coins` is an ordered list of 1–16 unique uppercase base symbols. Pair priority
uses that order; the single-asset EMA program uses only the first coin.
`data_dir` is relative to the repository root. Use
`-e BACKTEST_CONFIG=/workspace/my_config.json` before the Compose service name
to select another configuration. Use `-e BACKTEST_OFFLINE=1` to validate cached
coverage and fail on holes without downloading.

| Configuration name | Binary under build/release | Signal timeframe |
|---|---|---|
| 2EMA_crossover | backtest_double_EMA_float.exe | 1h |
| 2EMA_crossover_StochRSI | backtest_double_EMA_StochRSI_float_muti_pair.exe | 1h |
| BigWill | BigWill.exe | 1h |
| TRIX | backtest_TRIX_multi_pair.exe | 1h |
| BBTREND | BBTREND.exe | 2h |
| SuperTrend_EMA_ATR | SuperTrend_EMA_ATR.exe | 4h |
| EMA3_SRSI_ATR | 3EMA_SRSI_ATR.exe | 5m |
| SuperReversal_mtf | SuperReversal_mtf.exe | 5m, higher 1h |
| TEMPLATE_RSI_EMA | STRATEGY_TEMPLATE.exe | 1h |
| F_BigWill | F_BigWill.exe | 1h |
| F_BBTREND | F_BBTREND.exe | 2h |
| F_SuperReversal_mtf | F_SuperReversal_mtf.exe | 5m, higher 1h |
| F_EMA3_SRSI_ATR | F_3EMA_SRSI_ATR.exe | 5m |

Timeframes come from configuration; higher timeframes must be larger integral
multiples. Strategy parameter ranges and eligibility filters remain next to
their signal definitions in each source file. A short study may have no
candidate meeting its minimum trades, gain, or drawdown requirement.

## Data and automatic repair

Research data goes in ignored `data/research/`; downloaded archives and response
records go in ignored `.data_cache/`. The four original raw files are retained
under `data/fixtures/`, never automatically updated. They contain historical
coverage limitations and are not a complete research dataset.

Candles are assembled from checksum-verified Binance monthly ZIPs, daily ZIPs
for the current month, then REST for remaining holes. Existing and merged rows
are checked for timestamps, cadence, finite values, positive prices, OHLC
bounds, nonnegative volume, and conflicting duplicates. Writes use atomic
replacement and download preparation is locked across processes. An invalid
local candle file is preserved as a timestamped .invalid backup before rebuilding. Failed
transfers cannot masquerade as completed coverage. A gap still missing after
repair stops the run; there are no fabricated zero candles.

Futures funding uses published monthly settlement rates and intervals, plus
the paginated REST settlement ledger and settlement mark prices. The sources
are cross-checked to archive precision. Missing REST events present in the
archive require an exact-time mark-price candle to repair them. Missing marks,
source conflicts, and funding holes over eight hours stop the run. Intervals
shorter than eight hours are retained; the legacy `8h` filename is not used to
schedule payments. Cached ledgers record coverage and a content digest. When
archives are unavailable, completeness also depends on the REST source: local
hashes cannot prove that an upstream record was never omitted.

Funding keeps reported millisecond offsets. A settlement exactly at a candle
opening applies before opening orders; one millisecond after it applies after
those orders. Such a settlement is placed before unresolved intrabar stops,
because OHLC data cannot establish that finer ordering. Whole-second event
times select the coarsest supported execution candle that represents them,
automatically downloading finer candles if needed. Unsupported timestamps fail.

Data sources and conventions:
[Binance public data](https://github.com/binance/binance-public-data),
[funding history](https://developers.binance.com/docs/derivatives/usds-margined-futures/market-data/rest-api/Get-Funding-Rate-History).

## Execution and accounting

- Indicators see completed signal candles. Orders execute at the next opening.
  Higher-candle values become usable when that candle completes; warmup values
  never generate signals.
- Protective gaps fill at the opening. Otherwise stops and targets fill at their
  trigger price. A candle touching both uses the stop first and increments the
  ambiguity counter. Intrabar fill timestamps are the candle-close upper bound.
- ATR distances come from the signal candle and are anchored to the actual entry
  fill. Trailing stops update from a completed close and apply on the following
  candle. Holding-time limits start at the actual entry.
- Opening exits across all pairs occur before new allocations. Protective gap
  exits cannot reopen on that candle. No new position opens on the final candle;
  remaining positions close at its close.
- Initial capital is 1,000 USDT and commission is 0.1% per side. Futures entries
  reserve their fee inside the allocated budget. Short exit commission uses the
  actual exit notional. Funding is accounted for separately from commissions.
  A winner has positive net P&L after both fees and funding.
- Equity and drawdown are measured at every execution-candle close, including
  open positions. This cannot measure adverse movement within the candle.
  Futures use a 1x collateral/P&L model with no liquidation, spread, extra
  slippage, order book, or exchange minimum-size model. Nonpositive or
  nonfinite equity invalidates a candidate.

## Search and results

The default search samples 1,000 distinct valid grid combinations with seed 42.
Workers have separate caches; ties resolve by candidate order. `max_trials: 0`
requests exhaustive enumeration, limited to 10 million raw combinations to
bound memory. For larger grids use a finite budget.

The first 80% of execution candles selects parameters. The last 20% is evaluated
once for the selected parameters with a fresh wallet and no carried positions.
Preceding candles remain available for indicator history. Changing the holdout
fraction changes this split. This is one chronological holdout, not proof of
generalization across market regimes.

Each run writes a unique JSON file in `results/`. It records source identities,
bounds, seed, budget, parameters, eligibility thresholds, model assumptions,
build revision, elapsed evaluation time, current resident memory, and selected
training/holdout fills and equity. Intermediate candidates retain scalar
statistics only. `no_eligible_candidate` is a completed search with no qualifying
selection, not a successful trading strategy; it reports the maximum training
trade count and number of invalid candidates.

The selection score is gain / drawdown correction × net win rate. Drawdown
correction is `100 / (1 + DD_percent/100) - 100`; a zero-drawdown run has score
zero. Calmar is `CAGR / abs(DD_fraction)`, using elapsed seconds and a 365.25-day
year. Undefined Calmar is JSON `null`. There is no monthly Calmar or comparison
against stale `*_best.txt` files.

## Representative production check

```sh
docker compose run --rm -T backtest bash tools/smoke_all_strategies.sh --trials 3
```

This runs all configured strategies for 2024-01-01 through 2024-04-01 with real
data and automatic repair. It preserves the user configuration, requires every
run to finish, and reports no-eligible searches separately. Use `--start`,
`--end`, `--trials`, `--timeout`, or `--offline` to change the check.
Logs are in `.smoke_logs/`; full results remain in `results/`.

The independent EMA references accept `--datafile`, `--ema1` (slow), and
`--ema2` (fast), using `python/backtest_double_EMA_float.py` or its
`_numba.py` counterpart inside Compose. They evaluate fixed parameters and
do not tune on their holdout. See [CONTRIBUTING.md](CONTRIBUTING.md) for checks
and code organization.
