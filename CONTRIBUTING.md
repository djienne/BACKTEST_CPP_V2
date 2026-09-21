# Contributing

The model and its evidence matter more than preserving old output. Keep changes
small, reuse the existing executor and indicators, and avoid additional
frameworks or dependencies unless a concrete research need requires them.

## Where behavior belongs

| File | Responsibility |
|---|---|
| Strategy .cpp files | Parameter ranges, indicator choices, completed-bar signals, eligibility |
| strategy_runner.hh | Deterministic candidate selection, worker-local caches, chronological holdout, result records |
| trade_core.hh / .cpp | Orders, protection, fees, funding, equity and metrics |
| indicators.cpp / .hh | Indicator definitions, readiness, UTC resampling and projection |
| data_io.cpp / .hh | Strict loading, common timestamp interval, data preparation launch |
| tools/download_data.py | Source fetching, gap repair, funding coverage and atomic data writes |
| config.cpp / .hh | Shared configuration validation |
| tools.cpp / .hh | Small shared numeric/calendar helpers and result types |

A strategy returns an `Intent` for a completed signal candle. It does not write
its own portfolio loop. Compute indicators outside the candle loop and use
`Indicators::ready` through `strategy_runner::evaluate`. Preserve warmup
history without trading on unavailable values. Never resample by row count
without checking timestamps.

Position sizes, cash and accumulated statistics use double precision; OHLC and
published indicator arrays use float. Funding timestamps preserve milliseconds
in addition to the whole-second execution clock.

## Validation

Run verification through Docker Compose as shown in the README. A source
change should leave the smallest check that fails for the actual defect.
Prefer hand-computed cash flows, units, limiting cases, causal perturbations
and independent calculations over snapshots copied from the implementation.
Expected numerical values should come from a separate calculation, not another
wrapper around the same routine. Include nonzero cases so an all-zero stub fails.
For causal tests, verify that the perturbation changes the affected output while
leaving the earlier output unchanged. The number of passing assertions is not a
measure of model validity.

Current checks include:

- `tests.cpp`: analytic indicator cases, unit conventions, calendar boundaries,
  portfolio arithmetic, resampling alignment and higher-candle availability.
- `execution_tests.cpp`: the production executor with hand-calculated long and
  short accounts, fee-losing winners, stops, targets, gaps, trailing timing,
  funding order and marks, open-equity drawdown, insolvency, pair allocation,
  worker determinism and holdout isolation.
- `tools/test_download_data.py`: actual repair/merge/write logic with substituted
  upstream transport, including unavailable holes, corrupt ZIPs, conflicting
  candles, missing funding and interruption-safe writes. A repaired CSV is
  consumed by the real C++ loader.
- `tools/test_strategy_reference.py`: the production EMA evaluator against an
  independently written Python ledger, comparing individual fills, equity and
  metrics. Numba parity checks compilation of that same reference; it is not
  separate scientific evidence.
- `tools/smoke_all_strategies.sh`: bounded, completed runs of every production
  strategy on real Binance data. It requires network access on a cold cache.

Passing these checks establishes behavior under the stated model, not trading
edge. For a research claim, inspect fills and equity, evaluate untouched data,
vary fees and execution assumptions, and investigate sensitivity to period,
pair universe and ambiguous candles. Do not call a no-eligible search a trade
success, or a timeout a passed run.

## Builds and maintenance

`make BigWill` builds `build/release/BigWill.exe`; `make BUILD=asan BigWill`
builds a separate instrumented executable. Header dependencies are tracked.
Do not run an old root-level executable after rebuilding another mode.
The sanitizer gates include the concurrent search path in the execution tests.
TA-Lib is an external, uninstrumented library; sanitizer coverage is principally
the project code.

Use the checked-in clang-format style on files you change. `make format`
formats all root C++ sources and headers and is broader than most fixes need.
Run it through Compose when that scope is intentional.

Keep README commands, configuration examples, comments and docstrings aligned
with behavior. Remove obsolete examples instead of preserving misleading
benchmarks. The former generated-output golden files and `docker_easy.zip`
starter were retired; do not reintroduce a rebaseline switch. Raw historical
fixtures are retained unchanged and must not be overwritten by research
downloads. TA-Lib's bundled third-party documentation belongs to that source
distribution.

Before submitting, inspect the diff, run the relevant behavioral checks, and
record any unavailable verification explicitly. Do not publish credentials,
research downloads, build products, or unrelated worktree changes.
