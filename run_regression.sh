#!/usr/bin/env bash
#
# Full regression gate: build everything, run the unit suite, then compare the three
# numeric harnesses against the tracked fixtures in regression_expected/.
#
#   ./run_regression.sh              # verify against tracked fixtures
#   ./run_regression.sh --rebaseline # regenerate the fixtures from this build
#
# Unlike the previous PowerShell-only flow, the strategy list is not duplicated here:
# `make all` discovers strategies itself, so adding a strategy needs no edit to this
# script. run_regression.ps1 is a thin wrapper that calls into this file via WSL.

set -euo pipefail

cd "$(dirname "$0")"

REBASELINE=0
if [[ "${1:-}" == "--rebaseline" ]]; then
    REBASELINE=1
fi

TALIB_LIB="talib/talib_install/lib/libta_lib.so"
if [[ ! -f "$TALIB_LIB" ]]; then
    echo "FATAL: expected TA-Lib at $TALIB_LIB." >&2
    echo "Build it first:  unzip -o talib.zip && cd talib && ./configure --prefix=\$(pwd)/talib_install/ && make && make install" >&2
    echo "(build TA-Lib serially -- its own gen_code target races under 'make -j')" >&2
    exit 1
fi

TMP=".regression_tmp"
EXPECTED="regression_expected"
mkdir -p "$TMP"

# The binaries carry an rpath to the in-tree TA-Lib, but export the path too so a
# stale system libta_lib cannot be picked up instead.
export LD_LIBRARY_PATH="$PWD/talib/talib_install/lib:${LD_LIBRARY_PATH:-}"

strip_ansi() { sed -r 's/\x1B\[[0-9;]*m//g'; }

fail() {
    echo "REGRESSION FAILED: $*" >&2
    exit 1
}

# Compare a produced file against its fixture, or overwrite the fixture under
# --rebaseline. Prints a unified diff on mismatch so drift is readable at a glance.
check_or_rebaseline() {
    local label="$1" actual="$2" expected="$3"

    if [[ $REBASELINE -eq 1 ]]; then
        cp "$actual" "$expected"
        echo "$label: REBASELINED"
        return
    fi

    if [[ ! -f "$expected" ]]; then
        fail "$label: no fixture at $expected (run with --rebaseline to create it)"
    fi

    if ! diff -u "$expected" "$actual"; then
        fail "$label: output differs from $expected"
    fi
    echo "$label: OK"
}

echo "=== build ==="
make -j"$(nproc)" all

echo
echo "=== unit tests ==="
./tests.exe || fail "tests.exe returned non-zero"

echo
echo "=== downloader unit tests ==="
# Pure logic only -- no network, so this is safe to gate every push on.
if command -v python3 > /dev/null; then
    python3 -m unittest discover -s tools -p 'test_*.py' -q || fail "downloader tests failed"
else
    echo "  python3 not found, skipping"
fi

echo
echo "=== verification_regression ==="
./verification_regression.exe > "$TMP/verification_regression.txt"
strip_ansi < "$TMP/verification_regression.txt" \
    | grep -E '^(SPOT_ROWS|EMA11_TAIL_SUM|EMA448_TAIL_SUM|BBANDS_TAIL_SUM|VOLUME_TAIL_SUM|ALIGN_STARTS|ALIGN_LAST_TS|ALIGN_LAST_CLOSE_SUM|FUTURES_ROWS|FUNDING_ROWS|FUNDING_FIRST|FUNDING_SECOND|FUNDING_MISS)' \
    > "$TMP/verification_filtered.txt"
check_or_rebaseline "verification_regression" \
    "$TMP/verification_filtered.txt" "$EXPECTED/verification_regression.expected.txt"

echo
echo "=== strategy_regression ==="
./strategy_regression.exe > "$TMP/strategy_regression.txt"
strip_ansi < "$TMP/strategy_regression.txt" \
    | grep -E '^MULTI_SPOT_' \
    > "$TMP/strategy_filtered.txt"
check_or_rebaseline "strategy_regression" \
    "$TMP/strategy_filtered.txt" "$EXPECTED/strategy_regression.expected.txt"

echo
echo "=== backtest_double_EMA_float (full sweep) ==="
./backtest_double_EMA_float.exe > "$TMP/backtest_double_EMA_float.txt"
# The banner repeats as the sweep progresses; keep only the final block's metrics.
strip_ansi < "$TMP/backtest_double_EMA_float.txt" \
    | grep -E '^ (EMAs|Gain|Win rate|max DD|Gain over DDC|Score|Number of trades) +:' \
    | sed 's/^ //' \
    | tail -7 \
    > "$TMP/backtest_filtered.txt"
check_or_rebaseline "backtest_double_EMA_float" \
    "$TMP/backtest_filtered.txt" "$EXPECTED/backtest_double_EMA_float.expected.txt"

echo
echo "ALL REGRESSION CHECKS PASSED"
