#!/usr/bin/env bash
# Offline behavioral gate; expected results come from ledgers and independent methods.
set -euo pipefail
cd "$(dirname "$0")"
mode="${1:-release}"
case "$mode" in release|debug|asan|tsan) ;; *) echo "Use release, debug, asan, or tsan" >&2; exit 2;; esac
if [[ ! -f talib/talib_install/lib/libta_lib.so ]]; then
    bash install.sh
fi
if [[ "$mode" == release ]]; then
    make -j"${BUILD_JOBS:-2}" all
else
    make -j"${BUILD_JOBS:-2}" BUILD="$mode" tests execution_tests
fi
export ASAN_OPTIONS="detect_leaks=0:halt_on_error=1"
export UBSAN_OPTIONS="halt_on_error=1"
export TSAN_OPTIONS="halt_on_error=1"
timeout 120 "./build/$mode/tests.exe"
timeout 120 "./build/$mode/execution_tests.exe"
if [[ "$mode" == release ]]; then
    python3 -m unittest discover -s tools -p 'test_*.py' -q
fi
echo "$mode behavioral checks passed"
