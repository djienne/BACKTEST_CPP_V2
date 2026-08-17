#!/usr/bin/env bash
#
# One-shot setup: build the bundled TA-Lib, then build the backtester.
# Equivalent to the "Getting started" steps in README.md.

set -euo pipefail

cd "$(dirname "$0")"

if [ ! -d talib ]; then
    echo "==> extracting talib.zip"
    unzip -q -o talib.zip
fi

echo "==> building TA-Lib"
cd talib
./configure --prefix="$(pwd)/talib_install/"
# Serial on purpose: TA-Lib's own src/tools/gen_code target races under make -j
# and fails intermittently.
make
make install
cd ..

echo "==> building the backtester"
make -j"$(nproc)" all

echo
echo "Done. Run ./run_regression.sh to verify the build."
echo "No LD_LIBRARY_PATH needed: the binaries carry an rpath to talib/talib_install/lib."
