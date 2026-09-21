#!/usr/bin/env bash
# Run inside the Compose service. TA-Lib's code-generation build must be serial.
set -euo pipefail
cd "$(dirname "$0")"
if [[ ! -f talib/talib_install/lib/libta_lib.so ]]; then
    [[ -d talib ]] || unzip -q -o talib.zip
    (cd talib && ./configure --prefix="$PWD/talib_install" && make && make install)
fi
make -j"${BUILD_JOBS:-2}" all
echo "Built build/release/. Run bash run_regression.sh in the same service."
