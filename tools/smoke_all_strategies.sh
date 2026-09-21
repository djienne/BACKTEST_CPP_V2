#!/usr/bin/env bash
# Run in Compose; each strategy must finish. Timeouts count as failures.
set -euo pipefail
cd "$(dirname "$0")/.."
exec python3 tools/smoke_all_strategies.py "$@"
