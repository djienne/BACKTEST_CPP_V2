# Windows entry point; builds and verification run in Docker Compose.
param([ValidateSet("release", "debug", "asan", "tsan")][string]$Mode = "release")
$ErrorActionPreference = "Stop"
Push-Location $PSScriptRoot
try {
    if ($Mode -in @("asan", "tsan")) {
        docker compose -f compose.yaml -f compose.sanitizers.yaml run --build --rm -T backtest setarch x86_64 -R bash run_regression.sh $Mode
    } else {
        docker compose run --build --rm -T backtest bash run_regression.sh $Mode
    }
    if ($LASTEXITCODE -ne 0) { throw "Regression checks failed: $LASTEXITCODE" }
} finally { Pop-Location }
