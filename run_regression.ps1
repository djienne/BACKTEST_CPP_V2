$ErrorActionPreference = "Stop"

function Invoke-Bash {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Command
    )

    bash -lc $Command
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed: $Command"
    }
}

function Strip-Ansi {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Line
    )

    return ($Line -replace "\x1B\[[0-9;]*m", "").Trim()
}

function Assert-LinesEqual {
    param(
        [Parameter(Mandatory = $true)]
        [string[]]$Actual,
        [Parameter(Mandatory = $true)]
        [string[]]$Expected,
        [Parameter(Mandatory = $true)]
        [string]$Label
    )

    if ($Actual.Count -ne $Expected.Count) {
        throw "$Label line count mismatch. Expected $($Expected.Count), got $($Actual.Count)."
    }

    for ($i = 0; $i -lt $Expected.Count; $i++) {
        if ($Actual[$i] -ne $Expected[$i]) {
            throw "$Label mismatch at line $($i + 1). Expected '$($Expected[$i])', got '$($Actual[$i])'."
        }
    }
}

function Get-LastMetricValue {
    param(
        [Parameter(Mandatory = $true)]
        [object]$Content,
        [Parameter(Mandatory = $true)]
        [string]$Prefix,
        [Parameter(Mandatory = $true)]
        [string]$Pattern
    )

    $lines = @($Content) | Where-Object { $_ -ne $null }
    $line = $lines | Where-Object { $_ -like " $Prefix*" } | Select-Object -Last 1
    if (-not $line) {
        throw "Missing metric line for '$Prefix'."
    }

    $normalized = Strip-Ansi $line
    $normalized = $normalized.TrimStart()
    return ($normalized -replace $Pattern, "")
}

$RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$TempDir = Join-Path $RepoRoot ".regression_tmp"
New-Item -ItemType Directory -Force -Path $TempDir | Out-Null

$LinuxRepoRoot = "/mnt/" + $RepoRoot.Substring(0, 1).ToLower() + $RepoRoot.Substring(2).Replace("\", "/")

if (-not (Test-Path (Join-Path $RepoRoot "talib\talib_install\lib\libta_lib.so"))) {
    throw "Expected TA-Lib at talib/talib_install/lib/libta_lib.so. Extract/build TA-Lib first."
}

# Keep the build flags here in sync with Makefile's CXXFLAGS_OPT. We still invoke g++
# directly via WSL so a missing make or stale build state doesn't mask real failures.
$CXXFLAGS = "-O3 -fno-trapping-math -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare"
$COMMON   = "./custom_talib_wrapper.cpp ./tools.cpp ./trade_core.cpp"
$TALIB    = "-I./talib/talib_install/include -L./talib/talib_install/lib -Wl,-rpath,./talib/talib_install/lib -lta_lib -lpthread"

function Build-Target {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Source,
        [Parameter(Mandatory = $true)]
        [string]$Output
    )
    Invoke-Bash "cd '$LinuxRepoRoot' && g++ $CXXFLAGS -I./talib/talib_install/include $COMMON ./$Source -L./talib/talib_install/lib -Wl,-rpath,./talib/talib_install/lib -lta_lib -lpthread -o ./$Output"
}

# Tests first: fastest failure surface, no data-file dependency.
Build-Target -Source "tests.cpp" -Output "tests.exe"
Invoke-Bash "cd '$LinuxRepoRoot' && ./tests.exe"

# Regression drivers + all strategy binaries.
$AllSources = @(
    "verification_regression.cpp",
    "backtest_double_EMA_float.cpp",
    "backtest_double_EMA_StochRSI_float_muti_pair.cpp",
    "BigWill.cpp",
    "F_BigWill.cpp",
    "BBTREND.cpp",
    "F_BBTREND.cpp",
    "backtest_TRIX_multi_pair.cpp",
    "SuperReversal_mtf.cpp",
    "F_SuperReversal_mtf.cpp",
    "SuperTrend_EMA_ATR.cpp",
    "3EMA_SRSI_ATR.cpp",
    "F_3EMA_SRSI_ATR.cpp",
    "strategy_regression.cpp"
)
foreach ($src in $AllSources) {
    $exe = [IO.Path]::ChangeExtension($src, ".exe")
    Build-Target -Source $src -Output $exe
}

$VerificationOutput = Join-Path $TempDir "verification_regression.txt"
$BacktestOutput = Join-Path $TempDir "backtest_double_EMA_float.txt"
$StrategyRegressionOutput = Join-Path $TempDir "strategy_regression.txt"

Invoke-Bash "cd '$LinuxRepoRoot' && ./verification_regression.exe > ./.regression_tmp/verification_regression.txt"
Invoke-Bash "cd '$LinuxRepoRoot' && ./backtest_double_EMA_float.exe > ./.regression_tmp/backtest_double_EMA_float.txt"
Invoke-Bash "cd '$LinuxRepoRoot' && ./strategy_regression.exe > ./.regression_tmp/strategy_regression.txt"

$VerificationActual = Get-Content $VerificationOutput | ForEach-Object { Strip-Ansi $_ } | Where-Object {
    $_ -match "^(SPOT_ROWS|EMA11_TAIL_SUM|EMA448_TAIL_SUM|BBANDS_TAIL_SUM|ALIGN_STARTS|ALIGN_LAST_TS|ALIGN_LAST_CLOSE_SUM|FUTURES_ROWS|FUNDING_ROWS|FUNDING_FIRST|FUNDING_SECOND|FUNDING_MISS)"
}
$VerificationExpected = Get-Content (Join-Path $RepoRoot "regression_expected\verification_regression.expected.txt")
Assert-LinesEqual -Actual $VerificationActual -Expected $VerificationExpected -Label "verification_regression"

$BacktestContent = Get-Content $BacktestOutput
$BacktestActual = @(
    "EMAs             : " + (Get-LastMetricValue -Content $BacktestContent -Prefix "EMAs             :" -Pattern "^EMAs\s*:\s*")
    "Gain             : " + (Get-LastMetricValue -Content $BacktestContent -Prefix "Gain             :" -Pattern "^Gain\s*:\s*")
    "Win rate         : " + (Get-LastMetricValue -Content $BacktestContent -Prefix "Win rate         :" -Pattern "^Win rate\s*:\s*")
    "max DD           : " + (Get-LastMetricValue -Content $BacktestContent -Prefix "max DD           :" -Pattern "^max DD\s*:\s*")
    "Gain over DDC    : " + (Get-LastMetricValue -Content $BacktestContent -Prefix "Gain over DDC    :" -Pattern "^Gain over DDC\s*:\s*")
    "Score            : " + (Get-LastMetricValue -Content $BacktestContent -Prefix "Score            :" -Pattern "^Score\s*:\s*")
    "Number of trades : " + (Get-LastMetricValue -Content $BacktestContent -Prefix "Number of trades :" -Pattern "^Number of trades\s*:\s*")
)
$BacktestExpected = Get-Content (Join-Path $RepoRoot "regression_expected\backtest_double_EMA_float.expected.txt")
Assert-LinesEqual -Actual $BacktestActual -Expected $BacktestExpected -Label "backtest_double_EMA_float"

$StrategyRegressionActual = Get-Content $StrategyRegressionOutput | ForEach-Object { Strip-Ansi $_ } | Where-Object {
    $_ -match "^(MULTI_SPOT_EMA_SHORT|MULTI_SPOT_EMA_LONG|MULTI_SPOT_GAIN|MULTI_SPOT_WIN_RATE|MULTI_SPOT_MAX_DD|MULTI_SPOT_GAIN_OVER_DDC|MULTI_SPOT_SCORE|MULTI_SPOT_TRADES)"
}
$StrategyRegressionExpected = Get-Content (Join-Path $RepoRoot "regression_expected\strategy_regression.expected.txt")
Assert-LinesEqual -Actual $StrategyRegressionActual -Expected $StrategyRegressionExpected -Label "strategy_regression"

Write-Host "tests: OK"
Write-Host "verification_regression: OK"
Write-Host "backtest_double_EMA_float: OK"
Write-Host "strategy_regression: OK"
foreach ($src in $AllSources) {
    $exe = [IO.Path]::ChangeExtension($src, ".exe")
    Write-Host "$exe build: OK"
}
