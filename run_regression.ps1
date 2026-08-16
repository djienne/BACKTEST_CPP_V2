# Thin Windows wrapper around run_regression.sh.
#
# The build and all numeric comparisons live in run_regression.sh so there is exactly
# one source of truth; this script only translates the repo path into its WSL form and
# forwards any arguments (e.g. -Rebaseline).
#
#   powershell -File .\run_regression.ps1
#   powershell -File .\run_regression.ps1 -Rebaseline

param(
    [switch]$Rebaseline
)

$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$LinuxRepoRoot = "/mnt/" + $RepoRoot.Substring(0, 1).ToLower() + $RepoRoot.Substring(2).Replace("\", "/")

$ScriptArgs = @("$LinuxRepoRoot/run_regression.sh")
if ($Rebaseline) {
    $ScriptArgs += "--rebaseline"
}

wsl -e bash @ScriptArgs
if ($LASTEXITCODE -ne 0) {
    throw "run_regression.sh failed with exit code $LASTEXITCODE"
}
