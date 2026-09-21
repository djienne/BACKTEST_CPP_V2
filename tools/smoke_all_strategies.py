"""Bounded real-data production check. Downloads are enabled unless --offline is set.

Uses a disposable config, preserving the user's search budget and dates. A completed
search without an eligible candidate is reported honestly, not called a trading success.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT/"backtest_config.json")
    parser.add_argument("--start", default="2024-01-01")
    parser.add_argument("--end", default="2024-04-01")
    parser.add_argument("--trials", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    if args.trials < 1 or args.timeout < 1:
        parser.error("trials and timeout must be positive")
    config = json.loads(args.config.read_text())
    config["run"].update(start=args.start, end=args.end, max_trials=args.trials, seed=42)
    logs = ROOT/".smoke_logs"
    logs.mkdir(exist_ok=True)
    failures = 0
    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp)/"config.json"
        path.write_text(json.dumps(config))
        env = dict(os.environ, BACKTEST_CONFIG=str(path))
        if args.offline:
            env["BACKTEST_OFFLINE"] = "1"
        for name, strategy in config["strategies"].items():
            before = set((ROOT/"results").glob("*.json"))
            log = logs/(name+".log")
            with log.open("w") as out:
                try:
                    completed = subprocess.run([str(ROOT/"build/release"/(strategy["file"]+".exe"))],
                        cwd=ROOT, env=env, stdout=out, stderr=subprocess.STDOUT, timeout=args.timeout)
                    code = completed.returncode
                except subprocess.TimeoutExpired:
                    code = 124
            produced = set((ROOT/"results").glob("*.json"))-before
            matches = [p for p in produced if json.loads(p.read_text()).get("strategy") == name]
            if code or len(matches) != 1:
                print(f"{name}: FAILED ({code}); see {log}", flush=True)
                failures += 1
                continue
            report = json.loads(matches[0].read_text())
            if not 0 < report["trials"] <= args.trials or report["invalid_candidates"]:
                failures += 1
            if report["status"] == "complete" and not report["holdout"]["valid"]:
                failures += 1
            print(f'{name}: {report["status"]}; trials={report["trials"]}; '
                  f'max training trades={report["max_training_trades"]}; result={matches[0].name}', flush=True)
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
