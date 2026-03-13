import argparse
import os
import time

from double_ema_common import (
    MAX_EMA,
    MIN_EMA,
    RunResult,
    compare_results,
    compute_ema_matrix,
    evaluate_pair_python,
    generate_param_pairs,
    load_close_series,
    print_run_result,
    search_python,
)


abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

DATAFILE = "../data/data/binance/1h/BTC-USDT.csv"
STRAT_NAME = "2-EMA crossover simple (Python)"


def main() -> None:
    parser = argparse.ArgumentParser(description=STRAT_NAME)
    parser.add_argument("--datafile", default=DATAFILE)
    parser.add_argument("--sample-size", type=int, default=2048, help="Number of EMA pairs to test. Use 0 for the full grid.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--ema1", type=int, default=None, help="Evaluate a single EMA1 / EMA2 pair instead of searching.")
    parser.add_argument("--ema2", type=int, default=None)
    args = parser.parse_args()

    close = load_close_series(args.datafile)
    periods, pair_i, pair_j = generate_param_pairs(None if args.sample_size == 0 else args.sample_size, args.seed)
    ii_begin = MAX_EMA + 2

    print(f"Strategy: {STRAT_NAME}")
    print(f"Data file: {args.datafile}")
    print(f"Rows loaded: {close.shape[0]}")

    t_prepare = time.perf_counter()
    ema_matrix = compute_ema_matrix(close, periods)
    prep_seconds = time.perf_counter() - t_prepare
    print(f"EMA prep time: {prep_seconds:.3f}s")

    if (args.ema1 is None) != (args.ema2 is None):
        raise ValueError("Provide both --ema1 and --ema2 together.")

    t_search = time.perf_counter()
    if args.ema1 is not None:
        if args.ema1 < MIN_EMA or args.ema1 > MAX_EMA or args.ema2 < MIN_EMA or args.ema2 > MAX_EMA:
            raise ValueError("EMA values out of supported range.")
        idx_a = int(args.ema1 - MIN_EMA)
        idx_b = int(args.ema2 - MIN_EMA)
        result = evaluate_pair_python(close, ema_matrix[idx_a], ema_matrix[idx_b], ii_begin)
        result.ema1 = args.ema1
        result.ema2 = args.ema2
        mode = "single pair"
    else:
        result = search_python(close, ema_matrix, periods, pair_i, pair_j, ii_begin)
        mode = f"search over {pair_i.shape[0]} pairs"

    elapsed = time.perf_counter() - t_search

    if result is None:
        raise RuntimeError("No valid result found.")

    print(f"Mode: {mode}")
    print_run_result("Pure Python", result)
    print(f"Python runtime: {elapsed:.3f}s")


if __name__ == "__main__":
    main()
