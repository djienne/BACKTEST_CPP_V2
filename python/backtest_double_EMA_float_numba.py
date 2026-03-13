import argparse
import math
import os
import time

import numpy as np
from numba import njit

from double_ema_common import (
    FEE_PC,
    INITIAL_USDT,
    MAX_ALLOWED_DD,
    MAX_EMA,
    MIN_EMA,
    MIN_TRADES,
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
STRAT_NAME = "2-EMA crossover simple (Numba benchmark)"


@njit(cache=True)
def evaluate_pair_numba(close: np.ndarray, ema1: np.ndarray, ema2: np.ndarray, ii_begin: int, fee_pc: float, initial_usdt: float) -> tuple:
    in_position = False
    in_long = False
    nb_profit = 0
    nb_loss = 0
    nb_positions = 0
    usdt_amount = initial_usdt
    max_usdt_amount = usdt_amount
    max_drawdown = 0.0
    price_position_open = 0.0
    total_fees_paid = 0.0

    for ii in range(ii_begin, close.shape[0]):
        last_iteration = ii == close.shape[0] - 1
        current_price = close[ii]

        if usdt_amount > max_usdt_amount:
            max_usdt_amount = usdt_amount

        pc_change_with_max = (usdt_amount - max_usdt_amount) / max_usdt_amount * 100.0
        if pc_change_with_max < max_drawdown:
            max_drawdown = pc_change_with_max

        if in_position and in_long:
            usdt_amount_check = usdt_amount + (current_price - price_position_open) / price_position_open * usdt_amount
            if usdt_amount_check <= 0.0:
                return 0.0, 0.0, -100.0, 0.0, 0.0, nb_positions, 0.0, total_fees_paid

        open_long = ema2[ii] >= ema1[ii] and ema2[ii - 1] <= ema1[ii - 1]
        close_long = ema2[ii] <= ema1[ii] and ema2[ii - 1] >= ema1[ii - 1]

        if in_position and in_long and (close_long or last_iteration):
            amount_before = usdt_amount
            usdt_amount = usdt_amount + (current_price - price_position_open) / price_position_open * usdt_amount
            fee = current_price / price_position_open * amount_before * fee_pc / 100.0
            usdt_amount -= fee
            total_fees_paid += fee
            in_position = False
            in_long = False
            if current_price > price_position_open:
                nb_profit += 1
            else:
                nb_loss += 1

        if (not in_position) and open_long and (not last_iteration):
            price_position_open = current_price
            fee = usdt_amount * fee_pc / 100.0
            usdt_amount -= fee
            total_fees_paid += fee
            in_position = True
            in_long = True
            nb_positions += 1

    gain = (usdt_amount - initial_usdt) / initial_usdt * 100.0
    if nb_positions == 0:
        return gain, -1e30, 0.0, max_drawdown, -1e30, 0, usdt_amount, total_fees_paid

    win_rate = float(nb_profit) / float(nb_positions) * 100.0
    ddc = (1.0 / (1.0 + max_drawdown / 100.0) - 1.0) * 100.0
    gain_over_ddc = gain / ddc
    score = gain_over_ddc * win_rate
    return gain, gain_over_ddc, win_rate, max_drawdown, score, nb_positions, usdt_amount, total_fees_paid


@njit(cache=True)
def search_numba(close: np.ndarray, ema_matrix: np.ndarray, periods: np.ndarray, pair_i: np.ndarray, pair_j: np.ndarray, ii_begin: int) -> tuple:
    best_score = -1e30
    best_ema1 = -1
    best_ema2 = -1
    best_gain = 0.0
    best_gain_over_ddc = 0.0
    best_win_rate = 0.0
    best_max_dd = 0.0
    best_trades = 0
    best_wallet = 0.0
    best_fees = 0.0

    for idx in range(pair_i.shape[0]):
        idx_a = pair_i[idx]
        idx_b = pair_j[idx]
        gain, gain_over_ddc, win_rate, max_dd, score, trades, wallet, fees = evaluate_pair_numba(
            close, ema_matrix[idx_a], ema_matrix[idx_b], ii_begin, FEE_PC, INITIAL_USDT
        )
        if score > best_score and gain < 10000.0 and trades >= MIN_TRADES and max_dd > MAX_ALLOWED_DD:
            best_score = score
            best_ema1 = periods[idx_a]
            best_ema2 = periods[idx_b]
            best_gain = gain
            best_gain_over_ddc = gain_over_ddc
            best_win_rate = win_rate
            best_max_dd = max_dd
            best_trades = trades
            best_wallet = wallet
            best_fees = fees

    return (
        best_ema1,
        best_ema2,
        best_gain,
        best_gain_over_ddc,
        best_win_rate,
        best_max_dd,
        best_score,
        best_trades,
        best_wallet,
        best_fees,
    )


def to_result(values: tuple) -> RunResult:
    return RunResult(
        ema1=int(values[0]),
        ema2=int(values[1]),
        gain_pc=float(values[2]),
        gain_over_ddc=float(values[3]),
        win_rate=float(values[4]),
        max_dd=float(values[5]),
        score=float(values[6]),
        nb_posi_entered=int(values[7]),
        wallet_val_usdt=float(values[8]),
        total_fees_paid=float(values[9]),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=STRAT_NAME)
    parser.add_argument("--datafile", default=DATAFILE)
    parser.add_argument("--sample-size", type=int, default=2048, help="Number of EMA pairs to test. Use 0 for the full grid.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--ema1", type=int, default=None, help="Evaluate a single EMA1 / EMA2 pair in both backends.")
    parser.add_argument("--ema2", type=int, default=None)
    parser.add_argument("--skip-python", action="store_true", help="Only run the numba backend. Useful for the full grid.")
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

    if args.ema1 is not None:
        if args.ema1 < MIN_EMA or args.ema1 > MAX_EMA or args.ema2 < MIN_EMA or args.ema2 > MAX_EMA:
            raise ValueError("EMA values out of supported range.")
        idx_a = int(args.ema1 - MIN_EMA)
        idx_b = int(args.ema2 - MIN_EMA)

        python_result = None
        python_seconds = 0.0
        if not args.skip_python:
            t_python = time.perf_counter()
            python_result = evaluate_pair_python(close, ema_matrix[idx_a], ema_matrix[idx_b], ii_begin)
            python_seconds = time.perf_counter() - t_python
            python_result.ema1 = args.ema1
            python_result.ema2 = args.ema2

        evaluate_pair_numba(close, ema_matrix[idx_a], ema_matrix[idx_b], ii_begin, FEE_PC, INITIAL_USDT)
        t_numba = time.perf_counter()
        numba_result = to_result((args.ema1, args.ema2, *evaluate_pair_numba(close, ema_matrix[idx_a], ema_matrix[idx_b], ii_begin, FEE_PC, INITIAL_USDT)))
        numba_seconds = time.perf_counter() - t_numba
        mode = "single pair"
    else:
        python_result = None
        python_seconds = 0.0
        if not args.skip_python:
            t_python = time.perf_counter()
            python_result = search_python(close, ema_matrix, periods, pair_i, pair_j, ii_begin)
            python_seconds = time.perf_counter() - t_python

        search_numba(close, ema_matrix, periods, pair_i, pair_j, ii_begin)
        t_numba = time.perf_counter()
        numba_result = to_result(search_numba(close, ema_matrix, periods, pair_i, pair_j, ii_begin))
        numba_seconds = time.perf_counter() - t_numba
        mode = f"search over {pair_i.shape[0]} pairs"

    print(f"Mode: {mode}")
    print_run_result("Numba", numba_result)
    if python_result is not None:
        print_run_result("Pure Python", python_result)
        compare_results(python_result, numba_result, "python", "numba")
        search_speedup = python_seconds / numba_seconds if numba_seconds > 0.0 else math.inf
        total_python = prep_seconds + python_seconds
        total_numba = prep_seconds + numba_seconds
        total_speedup = total_python / total_numba if total_numba > 0.0 else math.inf
        print(f"Python search time : {python_seconds:.3f}s")
        print(f"Search speedup     : {search_speedup:.2f}x")
        print(f"Python total time  : {total_python:.3f}s")
        print(f"End-to-end speedup : {total_speedup:.2f}x")
    print(f"Numba search time  : {numba_seconds:.3f}s")
    print(f"Numba total time   : {prep_seconds + numba_seconds:.3f}s")


if __name__ == "__main__":
    main()
