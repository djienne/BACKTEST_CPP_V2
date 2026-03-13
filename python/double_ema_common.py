import math
from dataclasses import dataclass

import numpy as np


FEE_PC = 0.1
INITIAL_USDT = 1000.0
MIN_EMA = 3
MAX_EMA = 600
MIN_EMA_GAP = 7
MIN_TRADES = 100
MAX_ALLOWED_DD = -40.0


@dataclass
class RunResult:
    ema1: int
    ema2: int
    gain_pc: float
    gain_over_ddc: float
    win_rate: float
    max_dd: float
    score: float
    nb_posi_entered: int
    wallet_val_usdt: float
    total_fees_paid: float


def load_close_series(path: str) -> np.ndarray:
    raw = np.loadtxt(path, delimiter=",", skiprows=1, usecols=(4,))
    return raw.astype(np.float64)


def generate_periods() -> np.ndarray:
    return np.arange(MIN_EMA, MAX_EMA + 1, dtype=np.int32)


def generate_param_pairs(sample_size: int | None, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    periods = generate_periods()
    left = []
    right = []
    for i, ema1 in enumerate(periods):
        for j, ema2 in enumerate(periods):
            if abs(int(ema1) - int(ema2)) < MIN_EMA_GAP or ema1 < ema2:
                continue
            left.append(i)
            right.append(j)

    left = np.asarray(left, dtype=np.int32)
    right = np.asarray(right, dtype=np.int32)

    if sample_size is not None and sample_size < left.shape[0]:
        rng = np.random.default_rng(seed)
        selection = rng.choice(left.shape[0], size=sample_size, replace=False)
        left = left[selection]
        right = right[selection]

    return periods, left, right


def compute_ema(close: np.ndarray, period: int) -> np.ndarray:
    out = np.zeros(close.shape[0], dtype=np.float64)
    if period <= 0 or period > close.shape[0]:
        return out

    alpha = 2.0 / (period + 1.0)
    start = period - 1
    out[start] = np.mean(close[:period])
    for i in range(start + 1, close.shape[0]):
        out[i] = alpha * close[i] + (1.0 - alpha) * out[i - 1]
    return out


def compute_ema_matrix(close: np.ndarray, periods: np.ndarray) -> np.ndarray:
    out = np.empty((periods.shape[0], close.shape[0]), dtype=np.float64)
    for idx, period in enumerate(periods):
        out[idx, :] = compute_ema(close, int(period))
    return out


def evaluate_pair_python(close: np.ndarray, ema1: np.ndarray, ema2: np.ndarray, ii_begin: int) -> RunResult:
    in_position = False
    in_long = False
    nb_profit = 0
    nb_loss = 0
    nb_positions = 0
    usdt_amount = INITIAL_USDT
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
                return RunResult(0, 0, 0.0, 0.0, 0.0, -100.0, 0.0, nb_positions, 0.0, total_fees_paid)

        open_long = ema2[ii] >= ema1[ii] and ema2[ii - 1] <= ema1[ii - 1]
        close_long = ema2[ii] <= ema1[ii] and ema2[ii - 1] >= ema1[ii - 1]

        if in_position and in_long and (close_long or last_iteration):
            amount_before = usdt_amount
            usdt_amount = usdt_amount + (current_price - price_position_open) / price_position_open * usdt_amount
            fee = current_price / price_position_open * amount_before * FEE_PC / 100.0
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
            fee = usdt_amount * FEE_PC / 100.0
            usdt_amount -= fee
            total_fees_paid += fee
            in_position = True
            in_long = True
            nb_positions += 1

    gain = (usdt_amount - INITIAL_USDT) / INITIAL_USDT * 100.0
    if nb_positions == 0:
        return RunResult(0, 0, gain, -1e30, 0.0, max_drawdown, -1e30, 0, usdt_amount, total_fees_paid)

    win_rate = float(nb_profit) / float(nb_positions) * 100.0
    ddc = (1.0 / (1.0 + max_drawdown / 100.0) - 1.0) * 100.0
    gain_over_ddc = gain / ddc
    score = gain_over_ddc * win_rate
    return RunResult(0, 0, gain, gain_over_ddc, win_rate, max_drawdown, score, nb_positions, usdt_amount, total_fees_paid)


def search_python(close: np.ndarray, ema_matrix: np.ndarray, periods: np.ndarray, pair_i: np.ndarray, pair_j: np.ndarray, ii_begin: int) -> RunResult | None:
    best = None
    best_score = -1e30
    for idx_a, idx_b in zip(pair_i, pair_j):
        result = evaluate_pair_python(close, ema_matrix[idx_a], ema_matrix[idx_b], ii_begin)
        result.ema1 = int(periods[idx_a])
        result.ema2 = int(periods[idx_b])
        if result.score > best_score and result.gain_pc < 10000.0 and result.nb_posi_entered >= MIN_TRADES and result.max_dd > MAX_ALLOWED_DD:
            best = result
            best_score = result.score
    return best


def print_run_result(label: str, result: RunResult) -> None:
    print(f"{label}:")
    print(f"  EMA1 / EMA2 : {result.ema1} / {result.ema2}")
    print(f"  Gain %      : {result.gain_pc:.6f}")
    print(f"  Win rate %  : {result.win_rate:.6f}")
    print(f"  Max DD %    : {result.max_dd:.6f}")
    print(f"  Gain / DDC  : {result.gain_over_ddc:.6f}")
    print(f"  Score       : {result.score:.6f}")
    print(f"  Trades      : {result.nb_posi_entered}")
    print(f"  Wallet USDT : {result.wallet_val_usdt:.6f}")
    print(f"  Fees USDT   : {result.total_fees_paid:.6f}")


def compare_results(a: RunResult, b: RunResult, label_a: str, label_b: str) -> None:
    print(f"Result deltas ({label_a} - {label_b}):")
    print(f"  EMA1 / EMA2 : {a.ema1 - b.ema1} / {a.ema2 - b.ema2}")
    print(f"  Gain %      : {a.gain_pc - b.gain_pc:.10f}")
    print(f"  Win rate %  : {a.win_rate - b.win_rate:.10f}")
    print(f"  Max DD %    : {a.max_dd - b.max_dd:.10f}")
    print(f"  Gain / DDC  : {a.gain_over_ddc - b.gain_over_ddc:.10f}")
    print(f"  Score       : {a.score - b.score:.10f}")
    print(f"  Trades      : {a.nb_posi_entered - b.nb_posi_entered}")
    print(f"  Wallet USDT : {a.wallet_val_usdt - b.wallet_val_usdt:.10f}")
    print(f"  Fees USDT   : {a.total_fees_paid - b.total_fees_paid:.10f}")
