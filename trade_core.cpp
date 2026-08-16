#include "trade_core.hh"

namespace trade_core
{
void record_wallet_snapshot(const double wallet_val_usdt, const int64_t timestamp, double &max_wallet_val_usdt, double &max_drawdown, WalletTrace &trace)
{
    if (wallet_val_usdt > max_wallet_val_usdt)
    {
        max_wallet_val_usdt = wallet_val_usdt;
    }

    const double pc_change_with_max = (wallet_val_usdt - max_wallet_val_usdt) / max_wallet_val_usdt * 100.0;
    if (pc_change_with_max < max_drawdown)
    {
        max_drawdown = pc_change_with_max;
    }

    trace.wallet_values.push_back(wallet_val_usdt);
    trace.timestamps.push_back(timestamp);
}

ResultMetrics calculate_result_metrics(const double wallet_val_usdt, const double initial_usdt, const double max_drawdown, const TradeStats &stats)
{
    ResultMetrics metrics{};
    metrics.gain = (wallet_val_usdt - initial_usdt) / initial_usdt * 100.0;

    // A parameter set that never entered a position: win_rate would be 0/0. Report a
    // neutral, finite result instead of letting NaN reach the banner and score file.
    if (stats.nb_positions_entered == 0)
    {
        return metrics;
    }

    metrics.win_rate = double(stats.nb_profit) / double(stats.nb_positions_entered) * 100.0;

    // "Drawdown corrected": the gain needed to recover from max_drawdown. Zero drawdown
    // makes it 0, so guard the two divisions below rather than emitting +/-inf -- an inf
    // score would win the sweep outright over every genuinely-measured candidate.
    metrics.ddc = (1.0 / (1.0 + max_drawdown / 100.0) - 1.0) * 100.0;
    if (!(metrics.ddc > 0.0))
    {
        return metrics;
    }

    metrics.gain_over_ddc = metrics.gain / metrics.ddc;
    metrics.score = metrics.gain_over_ddc * metrics.win_rate;
    return metrics;
}

void populate_common_result(RUN_RESULTf &result, const ResultMetrics &metrics, const double wallet_val_usdt, const double max_drawdown, const double total_fees_paid_usdt, const TradeStats &stats, const uint max_open_trades)
{
    result.WALLET_VAL_USDT = wallet_val_usdt;
    result.gain_over_DDC = metrics.gain_over_ddc;
    result.gain_pc = metrics.gain;
    result.max_DD = max_drawdown;
    result.nb_posi_entered = stats.nb_positions_entered;
    result.win_rate = metrics.win_rate;
    result.score = metrics.score;
    result.total_fees_paid = total_fees_paid_usdt;
    result.max_open_trades = max_open_trades;
}
} // namespace trade_core
