#include "trade_core.hh"

namespace trade_core
{
void record_wallet_snapshot(const float wallet_val_usdt, const int timestamp, float &max_wallet_val_usdt, float &max_drawdown, WalletTrace &trace)
{
    if (wallet_val_usdt > max_wallet_val_usdt)
    {
        max_wallet_val_usdt = wallet_val_usdt;
    }

    const float pc_change_with_max = (wallet_val_usdt - max_wallet_val_usdt) / max_wallet_val_usdt * 100.0f;
    if (pc_change_with_max < max_drawdown)
    {
        max_drawdown = pc_change_with_max;
    }

    trace.wallet_values.push_back(wallet_val_usdt);
    trace.timestamps.push_back(timestamp);
}

ResultMetrics calculate_result_metrics(const float wallet_val_usdt, const float initial_usdt, const float max_drawdown, const TradeStats &stats)
{
    ResultMetrics metrics{};
    metrics.gain = (wallet_val_usdt - initial_usdt) / initial_usdt * 100.0f;
    metrics.win_rate = float(stats.nb_profit) / float(stats.nb_positions_entered) * 100.0f;
    metrics.ddc = (1.0f / (1.0f + max_drawdown / 100.0f) - 1.0f) * 100.0f;
    metrics.score = metrics.gain / metrics.ddc * metrics.win_rate;
    return metrics;
}

void populate_common_result(RUN_RESULTf &result, const ResultMetrics &metrics, const float wallet_val_usdt, const float max_drawdown, const float total_fees_paid_usdt, const TradeStats &stats, const uint max_open_trades)
{
    result.WALLET_VAL_USDT = wallet_val_usdt;
    result.gain_over_DDC = metrics.gain / metrics.ddc;
    result.gain_pc = metrics.gain;
    result.max_DD = max_drawdown;
    result.nb_posi_entered = stats.nb_positions_entered;
    result.win_rate = metrics.win_rate;
    result.score = metrics.score;
    result.total_fees_paid = total_fees_paid_usdt;
    result.max_open_trades = max_open_trades;
}
} // namespace trade_core
