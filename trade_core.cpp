#include "trade_core.hh"
namespace trade_core
{
ResultMetrics calculate_result_metrics(double final, double initial, double dd, const TradeStats &stats)
{
    ResultMetrics m;
    if (!(initial > 0) || !std::isfinite(final))
        return m;
    m.gain = 100 * (final / initial - 1);
    if (stats.nb_closed)
        m.win_rate = 100.0 * stats.nb_profit / stats.nb_closed;
    if (dd < 0 && dd > -100)
    {
        m.ddc = 100 * (1 / (1 + dd / 100) - 1);
        m.gain_over_ddc = m.gain / m.ddc;
        m.score = m.gain_over_ddc * m.win_rate;
    }
    return m;
}
void populate_common_result(RUN_RESULTf &r, const ResultMetrics &m, double wallet, double dd, double fees,
                            const TradeStats &s, uint limit)
{
    r.WALLET_VAL_USDT = wallet;
    r.gain_pc = m.gain;
    r.win_rate = m.win_rate;
    r.max_DD = dd;
    r.gain_over_DDC = m.gain_over_ddc;
    r.score = m.score;
    r.total_fees_paid = fees;
    r.nb_posi_entered = s.nb_positions_entered;
    r.max_open_trades = limit;
}
std::optional<double> calmar_ratio(double initial, double final, int64_t seconds, double dd)
{
    if (!(initial > 0 && final > 0 && seconds > 0 && dd < 0 && dd > -100))
        return std::nullopt;
    const double cagr = std::expm1(std::log(final / initial) * (365.25 * 86400) / double(seconds));
    const double ratio = cagr / (-dd / 100);
    return std::isfinite(ratio) ? std::optional<double>(ratio) : std::nullopt;
}
} // namespace trade_core
