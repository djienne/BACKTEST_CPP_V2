#pragma once
#include "tools.hh"
#include <cassert>
#include <stdexcept>

namespace trade_core
{
struct TradeStats
{
    uint nb_profit = 0, nb_loss = 0, nb_positions_entered = 0, nb_closed = 0;
};
struct ResultMetrics
{
    double gain = 0, win_rate = 0, ddc = 0, gain_over_ddc = 0, score = 0;
};
template <size_t N> struct PortfolioState
{
    double usdt_amount, max_wallet_val_usdt, wallet_val_usdt;
    double total_fees_paid_usdt = 0, net_funding = 0, max_drawdown = 0;
    uint active_positions = 0, nb_pairs;
    std::array<double, N> coin_amounts{}, price_position_open{}, entry_cost{}, position_funding{};
    explicit PortfolioState(double initial, uint pairs = N)
        : usdt_amount(initial), max_wallet_val_usdt(initial), wallet_val_usdt(initial), nb_pairs(pairs)
    {
        if (!(initial > 0) || !std::isfinite(initial) || pairs > N)
            throw std::runtime_error("Invalid portfolio");
    }
};
ResultMetrics calculate_result_metrics(double, double, double, const TradeStats &);
void populate_common_result(RUN_RESULTf &, const ResultMetrics &, double, double, double, const TradeStats &, uint);
std::optional<double> calmar_ratio(double initial, double final, int64_t elapsed_seconds, double max_drawdown);

template <size_t N>
double calculate_spot_wallet_val_usdt(const PortfolioState<N> &s, const std::array<float, N> &prices)
{
    double value = s.usdt_amount;
    for (uint i = 0; i < s.nb_pairs; ++i)
        value += s.coin_amounts[i] * prices[i];
    return value;
}
template <size_t N> double equity(const PortfolioState<N> &s, const std::array<float, N> &prices, bool futures)
{
    if (!futures)
        return calculate_spot_wallet_val_usdt(s, prices);
    return calculate_wallet_val_usdt<N>(s.usdt_amount, s.coin_amounts, prices, s.price_position_open, s.nb_pairs);
}
template <size_t N>
void open_position(PortfolioState<N> &s, TradeStats &stats, uint pair, double price, double fee_pc, uint limit,
                   int side, bool futures)
{
    if (pair >= s.nb_pairs || s.coin_amounts[pair] != 0 || s.active_positions >= limit || !(price > 0) ||
        !std::isfinite(price) || !(fee_pc >= 0 && fee_pc < 100) || !(s.usdt_amount > 0) || (side != 1 && side != -1) ||
        (!futures && side < 0))
        throw std::runtime_error("Invalid position entry");
    const double budget = s.usdt_amount / double(limit - s.active_positions);
    const double f = fee_pc / 100;
    const double notional = futures ? budget / (1 + f) : budget;
    const double fee = notional * f;
    s.coin_amounts[pair] = side * (futures ? notional : notional - fee) / price;
    s.usdt_amount -= budget;
    s.entry_cost[pair] = budget;
    s.position_funding[pair] = 0;
    s.total_fees_paid_usdt += fee;
    s.price_position_open[pair] = price;
    ++s.active_positions;
    ++stats.nb_positions_entered;
}
template <size_t N>
double close_position(PortfolioState<N> &s, TradeStats &stats, uint pair, double price, double fee_pc, bool futures)
{
    if (pair >= s.nb_pairs || s.coin_amounts[pair] == 0 || !(price > 0) || !std::isfinite(price))
        throw std::runtime_error("Invalid position exit");
    const double q = s.coin_amounts[pair];
    const double proceeds = q > 0 ? q * price : std::abs(q) * (2 * s.price_position_open[pair] - price);
    if (!futures && q < 0)
        throw std::runtime_error("Short in spot account");
    const double fee = std::abs(q) * price * fee_pc / 100;
    const double pnl = proceeds - fee - s.entry_cost[pair] - s.position_funding[pair];
    s.usdt_amount += proceeds - fee;
    s.total_fees_paid_usdt += fee;
    // Breakevens are closed trades, but not winners.
    if (pnl > 0)
        ++stats.nb_profit;
    else if (pnl < 0)
        ++stats.nb_loss;
    ++stats.nb_closed;
    --s.active_positions;
    s.coin_amounts[pair] = 0;
    return pnl;
}
template <size_t N> void open_spot_long(PortfolioState<N> &s, TradeStats &t, uint p, double c, double f, uint n)
{
    open_position(s, t, p, c, f, n, 1, false);
}
template <size_t N> void open_futures_long(PortfolioState<N> &s, TradeStats &t, uint p, double c, double f, uint n)
{
    open_position(s, t, p, c, f, n, 1, true);
}
template <size_t N> void open_futures_short(PortfolioState<N> &s, TradeStats &t, uint p, double c, double f, uint n)
{
    open_position(s, t, p, c, f, n, -1, true);
}
template <size_t N> void close_spot_long(PortfolioState<N> &s, TradeStats &t, uint p, double c, double f)
{
    close_position(s, t, p, c, f, false);
}
template <size_t N> void close_futures_long(PortfolioState<N> &s, TradeStats &t, uint p, double c, double f)
{
    close_position(s, t, p, c, f, true);
}
template <size_t N> void close_futures_short(PortfolioState<N> &s, TradeStats &t, uint p, double c, double f)
{
    close_position(s, t, p, c, f, true);
}
template <size_t N> void apply_funding_fee(PortfolioState<N> &s, uint p, double mark, double rate)
{
    const double amount = s.coin_amounts[p] * mark * rate;
    s.usdt_amount -= amount;
    s.position_funding[p] += amount;
    s.net_funding += amount;
}
struct Window
{
    size_t begin = 0, end = 0;
}; // execution-bar indices, end exclusive
struct Intent
{
    int entry = 0; // +1 long, -1 short, 0 no entry
    bool exit_long = false, exit_short = false;
    double stop_distance = 0, target_distance = 0, target_fraction = 0;
    bool trailing = false;
    int64_t max_hold_seconds = 0;
};

// Common OHLC execution path. Signals are observed at completed signal-bar closes.
// The coarser signal clock and funding-compatible execution clock stay separate.
template <class Signal>
RUN_RESULTf simulate(const MarketData &data, Window window, bool futures, uint max_open, Signal signal,
                     bool trace = false, double initial = 1000, double fee_pc = 0.1)
{
    constexpr size_t N = backtest_config::MAX_PAIRS;
    RUN_RESULTf result;
    const uint pairs = static_cast<uint>(data.execution.size());
    if (!pairs || pairs > N || !max_open || max_open > pairs || window.begin >= window.end ||
        window.end > data.execution[0].nb)
        throw std::runtime_error("Invalid evaluation window or position limit");
    PortfolioState<N> account(initial, pairs);
    TradeStats stats;
    std::array<Intent, N> pending{};
    std::array<double, N> stops{}, targets{}, trail_distance{};
    std::array<int64_t, N> deadline{};
    std::array<size_t, N> funding_index{};
    std::array<float, N> prices{};
    auto close = [&](uint p, double price, int64_t time, const std::string &why)
    {
        const double qty = account.coin_amounts[p], fees = account.total_fees_paid_usdt;
        const double pnl = close_position(account, stats, p, price, fee_pc, futures);
        if (trace)
            result.fills.push_back({time, p, why, price, qty, account.total_fees_paid_usdt - fees, pnl});
    };
    const auto first_time = data.execution[0].timestamp[window.begin];
    // A fresh evaluation can use a pre-window completed signal, never a carried position.
    auto set_signal = [&](int64_t close_time)
    {
        const auto &times = data.signal[0].timestamp;
        auto it = std::upper_bound(times.begin(), times.end(), close_time - data.signal_seconds);
        if (it == times.begin())
            return;
        const size_t i = static_cast<size_t>(it - times.begin() - 1);
        if (times[i] + data.signal_seconds != close_time)
            return;
        for (uint p = 0; p < pairs; ++p)
            pending[p] = signal(p, i);
    };
    set_signal(first_time);
    if (trace)
    {
        result.equity.push_back(initial);
        result.equity_times.push_back(first_time);
    }
    for (size_t i = window.begin; i < window.end; ++i)
    {
        const int64_t time = data.execution[0].timestamp[i], close_time = time + data.execution_seconds;
        const bool last = i + 1 == window.end;
        std::array<bool, N> protected_exit{};
        // Funding belongs to holders just before the timestamp, before open-time orders.
        for (uint p = 0; p < pairs; ++p)
        {
            const auto &fund = data.funding[p];
            while (funding_index[p] < fund.size() && fund[funding_index[p]].timestamp < time)
                ++funding_index[p];
            while (funding_index[p] < fund.size() && fund[funding_index[p]].timestamp == time &&
                   !fund[funding_index[p]].millisecond)
            {
                const auto &e = fund[funding_index[p]++];
                apply_funding_fee(account, p, e.mark_price, e.rate);
            }
        }
        // Complete all exits before any entry allocation, in configured coin order.
        for (uint p = 0; p < pairs; ++p)
        {
            const double q = account.coin_amounts[p], op = data.execution[p].open[i];
            if (!q)
                continue;
            const bool stop = stops[p] > 0 && (q > 0 ? op <= stops[p] : op >= stops[p]);
            const bool target = targets[p] > 0 && (q > 0 ? op >= targets[p] : op <= targets[p]);
            if (stop || target)
            {
                close(p, op, time, stop ? "gap_stop" : "gap_target");
                protected_exit[p] = true;
            }
            else if (deadline[p] && time >= deadline[p])
                close(p, op, time, "timeout");
            else if (q > 0 ? pending[p].exit_long : pending[p].exit_short)
                close(p, op, time, "signal_exit");
        }
        for (uint p = 0; p < pairs; ++p)
        {
            const auto &order = pending[p];
            if (!last && !protected_exit[p] && !account.coin_amounts[p] && order.entry &&
                account.active_positions < max_open && account.usdt_amount > 0)
            {
                const double price = data.execution[p].open[i], fees = account.total_fees_paid_usdt;
                open_position(account, stats, p, price, fee_pc, max_open, order.entry, futures);
                stops[p] = order.stop_distance > 0 ? price - order.entry * order.stop_distance : 0;
                double distance = order.target_distance;
                if (order.target_fraction > 0)
                    distance = distance > 0 ? std::min(distance, price * order.target_fraction)
                                            : price * order.target_fraction;
                targets[p] = distance > 0 ? price + order.entry * distance : 0;
                trail_distance[p] = order.trailing ? order.stop_distance : 0;
                deadline[p] = order.max_hold_seconds ? time + order.max_hold_seconds : 0;
                if (trace)
                    result.fills.push_back(
                        {time, p, "entry", price, account.coin_amounts[p], account.total_fees_paid_usdt - fees, 0});
            }
            pending[p] = {}; // a signal creates one next-open order, not a standing order
        }
        // Published settlements sometimes occur milliseconds after the opening.
        // Opening orders then precede funding. OHLC cannot order a stop within
        // those milliseconds: settlement precedes unresolved intrabar fills here.
        for (uint p = 0; p < pairs; ++p)
        {
            const auto &fund = data.funding[p];
            while (funding_index[p] < fund.size() && fund[funding_index[p]].timestamp == time)
            {
                const auto &e = fund[funding_index[p]++];
                apply_funding_fee(account, p, e.mark_price, e.rate);
            }
        }
        for (uint p = 0; p < pairs; ++p)
        {
            const auto &bar = data.execution[p];
            const double q = account.coin_amounts[p];
            if (q)
            {
                const bool stop = stops[p] > 0 && (q > 0 ? bar.low[i] <= stops[p] : bar.high[i] >= stops[p]);
                const bool target = targets[p] > 0 && (q > 0 ? bar.high[i] >= targets[p] : bar.low[i] <= targets[p]);
                if (stop && target)
                    ++result.ambiguous_bars;
                if (stop || target)
                    close(p, stop ? stops[p] : targets[p], close_time, stop ? "stop" : "target");
                else if (last)
                    close(p, bar.close[i], close_time, "final_exit");
                else if (trail_distance[p] > 0)
                {
                    const double next = bar.close[i] - (q > 0 ? 1 : -1) * trail_distance[p];
                    stops[p] = q > 0 ? std::max(stops[p], next) : std::min(stops[p], next);
                }
            }
            prices[p] = bar.close[i];
        }
        account.wallet_val_usdt = equity(account, prices, futures);
        if (!(account.wallet_val_usdt > 0) || !std::isfinite(account.wallet_val_usdt))
        {
            result.valid = false;
            result.invalid_reason = "non-positive or non-finite equity";
            break;
        }
        account.max_wallet_val_usdt = std::max(account.max_wallet_val_usdt, account.wallet_val_usdt);
        account.max_drawdown =
            std::min(account.max_drawdown, 100 * (account.wallet_val_usdt / account.max_wallet_val_usdt - 1));
        if (trace)
        {
            result.equity.push_back(account.wallet_val_usdt);
            result.equity_times.push_back(close_time);
        }
        if (!last)
            set_signal(close_time);
    }
    const auto metrics = calculate_result_metrics(account.wallet_val_usdt, initial, account.max_drawdown, stats);
    populate_common_result(result, metrics, account.wallet_val_usdt, account.max_drawdown, account.total_fees_paid_usdt,
                           stats, max_open);
    result.net_funding = account.net_funding;
    result.calmar_ratio = calmar_ratio(
        initial, account.wallet_val_usdt,
        data.execution[0].timestamp[window.end - 1] + data.execution_seconds - first_time, account.max_drawdown);
    if (!std::isfinite(result.score))
    {
        result.valid = false;
        result.invalid_reason = "non-finite score";
    }
    return result;
}
} // namespace trade_core
