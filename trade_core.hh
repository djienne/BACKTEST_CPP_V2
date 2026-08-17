#pragma once

#include <cassert>

#include "tools.hh"

namespace trade_core
{
struct TradeStats
{
    uint nb_profit = 0;
    uint nb_loss = 0;
    uint nb_positions_entered = 0;
};

struct WalletTrace
{
    std::vector<double> wallet_values{};
    std::vector<int64_t> timestamps{};
};

// All fields are guaranteed finite: calculate_result_metrics() special-cases the
// degenerate inputs (no trades, no drawdown) rather than letting 0/0 or x/0 through.
// This matters because the sweep ranks candidates with `score >`, and a single inf
// would beat every real result while a NaN would silently never be selected.
// All fields are guaranteed finite: calculate_result_metrics() special-cases the
// degenerate inputs (no trades, no drawdown) rather than letting 0/0 or x/0 through.
struct ResultMetrics
{
    double gain = 0.0;
    double win_rate = 0.0;
    double ddc = 0.0;
    double gain_over_ddc = 0.0;
    double score = 0.0;
};

// Wallet, fees and position sizes are double while prices and indicators stay float.
// These are accumulators: a run makes thousands of fee subtractions and position
// round-trips, and float's ~7 significant digits were being spent on a balance that
// routinely reaches five figures -- in the very numbers the sweep ranks on.
template <size_t N>
struct PortfolioState
{
    double usdt_amount = 0.0;
    double max_wallet_val_usdt = 0.0;
    double wallet_val_usdt = 0.0;
    double total_fees_paid_usdt = 0.0;
    double max_drawdown = 0.0;
    uint active_positions = 0;
    // Pairs actually traded. N is a compile-time ceiling (backtest_config::MAX_PAIRS) so
    // the arrays stay fixed-size and the hot loop allocation-free; the real count comes
    // from the config at runtime. Slots [nb_pairs, N) are never touched.
    uint nb_pairs = N;
    std::array<double, N> price_position_open{};
    std::array<double, N> coin_amounts{};

    explicit PortfolioState(const double initial_usdt, const uint active_pairs = N)
        : usdt_amount(initial_usdt), max_wallet_val_usdt(initial_usdt), wallet_val_usdt(initial_usdt),
          nb_pairs(active_pairs)
    {
        assert(active_pairs <= N);
        price_position_open.fill(0.0);
        coin_amounts.fill(0.0);
    }
};

void record_wallet_snapshot(const double wallet_val_usdt, const int64_t timestamp, double &max_wallet_val_usdt, double &max_drawdown, WalletTrace &trace);
ResultMetrics calculate_result_metrics(const double wallet_val_usdt, const double initial_usdt, const double max_drawdown, const TradeStats &stats);
void populate_common_result(RUN_RESULTf &result, const ResultMetrics &metrics, const double wallet_val_usdt, const double max_drawdown, const double total_fees_paid_usdt, const TradeStats &stats, const uint max_open_trades);

template <size_t N>
double calculate_spot_wallet_val_usdt(const PortfolioState<N> &state, const std::array<float, N> &current_prices)
{
    return state.usdt_amount + vector_product<N>(state.coin_amounts, current_prices, state.nb_pairs);
}

template <size_t N>
void record_spot_snapshot(PortfolioState<N> &state, WalletTrace &trace, const std::array<float, N> &current_prices, const int64_t timestamp)
{
    state.wallet_val_usdt = calculate_spot_wallet_val_usdt(state, current_prices);
    record_wallet_snapshot(state.wallet_val_usdt, timestamp, state.max_wallet_val_usdt, state.max_drawdown, trace);
}

template <size_t N>
void record_futures_snapshot(PortfolioState<N> &state, WalletTrace &trace, const std::array<float, N> &current_prices, const int64_t timestamp)
{
    state.wallet_val_usdt = calculate_wallet_val_usdt<N>(state.usdt_amount, state.coin_amounts, current_prices, state.price_position_open, state.nb_pairs);
    record_wallet_snapshot(state.wallet_val_usdt, timestamp, state.max_wallet_val_usdt, state.max_drawdown, trace);
}

template <size_t N>
void close_spot_long(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const double close_price, const double fee_pc)
{
    const double to_add = state.coin_amounts[pair_index] * close_price;
    state.usdt_amount += to_add;
    state.coin_amounts[pair_index] = 0.0;
    state.active_positions--;

    const double fee = to_add * fee_pc / 100.0;
    state.usdt_amount -= fee;
    state.total_fees_paid_usdt += fee;

    if (close_price >= state.price_position_open[pair_index])
    {
        stats.nb_profit++;
    }
    else
    {
        stats.nb_loss++;
    }
}

template <size_t N>
void open_spot_long(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const double close_price, const double fee_pc, const uint max_open_trades)
{
    assert(state.active_positions < max_open_trades);
    state.price_position_open[pair_index] = close_price;

    // Split the remaining cash across the slots still free, so a book that fills up
    // ends with equally-sized positions.
    const double usd_multiplier = 1.0 / double(max_open_trades - state.active_positions);
    state.coin_amounts[pair_index] = state.usdt_amount * usd_multiplier / close_price;
    state.usdt_amount -= state.usdt_amount * usd_multiplier;

    // Spot entry fee is taken in the purchased coin, not in USDT.
    const double fee = state.coin_amounts[pair_index] * fee_pc / 100.0;
    state.coin_amounts[pair_index] -= fee;
    state.total_fees_paid_usdt += fee * close_price;

    state.active_positions++;
    stats.nb_positions_entered++;
}

template <size_t N>
void apply_funding_fee(PortfolioState<N> &state, const uint pair_index, const double close_price, const double funding_fee)
{
    if (funding_fee != 0.0 && state.coin_amounts[pair_index] != 0.0)
    {
        const double fee = state.coin_amounts[pair_index] * close_price * funding_fee;
        state.total_fees_paid_usdt += fee;
        state.usdt_amount -= fee;
    }
}

// Futures model, applying to all four helpers below: leverage is 1 and there is no
// liquidation. A short is valued as its mirrored long, (2 * entry - price), which is
// exact at leverage 1 but has no floor -- a short that runs far enough against the book
// can drive the wallet negative and the sweep will keep trading it. Read futures results
// as an un-liquidated upper bound, not as an exchange-faithful simulation.
template <size_t N>
void close_futures_long(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const double close_price, const double fee_pc)
{
    const double to_add = state.coin_amounts[pair_index] * close_price;
    state.usdt_amount += to_add;

    const double fee = to_add * fee_pc / 100.0;
    state.usdt_amount -= std::abs(fee);
    state.total_fees_paid_usdt += std::abs(fee);

    if (close_price >= state.price_position_open[pair_index])
    {
        stats.nb_profit++;
    }
    else
    {
        stats.nb_loss++;
    }

    state.active_positions--;
    state.coin_amounts[pair_index] = 0.0;
}

template <size_t N>
void close_futures_short(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const double close_price, const double fee_pc)
{
    const double to_add = std::abs(state.coin_amounts[pair_index]) * (2.0 * state.price_position_open[pair_index] - close_price);
    state.usdt_amount += to_add;

    const double fee = to_add * fee_pc / 100.0;
    state.usdt_amount -= std::abs(fee);
    state.total_fees_paid_usdt += std::abs(fee);

    if (close_price <= state.price_position_open[pair_index])
    {
        stats.nb_profit++;
    }
    else
    {
        stats.nb_loss++;
    }

    state.active_positions--;
    state.coin_amounts[pair_index] = 0.0;
}

template <size_t N>
void open_futures_long(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const double close_price, const double fee_pc, const uint max_open_trades)
{
    assert(state.active_positions < max_open_trades);
    const double usd_multiplier = 1.0 / double(max_open_trades - state.active_positions);
    state.coin_amounts[pair_index] = state.usdt_amount * usd_multiplier / close_price;
    state.usdt_amount -= state.usdt_amount * usd_multiplier;

    const double fee = std::abs(state.coin_amounts[pair_index] * close_price * fee_pc / 100.0);
    state.usdt_amount -= fee;
    state.total_fees_paid_usdt += fee;

    state.price_position_open[pair_index] = close_price;
    state.active_positions++;
    stats.nb_positions_entered++;
}

template <size_t N>
void open_futures_short(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const double close_price, const double fee_pc, const uint max_open_trades)
{
    assert(state.active_positions < max_open_trades);
    const double usd_multiplier = 1.0 / double(max_open_trades - state.active_positions);
    // Negative coin amount marks the position as short.
    state.coin_amounts[pair_index] = -1.0 * state.usdt_amount * usd_multiplier / close_price;
    state.usdt_amount -= state.usdt_amount * usd_multiplier;

    const double fee = std::abs(state.coin_amounts[pair_index] * close_price * fee_pc / 100.0);
    state.usdt_amount -= fee;
    state.total_fees_paid_usdt += fee;

    state.price_position_open[pair_index] = close_price;
    state.active_positions++;
    stats.nb_positions_entered++;
}
} // namespace trade_core
