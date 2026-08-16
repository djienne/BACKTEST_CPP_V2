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
    std::vector<float> wallet_values{};
    std::vector<int> timestamps{};
};

// All fields are guaranteed finite: calculate_result_metrics() special-cases the
// degenerate inputs (no trades, no drawdown) rather than letting 0/0 or x/0 through.
// This matters because the sweep ranks candidates with `score >`, and a single inf
// would beat every real result while a NaN would silently never be selected.
struct ResultMetrics
{
    float gain = 0.0f;
    float win_rate = 0.0f;
    float ddc = 0.0f;
    float gain_over_ddc = 0.0f;
    float score = 0.0f;
};

template <size_t N>
struct PortfolioState
{
    float usdt_amount = 0.0f;
    float max_wallet_val_usdt = 0.0f;
    float wallet_val_usdt = 0.0f;
    float total_fees_paid_usdt = 0.0f;
    float max_drawdown = 0.0f;
    uint active_positions = 0;
    std::array<float, N> price_position_open{};
    std::array<float, N> coin_amounts{};

    explicit PortfolioState(const float initial_usdt)
        : usdt_amount(initial_usdt), max_wallet_val_usdt(initial_usdt), wallet_val_usdt(initial_usdt)
    {
        price_position_open.fill(0.0f);
        coin_amounts.fill(0.0f);
    }
};

void record_wallet_snapshot(const float wallet_val_usdt, const int timestamp, float &max_wallet_val_usdt, float &max_drawdown, WalletTrace &trace);
ResultMetrics calculate_result_metrics(const float wallet_val_usdt, const float initial_usdt, const float max_drawdown, const TradeStats &stats);
void populate_common_result(RUN_RESULTf &result, const ResultMetrics &metrics, const float wallet_val_usdt, const float max_drawdown, const float total_fees_paid_usdt, const TradeStats &stats, const uint max_open_trades);

template <size_t N>
float calculate_spot_wallet_val_usdt(const PortfolioState<N> &state, const std::array<float, N> &current_prices)
{
    return state.usdt_amount + vector_product<N>(state.coin_amounts, current_prices);
}

template <size_t N>
void record_spot_snapshot(PortfolioState<N> &state, WalletTrace &trace, const std::array<float, N> &current_prices, const int timestamp)
{
    state.wallet_val_usdt = calculate_spot_wallet_val_usdt(state, current_prices);
    record_wallet_snapshot(state.wallet_val_usdt, timestamp, state.max_wallet_val_usdt, state.max_drawdown, trace);
}

template <size_t N>
void record_futures_snapshot(PortfolioState<N> &state, WalletTrace &trace, const std::array<float, N> &current_prices, const int timestamp)
{
    state.wallet_val_usdt = calculate_wallet_val_usdt<N>(state.usdt_amount, state.coin_amounts, current_prices, state.price_position_open);
    record_wallet_snapshot(state.wallet_val_usdt, timestamp, state.max_wallet_val_usdt, state.max_drawdown, trace);
}

template <size_t N>
void close_spot_long(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const float close_price, const float fee_pc)
{
    const float to_add = state.coin_amounts[pair_index] * close_price;
    state.usdt_amount += to_add;
    state.coin_amounts[pair_index] = 0.0f;
    state.active_positions--;

    const float fee = to_add * fee_pc / 100.0f;
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
void open_spot_long(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const float close_price, const float fee_pc, const uint max_open_trades)
{
    assert(state.active_positions < max_open_trades);
    state.price_position_open[pair_index] = close_price;

    const float usd_multiplier = 1.0f / float(max_open_trades - state.active_positions);
    state.coin_amounts[pair_index] = state.usdt_amount * usd_multiplier / close_price;
    state.usdt_amount -= state.usdt_amount * usd_multiplier;

    const float fee = state.coin_amounts[pair_index] * fee_pc / 100.0f;
    state.coin_amounts[pair_index] -= fee;
    state.total_fees_paid_usdt += fee * close_price;

    state.active_positions++;
    stats.nb_positions_entered++;
}

template <size_t N>
void apply_funding_fee(PortfolioState<N> &state, const uint pair_index, const float close_price, const float funding_fee)
{
    if (funding_fee != 0.0f && state.coin_amounts[pair_index] != 0.0f)
    {
        const float fee = state.coin_amounts[pair_index] * close_price * funding_fee;
        state.total_fees_paid_usdt += fee;
        state.usdt_amount -= fee;
    }
}

template <size_t N>
void close_futures_long(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const float close_price, const float fee_pc)
{
    const float to_add = state.coin_amounts[pair_index] * close_price;
    state.usdt_amount += to_add;

    const float fee = to_add * fee_pc / 100.0f;
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
    state.coin_amounts[pair_index] = 0.0f;
}

template <size_t N>
void close_futures_short(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const float close_price, const float fee_pc)
{
    const float to_add = std::abs(state.coin_amounts[pair_index]) * (2.0f * state.price_position_open[pair_index] - close_price);
    state.usdt_amount += to_add;

    const float fee = to_add * fee_pc / 100.0f;
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
    state.coin_amounts[pair_index] = 0.0f;
}

template <size_t N>
void open_futures_long(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const float close_price, const float fee_pc, const uint max_open_trades)
{
    assert(state.active_positions < max_open_trades);
    const float usd_multiplier = 1.0f / float(max_open_trades - state.active_positions);
    state.coin_amounts[pair_index] = state.usdt_amount * usd_multiplier / close_price;
    state.usdt_amount -= state.usdt_amount * usd_multiplier;

    const float fee = std::abs(state.coin_amounts[pair_index] * close_price * fee_pc / 100.0f);
    state.usdt_amount -= fee;
    state.total_fees_paid_usdt += fee;

    state.price_position_open[pair_index] = close_price;
    state.active_positions++;
    stats.nb_positions_entered++;
}

template <size_t N>
void open_futures_short(PortfolioState<N> &state, TradeStats &stats, const uint pair_index, const float close_price, const float fee_pc, const uint max_open_trades)
{
    assert(state.active_positions < max_open_trades);
    const float usd_multiplier = 1.0f / float(max_open_trades - state.active_positions);
    state.coin_amounts[pair_index] = -1.0f * state.usdt_amount * usd_multiplier / close_price;
    state.usdt_amount -= state.usdt_amount * usd_multiplier;

    const float fee = std::abs(state.coin_amounts[pair_index] * close_price * fee_pc / 100.0f);
    state.usdt_amount -= fee;
    state.total_fees_paid_usdt += fee;

    state.price_position_open[pair_index] = close_price;
    state.active_positions++;
    stats.nb_positions_entered++;
}
} // namespace trade_core
