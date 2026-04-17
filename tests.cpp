// Minimal self-contained unit tests for the shared backtester machinery. Each test
// calls a named function that asserts its invariants; main() runs them all and
// reports. A test fails by printing "FAIL: <name> <expr>" and setting exit=1 — no
// external framework needed.
//
// Build:  make tests
// Run:    ./tests.exe
//
// Covered:
//   - trade_core::open_spot_long / close_spot_long fee accounting
//   - trade_core::calculate_result_metrics edge cases
//   - trade_core::record_wallet_snapshot drawdown math
//   - trade_core futures long/short PnL sign
//   - tools::integer_range boundary semantics
//   - tools::float_Nvalues_range N=1 edge case
//   - tools::find_max/find_min handle all-negative inputs
//   - tools::get_funding_fee_if_any timestamp filtering
//   - custom_talib_wrapper::TALIB_EMA warmup length

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "custom_talib_wrapper.hh"
#include "tools.hh"
#include "trade_core.hh"
#include <ta-lib/ta_libc.h>

namespace
{
int g_fail_count = 0;
int g_run_count = 0;

#define REQUIRE(expr)                                                                                          \
    do                                                                                                         \
    {                                                                                                          \
        ++g_run_count;                                                                                         \
        if (!(expr))                                                                                           \
        {                                                                                                      \
            ++g_fail_count;                                                                                    \
            std::cerr << "FAIL " << __FUNCTION__ << ": " #expr " at " << __FILE__ << ":" << __LINE__ << "\n";  \
        }                                                                                                      \
    } while (0)

#define REQUIRE_NEAR(a, b, eps)                                                                                 \
    do                                                                                                          \
    {                                                                                                           \
        ++g_run_count;                                                                                          \
        const double _a = static_cast<double>(a);                                                               \
        const double _b = static_cast<double>(b);                                                               \
        if (std::fabs(_a - _b) > (eps))                                                                         \
        {                                                                                                       \
            ++g_fail_count;                                                                                     \
            std::cerr << "FAIL " << __FUNCTION__ << ": " #a "=" << _a << " vs " #b "=" << _b                    \
                      << " (tol=" << (eps) << ") at " << __FILE__ << ":" << __LINE__ << "\n";                   \
        }                                                                                                       \
    } while (0)

void test_open_close_spot_long_fee_roundtrip()
{
    // Open then close at flat price: wallet should end at initial * (1-fee)^2.
    constexpr float initial = 1000.0f;
    constexpr float price = 100.0f;
    constexpr float fee = 0.1f;
    constexpr float f = fee / 100.0f;
    trade_core::PortfolioState<1> state(initial);
    trade_core::TradeStats stats{};
    trade_core::open_spot_long(state, stats, 0, price, fee, 1);
    REQUIRE_NEAR(state.usdt_amount, 0.0f, 1e-4);
    REQUIRE(state.coin_amounts[0] > 0.0f);
    trade_core::close_spot_long(state, stats, 0, price, fee);
    REQUIRE_NEAR(state.usdt_amount, initial * (1.0f - f) * (1.0f - f), 1e-2);
    REQUIRE_NEAR(state.coin_amounts[0], 0.0f, 1e-6);
    REQUIRE(stats.nb_positions_entered == 1);
    REQUIRE(stats.nb_profit + stats.nb_loss == 1);
}

void test_open_close_spot_long_price_up()
{
    // Price doubles during the hold: wallet ≈ 2 * initial * (1-fee)^2.
    constexpr float initial = 1000.0f;
    constexpr float fee = 0.1f;
    constexpr float f = fee / 100.0f;
    trade_core::PortfolioState<1> state(initial);
    trade_core::TradeStats stats{};
    trade_core::open_spot_long(state, stats, 0, 100.0f, fee, 1);
    trade_core::close_spot_long(state, stats, 0, 200.0f, fee);
    REQUIRE_NEAR(state.usdt_amount, 2.0f * initial * (1.0f - f) * (1.0f - f), 1e-2);
    REQUIRE(stats.nb_profit == 1);
    REQUIRE(stats.nb_loss == 0);
}

void test_calculate_result_metrics()
{
    trade_core::TradeStats stats{};
    stats.nb_positions_entered = 10;
    stats.nb_profit = 4;
    const trade_core::ResultMetrics m = trade_core::calculate_result_metrics(1500.0f, 1000.0f, -20.0f, stats);
    REQUIRE_NEAR(m.gain, 50.0f, 1e-4);
    REQUIRE_NEAR(m.win_rate, 40.0f, 1e-4);
    // ddc = (1 / (1 + -20/100) - 1) * 100 = (1/0.8 - 1)*100 = 25
    REQUIRE_NEAR(m.ddc, 25.0f, 1e-3);
    // score = 50 / 25 * 40 = 80
    REQUIRE_NEAR(m.score, 80.0f, 1e-3);
}

void test_record_wallet_snapshot_drawdown()
{
    trade_core::WalletTrace trace{};
    float max_val = 1000.0f;
    float max_dd = 0.0f;

    // ATH at 1000
    trade_core::record_wallet_snapshot(1000.0f, 1, max_val, max_dd, trace);
    REQUIRE_NEAR(max_val, 1000.0f, 1e-6);
    REQUIRE_NEAR(max_dd, 0.0f, 1e-6);

    // New ATH at 2000
    trade_core::record_wallet_snapshot(2000.0f, 2, max_val, max_dd, trace);
    REQUIRE_NEAR(max_val, 2000.0f, 1e-6);
    REQUIRE_NEAR(max_dd, 0.0f, 1e-6);

    // Drawdown to 1500 → -25%
    trade_core::record_wallet_snapshot(1500.0f, 3, max_val, max_dd, trace);
    REQUIRE_NEAR(max_val, 2000.0f, 1e-6);
    REQUIRE_NEAR(max_dd, -25.0f, 1e-3);

    // Recover to 1800, dd still -25 (deeper-DD is kept)
    trade_core::record_wallet_snapshot(1800.0f, 4, max_val, max_dd, trace);
    REQUIRE_NEAR(max_dd, -25.0f, 1e-3);

    // Crash to 1000 → -50%
    trade_core::record_wallet_snapshot(1000.0f, 5, max_val, max_dd, trace);
    REQUIRE_NEAR(max_dd, -50.0f, 1e-3);

    REQUIRE(trace.wallet_values.size() == 5);
}

void test_futures_long_close_sign()
{
    constexpr float initial = 1000.0f;
    constexpr float fee = 0.1f;
    trade_core::PortfolioState<1> state(initial);
    trade_core::TradeStats stats{};
    trade_core::open_futures_long(state, stats, 0, 100.0f, fee, 1);
    trade_core::close_futures_long(state, stats, 0, 120.0f, fee);
    REQUIRE(stats.nb_profit == 1);
    REQUIRE(stats.nb_loss == 0);
    REQUIRE(state.usdt_amount > initial * 1.15f); // +20% before fees
}

void test_futures_short_close_sign()
{
    constexpr float initial = 1000.0f;
    constexpr float fee = 0.1f;
    trade_core::PortfolioState<1> state(initial);
    trade_core::TradeStats stats{};
    trade_core::open_futures_short(state, stats, 0, 100.0f, fee, 1);
    trade_core::close_futures_short(state, stats, 0, 80.0f, fee); // price down -> short profits
    REQUIRE(stats.nb_profit == 1);
    REQUIRE(stats.nb_loss == 0);
    REQUIRE(state.usdt_amount > initial * 1.15f);
}

void test_integer_range()
{
    // Three-arg half-open: [3, 11) step 2 -> {3,5,7,9}
    const std::vector<int> r1 = integer_range(3, 11, 2);
    REQUIRE(r1.size() == 4);
    REQUIRE(r1.front() == 3);
    REQUIRE(r1.back() == 9);

    // Two-arg closed: [5, 8] -> {5,6,7,8}
    const std::vector<int> r2 = integer_range(5, 8);
    REQUIRE(r2.size() == 4);
    REQUIRE(r2.front() == 5);
    REQUIRE(r2.back() == 8);

    // min == max on two-arg: {min}
    const std::vector<int> r3 = integer_range(7, 7);
    REQUIRE(r3.size() == 1);
    REQUIRE(r3[0] == 7);
}

void test_float_range_N1()
{
    const std::vector<float> r = float_Nvalues_range(1.0f, 100.0f, 1);
    REQUIRE(r.size() == 1);
    REQUIRE_NEAR(r[0], 1.0f, 1e-6);

    const std::vector<float> r2 = float_Nvalues_range(0.0f, 10.0f, 11);
    REQUIRE(r2.size() == 11);
    REQUIRE_NEAR(r2.front(), 0.0f, 1e-6);
    REQUIRE_NEAR(r2.back(), 10.0f, 1e-6);
    REQUIRE_NEAR(r2[5], 5.0f, 1e-6);
}

void test_find_max_all_negative()
{
    // Regression: pre-fix find_max on an all-negative input returned 0 (worse than any
    // element in the data). Must now return the maximum (least-negative) element.
    const std::vector<float> v{-3.0f, -1.5f, -7.0f, -2.0f};
    const float m = find_max(v);
    REQUIRE_NEAR(m, -1.5f, 1e-6);

    const std::vector<float> empty{};
    const float me = find_max(empty);
    REQUIRE(me == std::numeric_limits<float>::lowest());
}

void test_find_min_all_positive()
{
    const std::vector<float> v{3.0f, 1.5f, 7.0f, 2.0f};
    REQUIRE_NEAR(find_min(v), 1.5f, 1e-6);
}

void test_get_funding_fee_timing()
{
    fundings FR{};
    // Build a synthetic funding series at 00:00, 08:00, 16:00 UTC on 2023-06-01.
    const long base = 1685577600; // 2023-06-01 00:00:00 UTC
    for (int i = 0; i < 3; ++i)
    {
        const long ts = base + i * 8 * 3600;
        FR.timestamp.push_back(ts);
        FR.funding.push_back(0.0001f * (i + 1));
        FR.funding_by_timestamp.emplace(ts, 0.0001f * (i + 1));
    }

    REQUIRE_NEAR(get_funding_fee_if_any(FR, base), 0.0001, 1e-8);
    REQUIRE_NEAR(get_funding_fee_if_any(FR, base + 8 * 3600), 0.0002, 1e-8);
    // Miss: 1 hour past a funding timestamp
    REQUIRE_NEAR(get_funding_fee_if_any(FR, base + 3600), 0.0, 1e-8);
    // Hit the slot but no matching entry in the map
    REQUIRE_NEAR(get_funding_fee_if_any(FR, base + 3 * 8 * 3600), 0.0, 1e-8);
}

void test_talib_ema_warmup()
{
    TA_RetCode rc = TA_Initialize();
    if (rc != TA_SUCCESS)
    {
        std::cerr << "WARN: TA-Lib init failed in tests, skipping TALIB_EMA sanity.\n";
        return;
    }
    std::vector<float> series;
    series.reserve(200);
    for (int i = 0; i < 200; ++i)
    {
        series.push_back(100.0f + static_cast<float>(i));
    }
    const std::vector<float> ema = TALIB_EMA(series, 20);
    REQUIRE(ema.size() == series.size());
    // First 19 points are warmup zeros; last EMA should be close to the linear-trend tail.
    REQUIRE(ema[0] == 0.0f);
    REQUIRE(ema.back() > 200.0f);
    REQUIRE(ema.back() < series.back());
    TA_Shutdown();
}

using TestFn = void (*)();
struct NamedTest
{
    const char *name;
    TestFn fn;
};

const NamedTest ALL_TESTS[] = {
    {"open_close_spot_long_fee_roundtrip", test_open_close_spot_long_fee_roundtrip},
    {"open_close_spot_long_price_up", test_open_close_spot_long_price_up},
    {"calculate_result_metrics", test_calculate_result_metrics},
    {"record_wallet_snapshot_drawdown", test_record_wallet_snapshot_drawdown},
    {"futures_long_close_sign", test_futures_long_close_sign},
    {"futures_short_close_sign", test_futures_short_close_sign},
    {"integer_range", test_integer_range},
    {"float_range_N1", test_float_range_N1},
    {"find_max_all_negative", test_find_max_all_negative},
    {"find_min_all_positive", test_find_min_all_positive},
    {"get_funding_fee_timing", test_get_funding_fee_timing},
    {"talib_ema_warmup", test_talib_ema_warmup},
};
} // namespace

int main()
{
    for (const NamedTest &t : ALL_TESTS)
    {
        const int failures_before = g_fail_count;
        t.fn();
        const bool ok = g_fail_count == failures_before;
        std::cout << (ok ? "PASS " : "FAIL ") << t.name << "\n";
    }

    std::cout << "\n"
              << (g_run_count - g_fail_count) << " / " << g_run_count << " checks passed\n";
    if (g_fail_count != 0)
    {
        std::cout << "TESTS FAILED: " << g_fail_count << "\n";
        return 1;
    }
    std::cout << "ALL TESTS PASSED\n";
    return 0;
}
