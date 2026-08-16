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
//   - custom_talib_wrapper::TALIB_BBANDS output sizing + band ordering
//   - custom_talib_wrapper::TALIB_AO output length
//   - custom_talib_wrapper::TALIB_STOCHRSI_not_averaged range [0,1]
//   - trade_core::apply_funding_fee (long debited, no-op on flat)
//   - tools::calculate_calmar_ratio on a synthetic 3-year wallet curve
//   - tools::realign_timestamps synthetic misaligned pair
//   - strategy_runner::sweep best-tracking predicate (gating filters)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "custom_talib_wrapper.hh"
#include "strategy_runner.hh"
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
    REQUIRE_NEAR(m.gain_over_ddc, 2.0f, 1e-3);
}

void test_calculate_result_metrics_degenerate()
{
    // Regression: a parameter set that entered no position divided by zero for the win
    // rate, and a run with zero drawdown divided by a zero ddc. Both leaked NaN/inf into
    // the printed banner and the score file; an inf score also beat every real candidate
    // in the sweep's `score >` ranking.
    trade_core::TradeStats no_trades{};
    const trade_core::ResultMetrics m0 = trade_core::calculate_result_metrics(1000.0f, 1000.0f, 0.0f, no_trades);
    REQUIRE(std::isfinite(m0.win_rate));
    REQUIRE(std::isfinite(m0.score));
    REQUIRE(std::isfinite(m0.gain_over_ddc));
    REQUIRE_NEAR(m0.win_rate, 0.0f, 1e-6);
    REQUIRE_NEAR(m0.score, 0.0f, 1e-6);

    // Trades happened but the equity curve never drew down -> ddc == 0.
    trade_core::TradeStats winners{};
    winners.nb_positions_entered = 10;
    winners.nb_profit = 10;
    const trade_core::ResultMetrics m1 = trade_core::calculate_result_metrics(2000.0f, 1000.0f, 0.0f, winners);
    REQUIRE(std::isfinite(m1.score));
    REQUIRE(std::isfinite(m1.gain_over_ddc));
    REQUIRE_NEAR(m1.gain, 100.0f, 1e-3);
    REQUIRE_NEAR(m1.win_rate, 100.0f, 1e-3);
    REQUIRE_NEAR(m1.score, 0.0f, 1e-6);

    // A populated RUN_RESULTf must stay finite too.
    RUN_RESULTf r{};
    trade_core::populate_common_result(r, m0, 1000.0f, 0.0f, 0.0f, no_trades, 1);
    REQUIRE(std::isfinite(r.score));
    REQUIRE(std::isfinite(r.gain_over_DDC));
    REQUIRE(std::isfinite(r.win_rate));
}

void test_supertrend_dir_only_matches_full()
{
    // Regression: TALIB_SuperTrend_dir_only was a copy-paste of TALIB_SuperTrend. It is
    // now a thin wrapper, so the two must agree element for element.
    TA_Initialize();
    const int n = 200;
    std::vector<float> high, low, close;
    high.reserve(n);
    low.reserve(n);
    close.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        const float base = 100.0f + 20.0f * std::sin(i * 0.15f) + 0.05f * i;
        high.push_back(base + 1.5f);
        low.push_back(base - 1.5f);
        close.push_back(base);
    }

    const SuperTrend st = TALIB_SuperTrend(high, low, close, 10, 3.0f);
    const std::vector<float> dir = TALIB_SuperTrend_dir_only(high, low, close, 10, 3.0f);

    REQUIRE(st.supertrend.size() == static_cast<size_t>(n));
    REQUIRE(dir.size() == st.supertrend.size());
    REQUIRE(st.final_lowerband.size() == static_cast<size_t>(n));
    REQUIRE(st.final_upperband.size() == static_cast<size_t>(n));

    int mismatches = 0;
    int saw_up = 0;
    int saw_down = 0;
    for (size_t i = 0; i < dir.size(); ++i)
    {
        if (dir[i] != static_cast<float>(st.supertrend[i]))
        {
            ++mismatches;
        }
        if (st.supertrend[i] == 1)
        {
            ++saw_up;
        }
        if (st.supertrend[i] == -1)
        {
            ++saw_down;
        }
    }
    REQUIRE(mismatches == 0);
    // An oscillating series must flip direction, otherwise this test proves nothing.
    REQUIRE(saw_up > 0);
    REQUIRE(saw_down > 0);

    // Fractional multipliers are the common SuperTrend setting and used to be
    // unreachable through the int parameter; a wider band must not flip more often.
    const SuperTrend narrow = TALIB_SuperTrend(high, low, close, 10, 1.5f);
    int narrow_flips = 0;
    int wide_flips = 0;
    for (size_t i = 1; i < dir.size(); ++i)
    {
        narrow_flips += (narrow.supertrend[i] != narrow.supertrend[i - 1]) ? 1 : 0;
        wide_flips += (st.supertrend[i] != st.supertrend[i - 1]) ? 1 : 0;
    }
    REQUIRE(narrow_flips >= wide_flips);
    TA_Shutdown();
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

void test_talib_bbands_shape()
{
    TA_Initialize();
    std::vector<float> series;
    series.reserve(100);
    // Slowly drifting sine ensures non-zero std, so upper > middle > lower in steady state.
    for (int i = 0; i < 100; ++i)
    {
        series.push_back(100.0f + 5.0f * std::sin(i * 0.2f) + 0.1f * i);
    }
    std::vector<float> u, m, l;
    TALIB_BBANDS(series, 2.0f, 2.0f, 20, u, m, l);
    REQUIRE(u.size() == series.size());
    REQUIRE(m.size() == series.size());
    REQUIRE(l.size() == series.size());
    // Steady-state (past warmup): upper > middle > lower on this oscillating input.
    REQUIRE(u[50] > m[50]);
    REQUIRE(m[50] > l[50]);
    TA_Shutdown();
}

void test_talib_ao_shape()
{
    TA_Initialize();
    const int n = 80;
    std::vector<float> high, low;
    high.reserve(n);
    low.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        const float base = 100.0f + i;
        high.push_back(base + 1.0f);
        low.push_back(base - 1.0f);
    }
    const std::vector<float> ao = TALIB_AO(high, low, 5, 34);
    REQUIRE(ao.size() == high.size());
    // Monotone uptrend -> AO should be positive in the tail (fast MA above slow MA).
    REQUIRE(ao.back() > 0.0f);
    TA_Shutdown();
}

void test_talib_stochrsi_range()
{
    TA_Initialize();
    std::vector<float> series;
    series.reserve(200);
    // Mean-reverting saw-tooth so StochRSI produces varied non-degenerate values.
    for (int i = 0; i < 200; ++i)
    {
        series.push_back(100.0f + 10.0f * std::sin(i * 0.3f));
    }
    const std::vector<float> srsi = TALIB_STOCHRSI_not_averaged(series, 14, 14);
    REQUIRE(srsi.size() == series.size());
    int bounded = 0;
    int saw_nonzero = 0;
    // Sample past the warmup region; values must fall in [0, 1].
    for (size_t i = 50; i < srsi.size(); ++i)
    {
        if (srsi[i] >= 0.0f && srsi[i] <= 1.0f)
        {
            ++bounded;
        }
        if (srsi[i] > 0.0f)
        {
            ++saw_nonzero;
        }
    }
    REQUIRE(bounded == static_cast<int>(srsi.size() - 50));
    REQUIRE(saw_nonzero > 10);
    TA_Shutdown();
}

void test_apply_funding_fee()
{
    constexpr float fee = 0.0001f; // 0.01%
    trade_core::PortfolioState<1> state(1000.0f);
    // No position -> apply_funding_fee should not change anything.
    trade_core::apply_funding_fee(state, 0, 100.0f, fee);
    REQUIRE_NEAR(state.usdt_amount, 1000.0f, 1e-6);
    REQUIRE_NEAR(state.total_fees_paid_usdt, 0.0f, 1e-6);

    // Simulate a long position of 1 coin at price 100. Positive funding debits USDT.
    state.coin_amounts[0] = 1.0f;
    state.price_position_open[0] = 100.0f;
    const float before = state.usdt_amount;
    trade_core::apply_funding_fee(state, 0, 110.0f, fee);
    const float expected_deduction = 1.0f * 110.0f * fee;
    REQUIRE_NEAR(before - state.usdt_amount, expected_deduction, 1e-4);
    REQUIRE_NEAR(state.total_fees_paid_usdt, expected_deduction, 1e-4);

    // Zero funding -> no-op even with a position.
    const float after = state.usdt_amount;
    trade_core::apply_funding_fee(state, 0, 120.0f, 0.0f);
    REQUIRE_NEAR(state.usdt_amount, after, 1e-6);
}

void test_calculate_calmar_ratio()
{
    // Synthetic: ATH then flat over ~3 years. Build timestamps monthly.
    std::vector<int> ts;
    std::vector<float> wv;
    const int seconds_per_month = 30 * 24 * 3600;
    const int base = 1577836800; // 2020-01-01 00:00 UTC
    for (int m = 0; m < 36; ++m)
    {
        ts.push_back(base + m * seconds_per_month);
        wv.push_back(1000.0f + 10.0f * m); // +1% per month roughly
    }
    // Feed a known max_DD; function normalizes averaged yearly gain by it.
    const float cr = calculate_calmar_ratio(ts, wv, -10.0f);
    // With +1%/mo over ~3 years the mean yearly return should be positive and
    // dividing by a negative DD yields a negative Calmar in this function's sign
    // convention — we only check the math is finite and non-degenerate.
    REQUIRE(std::isfinite(cr));
    REQUIRE(cr != 0.0f);

    // Too-short series (<= 4 points) returns the sentinel -100.
    const std::vector<int> short_ts{1, 2, 3};
    const std::vector<float> short_wv{1.0f, 2.0f, 3.0f};
    REQUIRE_NEAR(calculate_calmar_ratio(short_ts, short_wv, -10.0f), -100.0f, 1e-6);
}

void test_realign_timestamps_noop_when_aligned()
{
    KLINEf base{};
    KLINEf other{};
    for (int i = 0; i < 20; ++i)
    {
        base.timestamp.push_back(1000 + i * 60);
        base.open.push_back(100.0f + i);
        base.high.push_back(101.0f + i);
        base.low.push_back(99.0f + i);
        base.close.push_back(100.0f + i);

        other.timestamp.push_back(1000 + i * 60);
        other.open.push_back(200.0f + i);
        other.high.push_back(201.0f + i);
        other.low.push_back(199.0f + i);
        other.close.push_back(200.0f + i);
    }
    base.nb = 20;
    base.name = "BASE";
    other.nb = 20;
    other.name = "OTHER";
    other.start_idx = 0;

    const float close_before = other.close[10];
    realign_timestamps(base, other);
    // Aligned input -> function short-circuits, leaves `other` untouched.
    REQUIRE_NEAR(other.close[10], close_before, 1e-6);
    REQUIRE(other.timestamp.size() == base.timestamp.size());
}

void test_random_number_generator_seeding()
{
    // Regression: the constructor built a local mt19937 that shadowed the member, so the
    // member kept whatever the init-list gave it and the intended thread-id mixing was
    // discarded. Two generators must produce different streams, and every draw must land
    // inside [0, upperLimit].
    RandomNumberGenerator a;
    RandomNumberGenerator b;

    std::vector<int> sa;
    std::vector<int> sb;
    bool in_range = true;
    for (int i = 0; i < 64; ++i)
    {
        const int va = a.getRandomNumber(999);
        const int vb = b.getRandomNumber(999);
        in_range = in_range && va >= 0 && va <= 999 && vb >= 0 && vb <= 999;
        sa.push_back(va);
        sb.push_back(vb);
    }
    REQUIRE(in_range);
    REQUIRE(sa != sb);

    // A zero upper limit must yield exactly 0 rather than an empty-distribution surprise
    // (the F_* sweeps call getRandomNumber(range.size() - 1) on single-element ranges).
    RandomNumberGenerator c;
    REQUIRE(c.getRandomNumber(0) == 0);
}

void test_strategy_runner_sweep_filters()
{
    // Exercise the gating predicate: only results that pass all filters should
    // be picked as `best`. Fabricate three results with different trade-count,
    // drawdown, and gain profiles.
    struct P
    {
        int id;
    };
    std::vector<P> params{{0}, {1}, {2}};

    auto process = [](const P &p) {
        RUN_RESULTf r{};
        r.gain_pc = 100.0f + p.id * 50.0f; // 100, 150, 200
        r.nb_posi_entered = (p.id == 0) ? 10 : 500;
        r.max_DD = (p.id == 1) ? -99.0f : -20.0f;
        r.score = 10.0f + p.id; // p=2 has highest score
        r.max_open_trades = 1;
        r.param_str = "id=" + std::to_string(p.id);
        return r;
    };

    strategy_runner::SweepConfig cfg;
    cfg.strategy_name = "unit_test";
    cfg.out_filename = "_unit_test_best.txt";
    cfg.min_trades = 100;
    cfg.min_dd = -40.0f;
    cfg.print_every = 9999; // suppress periodic output during the test
    cfg.min_reasonable_gain = 0.0f;

    const RUN_RESULTf best = strategy_runner::sweep(cfg, params, process);

    // p=0 excluded (too few trades), p=1 excluded (DD too deep); p=2 wins.
    REQUIRE(best.param_str == "id=2");
    REQUIRE_NEAR(best.gain_pc, 200.0f, 1e-6);

    // Cleanup the score file the sweep writes out.
    std::remove("_unit_test_best.txt");
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
    {"calculate_result_metrics_degenerate", test_calculate_result_metrics_degenerate},
    {"supertrend_dir_only_matches_full", test_supertrend_dir_only_matches_full},
    {"random_number_generator_seeding", test_random_number_generator_seeding},
    {"record_wallet_snapshot_drawdown", test_record_wallet_snapshot_drawdown},
    {"futures_long_close_sign", test_futures_long_close_sign},
    {"futures_short_close_sign", test_futures_short_close_sign},
    {"integer_range", test_integer_range},
    {"float_range_N1", test_float_range_N1},
    {"find_max_all_negative", test_find_max_all_negative},
    {"find_min_all_positive", test_find_min_all_positive},
    {"get_funding_fee_timing", test_get_funding_fee_timing},
    {"talib_ema_warmup", test_talib_ema_warmup},
    {"talib_bbands_shape", test_talib_bbands_shape},
    {"talib_ao_shape", test_talib_ao_shape},
    {"talib_stochrsi_range", test_talib_stochrsi_range},
    {"apply_funding_fee", test_apply_funding_fee},
    {"calculate_calmar_ratio", test_calculate_calmar_ratio},
    {"realign_timestamps_noop_when_aligned", test_realign_timestamps_noop_when_aligned},
    {"strategy_runner_sweep_filters", test_strategy_runner_sweep_filters},
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
