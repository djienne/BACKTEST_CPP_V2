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
    // The tolerances here are 1e-9 rather than the 1e-2 the float wallet needed --
    // that gap is the point of accumulating money in double.
    constexpr double initial = 1000.0;
    constexpr double price = 100.0;
    constexpr double fee = 0.1;
    constexpr double f = fee / 100.0;
    trade_core::PortfolioState<1> state(initial);
    trade_core::TradeStats stats{};
    trade_core::open_spot_long(state, stats, 0, price, fee, 1);
    REQUIRE_NEAR(state.usdt_amount, 0.0, 1e-9);
    REQUIRE(state.coin_amounts[0] > 0.0);
    trade_core::close_spot_long(state, stats, 0, price, fee);
    REQUIRE_NEAR(state.usdt_amount, initial * (1.0 - f) * (1.0 - f), 1e-9);
    REQUIRE_NEAR(state.coin_amounts[0], 0.0, 1e-12);
    REQUIRE(stats.nb_positions_entered == 1);
    REQUIRE(stats.nb_profit + stats.nb_loss == 1);
}

void test_open_close_spot_long_price_up()
{
    // Price doubles during the hold: wallet = 2 * initial * (1-fee)^2.
    constexpr double initial = 1000.0;
    constexpr double fee = 0.1;
    constexpr double f = fee / 100.0;
    trade_core::PortfolioState<1> state(initial);
    trade_core::TradeStats stats{};
    trade_core::open_spot_long(state, stats, 0, 100.0, fee, 1);
    trade_core::close_spot_long(state, stats, 0, 200.0, fee);
    REQUIRE_NEAR(state.usdt_amount, 2.0 * initial * (1.0 - f) * (1.0 - f), 1e-9);
    REQUIRE(stats.nb_profit == 1);
    REQUIRE(stats.nb_loss == 0);
}

void test_wallet_precision_over_many_roundtrips()
{
    // Regression for the float wallet: 5000 flat-price round-trips at a 0.1% fee. The
    // closed form is initial * (1-f)^(2*n), and with a float accumulator the running
    // balance drifted off it well beyond any sane tolerance. Double must track it.
    constexpr double initial = 1000.0;
    constexpr double fee = 0.1;
    constexpr double f = fee / 100.0;
    constexpr int n = 5000;

    trade_core::PortfolioState<1> state(initial);
    trade_core::TradeStats stats{};
    for (int i = 0; i < n; ++i)
    {
        trade_core::open_spot_long(state, stats, 0, 100.0, fee, 1);
        trade_core::close_spot_long(state, stats, 0, 100.0, fee);
    }

    const double expected = initial * std::pow(1.0 - f, 2 * n);
    REQUIRE(stats.nb_positions_entered == n);
    // Relative error must stay near machine epsilon rather than float's ~1e-7 per op.
    REQUIRE(std::fabs(state.usdt_amount - expected) / expected < 1e-12);
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
    double max_val = 1000.0;
    double max_dd = 0.0;

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
    constexpr double fee = 0.0001; // 0.01%
    trade_core::PortfolioState<1> state(1000.0);
    // No position -> apply_funding_fee should not change anything.
    trade_core::apply_funding_fee(state, 0, 100.0, fee);
    REQUIRE_NEAR(state.usdt_amount, 1000.0, 1e-12);
    REQUIRE_NEAR(state.total_fees_paid_usdt, 0.0, 1e-12);

    // Simulate a long position of 1 coin at price 100. Positive funding debits USDT.
    state.coin_amounts[0] = 1.0;
    state.price_position_open[0] = 100.0;
    // double, matching PortfolioState: capturing the balance in a float here would
    // reintroduce exactly the rounding this change removes.
    const double before = state.usdt_amount;
    trade_core::apply_funding_fee(state, 0, 110.0, fee);
    const double expected_deduction = 1.0 * 110.0 * fee;
    REQUIRE_NEAR(before - state.usdt_amount, expected_deduction, 1e-9);
    REQUIRE_NEAR(state.total_fees_paid_usdt, expected_deduction, 1e-9);

    // Zero funding -> no-op even with a position.
    const double after = state.usdt_amount;
    trade_core::apply_funding_fee(state, 0, 120.0, 0.0);
    REQUIRE_NEAR(state.usdt_amount, after, 1e-12);
}

void test_calculate_calmar_ratio()
{
    // Synthetic: ATH then flat over ~3 years. Build timestamps monthly.
    std::vector<int64_t> ts;
    std::vector<double> wv;
    const int64_t seconds_per_month = 30 * 24 * 3600;
    const int64_t base = 1577836800; // 2020-01-01 00:00 UTC
    for (int m = 0; m < 36; ++m)
    {
        ts.push_back(base + m * seconds_per_month);
        wv.push_back(1000.0 + 10.0 * m); // +1% per month roughly
    }
    // Callers pass ResultMetrics::ddc (positive) as the normalizer.
    const double cr = calculate_calmar_ratio(ts, wv, 10.0);
    REQUIRE(std::isfinite(cr));
    REQUIRE(cr > 0.0);

    // Too-short series (<= 4 points) returns the sentinel -100.
    const std::vector<int64_t> short_ts{1, 2, 3};
    const std::vector<double> short_wv{1.0, 2.0, 3.0};
    REQUIRE_NEAR(calculate_calmar_ratio(short_ts, short_wv, 10.0), -100.0, 1e-6);

    // A non-positive normalizer (zero drawdown) must not divide by zero.
    REQUIRE_NEAR(calculate_calmar_ratio(ts, wv, 0.0), -100.0, 1e-6);
    REQUIRE(std::isfinite(calculate_calmar_ratio_monthly(ts, wv, 0.0)));
}

void test_convert_to_unix_timestamp_is_utc()
{
    // Regression: this used std::mktime, which reads the tm as *local* time. Since
    // read_input_data_f uses the result as a row cutoff, the number of rows loaded from
    // a futures file depended on the machine's timezone -- CI on UTC loaded one more row
    // than a UTC+2 laptop, which is how it surfaced.
    //
    // 2023-06-18 00:00:00 UTC.
    REQUIRE(convertToUnixTimestamp("2023-06-18") == 1687046400);
    // 1970-01-01 is the epoch itself, so any local-time offset shows up as non-zero.
    REQUIRE(convertToUnixTimestamp("1970-01-01") == 0);
    // Past the signed-32-bit boundary: 2040-01-01 00:00:00 UTC.
    REQUIRE(convertToUnixTimestamp("2040-01-01") == 2208988800);
    // Unparseable input keeps its -1 sentinel.
    REQUIRE(convertToUnixTimestamp("not-a-date") == -1);

    // Round-trips through the UTC calendar helpers.
    const int64_t ts = convertToUnixTimestamp("2023-06-18");
    REQUIRE(get_year_from_timestamp(ts) == 2023);
    REQUIRE(get_month_from_timestamp(ts) == 6);
    REQUIRE(get_day_from_timestamp(ts) == 18);
    REQUIRE(get_hour_from_timestamp(ts) == 0);

    // getCurrentDateMinusTwoDays must parse back and land two days before today (UTC),
    // which is only true if both functions agree on the timezone.
    const std::string two_days_ago = getCurrentDateMinusTwoDays();
    REQUIRE(two_days_ago.size() == 10);
    const int64_t parsed = convertToUnixTimestamp(two_days_ago);
    REQUIRE(parsed > 0);
    const int64_t now = static_cast<int64_t>(std::time(nullptr));
    const int64_t age_days = (now - parsed) / 86400;
    REQUIRE(age_days >= 2);
    REQUIRE(age_days <= 3); // 2 whole days plus however far into today we are
}

void test_generate_range_int()
{
    // Regression: the step was an integer division, so the range never reached vmax.
    // generateRange_int(3, 600, 300) used to yield 3..302 -- half the intended sweep
    // space, silently -- and N == 1 divided by zero.
    const std::vector<int> r = generateRange_int(3, 600, 300);
    REQUIRE(r.front() == 3);
    REQUIRE(r.back() == 600);
    REQUIRE(r.size() == 300);

    const std::vector<int> r2 = generateRange_int(5, 400, 120);
    REQUIRE(r2.front() == 5);
    REQUIRE(r2.back() == 400);

    // Monotonically non-decreasing and strictly increasing after dedup.
    bool increasing = true;
    for (size_t i = 1; i < r.size(); ++i)
    {
        increasing = increasing && r[i] > r[i - 1];
    }
    REQUIRE(increasing);

    // N == 1 is the old divide-by-zero.
    const std::vector<int> r3 = generateRange_int(7, 99, 1);
    REQUIRE(r3.size() == 1);
    REQUIRE(r3[0] == 7);

    // More points requested than the interval holds: deduplicated, endpoints kept.
    const std::vector<int> r4 = generateRange_int(1, 3, 10);
    REQUIRE(r4.front() == 1);
    REQUIRE(r4.back() == 3);
    REQUIRE(r4.size() == 3);
}

void test_utc_timestamp_helpers()
{
    // Regression: these used localtime(), so results depended on the host timezone and
    // raced on the shared static tm from the F_* worker threads. Exchange klines are
    // UTC, so these must read UTC regardless of TZ.
    // 2023-06-18 16:00:00 UTC.
    const int64_t ts = 1687104000;
    REQUIRE(get_year_from_timestamp(ts) == 2023);
    REQUIRE(get_month_from_timestamp(ts) == 6);
    REQUIRE(get_day_from_timestamp(ts) == 18);
    REQUIRE(get_hour_from_timestamp(ts) == 16);

    // A timestamp beyond the 2038 signed-32-bit boundary must still resolve, which the
    // previous `int` parameters made impossible.
    // 2040-01-01 00:00:00 UTC.
    const int64_t ts_2040 = 2208988800;
    REQUIRE(get_year_from_timestamp(ts_2040) == 2040);
    REQUIRE(get_month_from_timestamp(ts_2040) == 1);
    REQUIRE(get_day_from_timestamp(ts_2040) == 1);
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

void test_strategy_param_structs_are_fully_initialized()
{
    // Regression for the bug that stopped 3EMA_SRSI_ATR opening any position at all:
    // EMA3_params has eight members and its initializer supplied seven, so the
    // max-open-trades value landed in SRSIU and max_open_trades defaulted to 0 --
    // making `active_positions < MAX_OPEN_TRADES` false on every bar.
    //
    // Verified here rather than by running the strategy: its sweep is 83.5 million
    // parameter sets over 5m data, and it only prints a result block at the end.
    const EMA3_params p{3, 50, 200, 5.0f, 5.0f, 0.2f, 0.8f, 4};
    REQUIRE(p.ema1 == 3);
    REQUIRE(p.ema2 == 50);
    REQUIRE(p.ema3 == 200);
    REQUIRE_NEAR(p.up, 5.0f, 1e-6);
    REQUIRE_NEAR(p.down, 5.0f, 1e-6);
    REQUIRE_NEAR(p.SRSIL, 0.2f, 1e-6);
    REQUIRE_NEAR(p.SRSIU, 0.8f, 1e-6);
    // The one that mattered: a zero here means the strategy cannot trade.
    REQUIRE(p.max_open_trades == 4);
    REQUIRE(p.max_open_trades != 0);

    // A zero cap makes open_spot_long unreachable, which is exactly what the broken
    // initializer produced. Show the guard the strategies use rejects it.
    trade_core::PortfolioState<2> state(1000.0, 2);
    const uint broken_cap = 0;
    REQUIRE(!(state.active_positions < broken_cap)); // no position can ever open
    const uint good_cap = 2;
    REQUIRE(state.active_positions < good_cap);

    // The other param structs must round-trip their members too.
    const BigWill_params b{5, 34, 20, 200, 3};
    REQUIRE(b.AO_fast == 5);
    REQUIRE(b.max_open_trades == 3);
    const trix_params t{100, 9, 21, 2};
    REQUIRE(t.max_open_trades == 2);
    const BBTREND_params bb{200, 20, 2.0f, 5};
    REQUIRE(bb.max_open_trades == 5);
    const ST_EMA_ATR_params st{100, 3.0f, 3.0f, 6};
    REQUIRE(st.max_open_trades == 6);
    const SR_params sr{10, 100, 4};
    REQUIRE(sr.max_open_trades == 4);
}

void test_first_tradable_index()
{
    // Takes the worst warmup among the indicators in use...
    REQUIRE(strategy_runner::first_tradable_index({10, 250, 33}) == 250);
    // ...but never starts before the pair's own data does...
    REQUIRE(strategy_runner::first_tradable_index({10, 250, 33}, 600) == 600);
    REQUIRE(strategy_runner::first_tradable_index({10, 250, 33}, 100) == 250);
    // ...and leaves room for strategies that read backwards from the current bar.
    REQUIRE(strategy_runner::first_tradable_index({10, 250}, 100, 2) == 252);
    REQUIRE(strategy_runner::first_tradable_index({}, 0, 0) == 0);

    // The point of it: a real warmup drives the start index, so raising an indicator
    // period cannot quietly outgrow a hard-coded constant.
    TA_Initialize();
    std::vector<float> series;
    for (int i = 0; i < 900; ++i)
    {
        series.push_back(100.0f + static_cast<float>(i));
    }
    size_t warm_short = 0;
    size_t warm_long = 0;
    TALIB_EMA(series, 50, &warm_short);
    TALIB_EMA(series, 800, &warm_long);
    const uint begin = strategy_runner::first_tradable_index({warm_short, warm_long}, 600);
    // The 800-period EMA warms up past the old hard-coded 600.
    REQUIRE(begin == static_cast<uint>(warm_long));
    REQUIRE(begin > 600);
    TA_Shutdown();
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

// ---------------------------------------------------------------------------
//  Indicator library
//
//  Each test pins output length, warmup position, and at least one value with a
//  known analytic answer -- shape-only assertions would pass on a wrapper plumbed
//  to the wrong TA-Lib function.
// ---------------------------------------------------------------------------

// A deterministic oscillating-with-drift OHLCV series, enough bars to clear any
// warmup used below.
struct Ohlcv
{
    std::vector<float> open, high, low, close, volume;
};

Ohlcv make_series(const int n = 300)
{
    Ohlcv s;
    for (int i = 0; i < n; ++i)
    {
        const float base = 100.0f + 15.0f * std::sin(i * 0.13f) + 0.08f * i;
        s.open.push_back(base - 0.4f);
        s.close.push_back(base);
        s.high.push_back(base + 1.2f);
        s.low.push_back(base - 1.2f);
        s.volume.push_back(1000.0f + 250.0f * std::cos(i * 0.21f));
    }
    return s;
}

// A ramp has an exact closed form for the window-average families.
std::vector<float> make_ramp(const int n = 100)
{
    std::vector<float> v;
    for (int i = 0; i < n; ++i)
    {
        v.push_back(static_cast<float>(i));
    }
    return v;
}

void test_moving_averages()
{
    TA_Initialize();
    const std::vector<float> ramp = make_ramp(100);

    // On a unit ramp, SMA(n) at index i is the mean of i-n+1..i = i - (n-1)/2.
    size_t warm_sma = 0;
    const std::vector<float> sma = TALIB_SMA(ramp, 10, &warm_sma);
    REQUIRE(sma.size() == ramp.size());
    REQUIRE(warm_sma == 9);
    REQUIRE_NEAR(sma[9], 4.5f, 1e-3);
    REQUIRE_NEAR(sma[50], 45.5f, 1e-3);
    REQUIRE(sma[8] == 0.0f); // warmup is zero-padded

    // WMA weights linearly, so on a ramp it sits closer to the newest value:
    // WMA(n) at i = i - (n-1)/3.
    size_t warm_wma = 0;
    const std::vector<float> wma = TALIB_WMA(ramp, 10, &warm_wma);
    REQUIRE(warm_wma == 9);
    REQUIRE_NEAR(wma[50], 50.0f - 3.0f, 1e-3);

    // On a pure linear trend every EMA-family average is a lag of the input; DEMA and
    // TEMA reduce that lag, so they must sit strictly closer to the current value.
    const std::vector<float> ema = TALIB_EMA(ramp, 10);
    const std::vector<float> dema = TALIB_DEMA(ramp, 10);
    const std::vector<float> tema = TALIB_TEMA(ramp, 10);
    const float current = ramp[90];
    REQUIRE(std::fabs(dema[90] - current) < std::fabs(ema[90] - current));
    REQUIRE(std::fabs(tema[90] - current) < std::fabs(ema[90] - current));

    // HMA is composed here rather than taken from TA-Lib; it must reduce lag too and
    // report a warmup that actually covers its zero-padded head.
    size_t warm_hma = 0;
    const std::vector<float> hma = TALIB_HMA(ramp, 16, &warm_hma);
    REQUIRE(hma.size() == ramp.size());
    REQUIRE(warm_hma > 0);
    REQUIRE(hma[warm_hma - 1] == 0.0f);
    REQUIRE(std::fabs(hma[90] - current) < std::fabs(ema[90] - current));

    size_t warm_kama = 0;
    const std::vector<float> kama = TALIB_KAMA(ramp, 10, &warm_kama);
    REQUIRE(kama.size() == ramp.size());
    REQUIRE(kama[90] > 0.0f);
    TA_Shutdown();
}

void test_macd()
{
    TA_Initialize();
    const Ohlcv s = make_series();
    const MACDResult m = TALIB_MACD(s.close, 12, 26, 9);

    REQUIRE(m.macd.size() == s.close.size());
    REQUIRE(m.signal.size() == s.close.size());
    REQUIRE(m.histogram.size() == s.close.size());
    REQUIRE(m.warmup > 0);
    REQUIRE(m.macd[m.warmup - 1] == 0.0f);

    // The defining identity: histogram == macd - signal, past warmup.
    bool identity_holds = true;
    bool saw_both_signs = false;
    int pos = 0, neg = 0;
    for (size_t i = m.warmup; i < m.macd.size(); ++i)
    {
        identity_holds = identity_holds && std::fabs(m.histogram[i] - (m.macd[i] - m.signal[i])) < 1e-3;
        pos += (m.histogram[i] > 0.0f) ? 1 : 0;
        neg += (m.histogram[i] < 0.0f) ? 1 : 0;
    }
    saw_both_signs = pos > 0 && neg > 0;
    REQUIRE(identity_holds);
    REQUIRE(saw_both_signs); // an oscillating input must cross zero
    TA_Shutdown();
}

void test_stoch_and_aroon()
{
    TA_Initialize();
    const Ohlcv s = make_series();

    const StochResult st = TALIB_STOCH(s.high, s.low, s.close, 14, 3, 3);
    REQUIRE(st.k.size() == s.close.size());
    REQUIRE(st.d.size() == s.close.size());
    REQUIRE(st.warmup > 0);
    // %K and %D are percentages, so they must stay within [0, 100].
    bool bounded = true;
    for (size_t i = st.warmup; i < st.k.size(); ++i)
    {
        bounded = bounded && st.k[i] >= -1e-3f && st.k[i] <= 100.001f && st.d[i] >= -1e-3f && st.d[i] <= 100.001f;
    }
    REQUIRE(bounded);

    const AroonResult ar = TALIB_AROON(s.high, s.low, 14);
    REQUIRE(ar.up.size() == s.close.size());
    REQUIRE(ar.down.size() == s.close.size());
    bool aroon_bounded = true;
    for (size_t i = ar.warmup; i < ar.up.size(); ++i)
    {
        aroon_bounded = aroon_bounded && ar.up[i] >= -1e-3f && ar.up[i] <= 100.001f &&
                        ar.down[i] >= -1e-3f && ar.down[i] <= 100.001f;
    }
    REQUIRE(aroon_bounded);
    TA_Shutdown();
}

void test_momentum_family()
{
    TA_Initialize();
    const std::vector<float> ramp = make_ramp(100);

    // MOM(n) on a unit ramp is exactly n.
    size_t warm_mom = 0;
    const std::vector<float> mom = TALIB_MOM(ramp, 10, &warm_mom);
    REQUIRE(warm_mom == 10);
    REQUIRE_NEAR(mom[50], 10.0f, 1e-3);

    // ROC(n) at i = (v[i] - v[i-n]) / v[i-n] * 100 = 10/40*100 at i=50.
    const std::vector<float> roc = TALIB_ROC(ramp, 10);
    REQUIRE_NEAR(roc[50], 25.0f, 1e-3);

    const Ohlcv s = make_series();
    size_t warm_cci = 0;
    const std::vector<float> cci = TALIB_CCI(s.high, s.low, s.close, 20, &warm_cci);
    REQUIRE(cci.size() == s.close.size());
    REQUIRE(warm_cci == 19);

    size_t warm_u = 0;
    const std::vector<float> ult = TALIB_ULTOSC(s.high, s.low, s.close, 7, 14, 28, &warm_u);
    REQUIRE(ult.size() == s.close.size());
    bool bounded = true;
    for (size_t i = warm_u; i < ult.size(); ++i)
    {
        bounded = bounded && ult[i] >= -1e-3f && ult[i] <= 100.001f;
    }
    REQUIRE(bounded);
    TA_Shutdown();
}

void test_trend_strength_family()
{
    TA_Initialize();
    const Ohlcv s = make_series();

    const DirectionalResult dmi = TALIB_DMI(s.high, s.low, s.close, 14);
    REQUIRE(dmi.adx.size() == s.close.size());
    REQUIRE(dmi.plus_di.size() == s.close.size());
    REQUIRE(dmi.minus_di.size() == s.close.size());
    REQUIRE(dmi.warmup > 0);
    bool bounded = true;
    for (size_t i = dmi.warmup; i < dmi.adx.size(); ++i)
    {
        bounded = bounded && dmi.adx[i] >= -1e-3f && dmi.adx[i] <= 100.001f;
    }
    REQUIRE(bounded);

    // A strictly rising market must show +DI above -DI.
    std::vector<float> up_h, up_l, up_c;
    for (int i = 0; i < 120; ++i)
    {
        up_c.push_back(100.0f + 2.0f * i);
        up_h.push_back(101.0f + 2.0f * i);
        up_l.push_back(99.0f + 2.0f * i);
    }
    const DirectionalResult up = TALIB_DMI(up_h, up_l, up_c, 14);
    REQUIRE(up.plus_di[100] > up.minus_di[100]);

    // Parabolic SAR trails an uptrend from below.
    size_t warm_sar = 0;
    const std::vector<float> sar = TALIB_SAR(up_h, up_l, 0.02, 0.2, &warm_sar);
    REQUIRE(sar.size() == up_c.size());
    REQUIRE(sar[100] < up_c[100]);
    TA_Shutdown();
}

void test_volatility_family()
{
    TA_Initialize();
    const Ohlcv s = make_series();

    size_t warm_natr = 0;
    const std::vector<float> natr = TALIB_NATR(s.high, s.low, s.close, 14, &warm_natr);
    REQUIRE(natr.size() == s.close.size());
    REQUIRE(natr[100] > 0.0f);

    // TRANGE on this series is a constant 2.4 (high-low), since gaps never dominate.
    size_t warm_tr = 0;
    const std::vector<float> tr = TALIB_TRANGE(s.high, s.low, s.close, &warm_tr);
    REQUIRE(tr.size() == s.close.size());
    REQUIRE(tr[100] >= 2.4f - 1e-3f);

    // StdDev of a constant series is 0.
    const std::vector<float> flat(80, 42.0f);
    const std::vector<float> sd = TALIB_STDDEV(flat, 20, 1.0);
    REQUIRE_NEAR(sd[50], 0.0f, 1e-4);

    // Keltner: bands straddle the EMA centreline, upper above lower.
    const BandsResult kc = KELTNER_CHANNELS(s.high, s.low, s.close, 20, 10, 2.0f);
    REQUIRE(kc.upper.size() == s.close.size());
    REQUIRE(kc.warmup > 0);
    REQUIRE(kc.upper[100] > kc.middle[100]);
    REQUIRE(kc.middle[100] > kc.lower[100]);

    // Donchian: upper is the rolling high, lower the rolling low, middle the midpoint,
    // and price must lie inside the envelope.
    const BandsResult dc = DONCHIAN_CHANNELS(s.high, s.low, 20);
    REQUIRE(dc.upper.size() == s.close.size());
    REQUIRE_NEAR(dc.middle[100], 0.5f * (dc.upper[100] + dc.lower[100]), 1e-3);
    REQUIRE(dc.upper[100] >= s.high[100]);
    REQUIRE(dc.lower[100] <= s.low[100]);
    TA_Shutdown();
}

void test_volume_family()
{
    TA_Initialize();

    // OBV has an exact hand-checkable definition: add volume on an up close, subtract
    // on a down close, carry on unchanged.
    const std::vector<float> close{10.0f, 11.0f, 10.5f, 10.5f, 12.0f};
    const std::vector<float> vol{100.0f, 200.0f, 300.0f, 400.0f, 500.0f};
    const std::vector<float> obv = TALIB_OBV(close, vol);
    REQUIRE(obv.size() == close.size());
    REQUIRE_NEAR(obv[0], 100.0f, 1e-3);  // seed
    REQUIRE_NEAR(obv[1], 300.0f, 1e-3);  // up   -> +200
    REQUIRE_NEAR(obv[2], 0.0f, 1e-3);    // down -> -300
    REQUIRE_NEAR(obv[3], 0.0f, 1e-3);    // flat -> unchanged
    REQUIRE_NEAR(obv[4], 500.0f, 1e-3);  // up   -> +500

    const Ohlcv s = make_series();
    size_t warm_mfi = 0;
    const std::vector<float> mfi = TALIB_MFI(s.high, s.low, s.close, s.volume, 14, &warm_mfi);
    REQUIRE(mfi.size() == s.close.size());
    bool bounded = true;
    for (size_t i = warm_mfi; i < mfi.size(); ++i)
    {
        bounded = bounded && mfi[i] >= -1e-3f && mfi[i] <= 100.001f;
    }
    REQUIRE(bounded);

    const std::vector<float> ad = TALIB_AD(s.high, s.low, s.close, s.volume);
    REQUIRE(ad.size() == s.close.size());
    const std::vector<float> adosc = TALIB_ADOSC(s.high, s.low, s.close, s.volume, 3, 10);
    REQUIRE(adosc.size() == s.close.size());

    // Constant price and volume make rolling VWAP exactly the typical price.
    const std::vector<float> ch(60, 50.0f), hh(60, 52.0f), lh(60, 48.0f), vv(60, 10.0f);
    size_t warm_vwap = 0;
    const std::vector<float> vwap = VWAP_ROLLING(hh, lh, ch, vv, 20, &warm_vwap);
    REQUIRE(warm_vwap == 19);
    REQUIRE_NEAR(vwap[40], 50.0f, 1e-3); // hlc3 of (52,48,50) == 50
    REQUIRE(vwap[18] == 0.0f);           // warmup zero-padded

    // Relative volume is 1.0 when volume equals its own average.
    size_t warm_rv = 0;
    const std::vector<float> rv = RELATIVE_VOLUME(vv, 20, &warm_rv);
    REQUIRE_NEAR(rv[40], 1.0f, 1e-4);
    // Double the volume on the last bar -> above 1.
    std::vector<float> spiky = vv;
    spiky[50] = 100.0f;
    const std::vector<float> rv2 = RELATIVE_VOLUME(spiky, 20, nullptr);
    REQUIRE(rv2[50] > 1.0f);
    TA_Shutdown();
}

void test_price_transforms_and_heikin_ashi()
{
    const std::vector<float> o{10.0f, 12.0f}, h{14.0f, 16.0f}, l{8.0f, 9.0f}, c{12.0f, 15.0f};

    const std::vector<float> hl2 = PRICE_HL2(h, l);
    REQUIRE_NEAR(hl2[0], 11.0f, 1e-6);
    const std::vector<float> hlc3 = PRICE_HLC3(h, l, c);
    REQUIRE_NEAR(hlc3[0], (14.0f + 8.0f + 12.0f) / 3.0f, 1e-5);
    const std::vector<float> ohlc4 = PRICE_OHLC4(o, h, l, c);
    REQUIRE_NEAR(ohlc4[0], (10.0f + 14.0f + 8.0f + 12.0f) / 4.0f, 1e-6);

    const HeikinAshi ha = HEIKIN_ASHI(o, h, l, c);
    REQUIRE(ha.close.size() == 2);
    // HA close is the OHLC average; HA open seeds from (open+close)/2 on bar 0 and is
    // the running average of the previous HA candle afterwards.
    REQUIRE_NEAR(ha.close[0], 11.0f, 1e-5);
    REQUIRE_NEAR(ha.open[0], 11.0f, 1e-5);
    REQUIRE_NEAR(ha.close[1], (12.0f + 16.0f + 9.0f + 15.0f) / 4.0f, 1e-5);
    REQUIRE_NEAR(ha.open[1], 0.5f * (ha.open[0] + ha.close[0]), 1e-5);
    // HA high/low must envelope the HA body.
    REQUIRE(ha.high[1] >= std::max(ha.open[1], ha.close[1]));
    REQUIRE(ha.low[1] <= std::min(ha.open[1], ha.close[1]));
    REQUIRE(ha.warmup == 1);
}

void test_resample_timeframe()
{
    // 24 five-minute candles starting exactly on an hour -> 2 one-hour candles.
    KLINEf in{};
    const int64_t base = 1687104000; // 2023-06-18 16:00:00 UTC, an exact hour boundary
    for (int i = 0; i < 24; ++i)
    {
        in.timestamp.push_back(base + i * 300);
        in.open.push_back(100.0f + i);
        in.close.push_back(100.5f + i);
        in.high.push_back(101.0f + i);
        in.low.push_back(99.0f + i);
        in.volume.push_back(10.0f);
    }
    in.nb = 24;
    in.name = "TEST";

    const Resampled r = RESAMPLE_TIMEFRAME(in, 12, 5, 60);
    const KLINEf &out = r.kline;
    REQUIRE(r.ltf_offset == 0);
    REQUIRE(out.nb == 2);
    REQUIRE(out.timestamp[0] == base);
    REQUIRE(out.timestamp[1] == base + 3600);
    // Open of the group, close of the last bar, extremes over the whole group,
    // volume summed.
    REQUIRE_NEAR(out.open[0], 100.0f, 1e-5);
    REQUIRE_NEAR(out.close[0], 100.5f + 11.0f, 1e-5);
    REQUIRE_NEAR(out.high[0], 101.0f + 11.0f, 1e-5);
    REQUIRE_NEAR(out.low[0], 99.0f, 1e-5);
    REQUIRE_NEAR(out.volume[0], 120.0f, 1e-4);
    REQUIRE_NEAR(out.open[1], 112.0f, 1e-5);
}

void test_resample_timeframe_off_boundary_start()
{
    // Regression for the real bundled data: BTC 5m futures starts at 16:55, so slicing
    // from index 0 built "hourly" candles spanning 16:55 to 17:50. The resampler must
    // skip forward to the first true hour boundary and report that offset.
    KLINEf in{};
    const int64_t hour = 1687104000;   // an exact hour
    const int64_t start = hour - 300;  // 5 minutes earlier: 16:55
    for (int i = 0; i < 26; ++i)
    {
        in.timestamp.push_back(start + i * 300);
        in.open.push_back(100.0f + i);
        in.close.push_back(100.5f + i);
        in.high.push_back(101.0f + i);
        in.low.push_back(99.0f + i);
        in.volume.push_back(10.0f);
    }
    in.nb = 26;
    in.name = "OFFSET";

    const Resampled r = RESAMPLE_TIMEFRAME(in, 12, 5, 60);
    // Bar 0 is 16:55; bar 1 is 17:00 and starts the first whole hour.
    REQUIRE(r.ltf_offset == 1);
    REQUIRE(r.kline.nb == 2);
    REQUIRE(r.kline.timestamp[0] == hour);
    REQUIRE(r.kline.timestamp[1] == hour + 3600);
    // The aggregated candle must be built from bars 1..12, not 0..11.
    REQUIRE_NEAR(r.kline.open[0], 101.0f, 1e-5);
    REQUIRE_NEAR(r.kline.close[0], 100.5f + 12.0f, 1e-5);

    // Projecting back must respect the same offset.
    const std::vector<float> proj = PROJECT_HTF_TO_LTF({7.0f, 8.0f}, 12, in.close.size(), r.ltf_offset, -1.0f);
    REQUIRE(proj.size() == 26);
    // Hour 0 spans bars 1..12 and closes at the end of bar 12, so bars 0..12 predate
    // any completed hour and carry the fill.
    REQUIRE(proj[0] == -1.0f);
    REQUIRE(proj[12] == -1.0f);
    // Hour 0's value is visible from bar 13 through bar 24; hour 1 closes at bar 24 and
    // becomes visible at bar 25.
    REQUIRE(proj[13] == 7.0f);
    REQUIRE(proj[24] == 7.0f);
    REQUIRE(proj[25] == 8.0f);
}

void test_project_htf_to_ltf_has_no_lookahead()
{
    // Three higher-timeframe values projected onto 36 lower-timeframe bars.
    const std::vector<float> htf{10.0f, 20.0f, 30.0f};
    const std::vector<float> ltf = PROJECT_HTF_TO_LTF(htf, 12, 36, 0, -777.0f);

    REQUIRE(ltf.size() == 36);

    // The critical property: group 0's value must NOT be visible during group 0 --
    // that candle has not closed yet. Bars 0..11 carry the fill sentinel.
    bool head_is_fill = true;
    for (size_t i = 0; i < 12; ++i)
    {
        head_is_fill = head_is_fill && ltf[i] == -777.0f;
    }
    REQUIRE(head_is_fill);

    // Group 0's value becomes visible across group 1, and group 1's across group 2.
    bool block1 = true, block2 = true;
    for (size_t i = 12; i < 24; ++i)
    {
        block1 = block1 && ltf[i] == 10.0f;
    }
    for (size_t i = 24; i < 36; ++i)
    {
        block2 = block2 && ltf[i] == 20.0f;
    }
    REQUIRE(block1);
    REQUIRE(block2);
    // htf[2] would only appear at bar 36, which is past the end -- never leaked.
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
    {"wallet_precision_over_many_roundtrips", test_wallet_precision_over_many_roundtrips},
    {"calculate_result_metrics", test_calculate_result_metrics},
    {"calculate_result_metrics_degenerate", test_calculate_result_metrics_degenerate},
    {"supertrend_dir_only_matches_full", test_supertrend_dir_only_matches_full},
    {"moving_averages", test_moving_averages},
    {"macd", test_macd},
    {"stoch_and_aroon", test_stoch_and_aroon},
    {"momentum_family", test_momentum_family},
    {"trend_strength_family", test_trend_strength_family},
    {"volatility_family", test_volatility_family},
    {"volume_family", test_volume_family},
    {"price_transforms_and_heikin_ashi", test_price_transforms_and_heikin_ashi},
    {"resample_timeframe", test_resample_timeframe},
    {"resample_timeframe_off_boundary_start", test_resample_timeframe_off_boundary_start},
    {"project_htf_to_ltf_has_no_lookahead", test_project_htf_to_ltf_has_no_lookahead},
    {"first_tradable_index", test_first_tradable_index},
    {"strategy_param_structs_are_fully_initialized", test_strategy_param_structs_are_fully_initialized},
    {"random_number_generator_seeding", test_random_number_generator_seeding},
    {"record_wallet_snapshot_drawdown", test_record_wallet_snapshot_drawdown},
    {"futures_long_close_sign", test_futures_long_close_sign},
    {"futures_short_close_sign", test_futures_short_close_sign},
    {"integer_range", test_integer_range},
    {"generate_range_int", test_generate_range_int},
    {"utc_timestamp_helpers", test_utc_timestamp_helpers},
    {"convert_to_unix_timestamp_is_utc", test_convert_to_unix_timestamp_is_utc},
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
