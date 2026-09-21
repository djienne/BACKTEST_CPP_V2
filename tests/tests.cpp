// Analytic component checks and account arithmetic. Run in Docker Compose with
// bash run_regression.sh; production execution is checked in execution_tests.cpp.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "engine/indicators.hh"
#include "engine/strategy_runner.hh"
#include "engine/tools.hh"
#include "engine/trade_core.hh"
#include <ta-lib/ta_libc.h>

namespace
{
int g_fail_count = 0;
int g_run_count = 0;

#define REQUIRE(expr)                                                                                                  \
    do                                                                                                                 \
    {                                                                                                                  \
        ++g_run_count;                                                                                                 \
        if (!(expr))                                                                                                   \
        {                                                                                                              \
            ++g_fail_count;                                                                                            \
            std::cerr << "FAIL " << __FUNCTION__ << ": " #expr " at " << __FILE__ << ":" << __LINE__ << "\n";          \
        }                                                                                                              \
    } while (0)

#define REQUIRE_NEAR(a, b, eps)                                                                                        \
    do                                                                                                                 \
    {                                                                                                                  \
        ++g_run_count;                                                                                                 \
        const double _a = static_cast<double>(a);                                                                      \
        const double _b = static_cast<double>(b);                                                                      \
        if (!std::isfinite(_a) || !std::isfinite(_b) || std::fabs(_a - _b) > (eps))                                    \
        {                                                                                                              \
            ++g_fail_count;                                                                                            \
            std::cerr << "FAIL " << __FUNCTION__ << ": " #a "=" << _a << " vs " #b "=" << _b << " (tol=" << (eps)      \
                      << ") at " << __FILE__ << ":" << __LINE__ << "\n";                                               \
        }                                                                                                              \
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
    stats.nb_closed = 10;
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
    winners.nb_closed = 10;
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

void test_supertrend_known_reversal()
{
    const std::vector<float> close{100, 100, 100, 103, 97};
    const std::vector<float> high{101, 101, 101, 104, 98}, low{99, 99, 99, 102, 96};
    const auto st = TALIB_SuperTrend(high, low, close, 2, 1.0f);
    // ATR at index 2 is 2; at index 3 it is (2+4)/2=3. The lower
    // band rises from 98 to 100; the next close at 97 breaks it downward.
    REQUIRE(st.warmup == 2);
    REQUIRE_NEAR(st.final_lowerband[2], 98, 1e-6);
    REQUIRE_NEAR(st.final_lowerband[3], 100, 1e-6);
    REQUIRE(st.supertrend[2] == 1 && st.supertrend[3] == 1 && st.supertrend[4] == -1);
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
    std::vector<float> series;
    for (int i = 0; i < 200; ++i)
        series.push_back(100.0f + i);
    size_t warm = 0;
    const auto ema = TALIB_EMA(series, 20, &warm);
    REQUIRE(ema.size() == 200 && warm == 19);
    REQUIRE_NEAR(ema[18], 0, 1e-6);
    // SMA seed is 109.5; a unit ramp then keeps its exact 9.5-bar lag.
    REQUIRE_NEAR(ema[19], 109.5, 1e-6);
    REQUIRE_NEAR(ema[199], 289.5, 1e-6);
}

void test_talib_bbands_known_values()
{
    const auto b = TALIB_BBANDS_R({1, 2, 3, 4, 5}, 2, 2, 3);
    // Last window is [3,4,5]: mean 4, population variance 2/3.
    REQUIRE(b.warmup == 2 && b.upper.size() == 5);
    REQUIRE_NEAR(b.middle[4], 4, 1e-6);
    REQUIRE_NEAR(b.upper[4], 5.63299316185545, 1e-6);
    REQUIRE_NEAR(b.lower[4], 2.36700683814455, 1e-6);
}

void test_talib_ao_ramp()
{
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
    // Difference between 5- and 34-point averages of a unit ramp: (34-5)/2.
    REQUIRE_NEAR(ao.back(), 14.5, 1e-6);
}

void test_talib_stochrsi_known_values()
{
    const auto raw = TALIB_STOCHRSI_not_averaged({100, 101, 100, 102, 100, 103, 100}, 3, 2);
    // Wilder RSI at indices 3,4,5 is 83.333333,35.714286,76.315789.
    // Normalizing the last value in that window gives 0.8526316, rounded to .853.
    REQUIRE_NEAR(raw[4], 0, 1e-6);
    REQUIRE_NEAR(raw[5], 0.853, 1e-6);
    REQUIRE_NEAR(raw[6], 0, 1e-6);
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
    REQUIRE_NEAR(state.net_funding, expected_deduction, 1e-9);

    // Zero funding -> no-op even with a position.
    const double after = state.usdt_amount;
    trade_core::apply_funding_fee(state, 0, 120.0, 0.0);
    REQUIRE_NEAR(state.usdt_amount, after, 1e-12);
}

void test_calculate_calmar_ratio()
{
    // Doubling over two 365.25-day years gives sqrt(2)-1 annual growth;
    // divided by a 25% drawdown, the ratio is 1.65685424949238.
    auto r = trade_core::calmar_ratio(1000, 2000, 63115200, -25);
    REQUIRE(r.has_value());
    REQUIRE_NEAR(*r, 1.65685424949238, 1e-12);
    REQUIRE(!trade_core::calmar_ratio(1000, 2000, 63115200, 0));
    REQUIRE(*trade_core::calmar_ratio(1000, 1140, 14 * 86400, -10) > 0);
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
        base.timestamp.push_back(1020 + i * 60);
        base.open.push_back(100.0f + i);
        base.high.push_back(101.0f + i);
        base.low.push_back(99.0f + i);
        base.close.push_back(100.0f + i);
        base.volume.push_back(1);

        other.timestamp.push_back(1020 + i * 60);
        other.open.push_back(200.0f + i);
        other.high.push_back(201.0f + i);
        other.low.push_back(199.0f + i);
        other.close.push_back(200.0f + i);
        other.volume.push_back(1);
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

// Mirrors the long take-profit / stop-loss test SuperTrend_EMA_ATR and F_3EMA_SRSI_ATR
// both use, so the convention is pinned in one place.
// ---------------------------------------------------------------------------
//  Indicator library
//
//  Hand-computed values and limiting cases check numerical meaning. Shape and
//  readiness checks supplement them; neither alone establishes correctness.
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

    // DEMA and TEMA cancel a linear trend's lag after their SMA-seeded warmups.
    const std::vector<float> dema = TALIB_DEMA(ramp, 10);
    const std::vector<float> tema = TALIB_TEMA(ramp, 10);
    REQUIRE_NEAR(dema[90], 90, 1e-5);
    REQUIRE_NEAR(tema[90], 90, 1e-5);

    // HMA is composed here rather than taken from TA-Lib; it must reduce lag too and
    // report a warmup that actually covers its zero-padded head.
    size_t warm_hma = 0;
    const std::vector<float> hma = TALIB_HMA(ramp, 16, &warm_hma);
    REQUIRE(hma.size() == ramp.size());
    REQUIRE(warm_hma > 0);
    REQUIRE(hma[warm_hma - 1] == 0.0f);
    REQUIRE_NEAR(hma[90], 89.3333333333333, 1e-4); // HMA(16) ramp lag is 2/3

    size_t warm_kama = 0;
    const std::vector<float> kama = TALIB_KAMA(ramp, 10, &warm_kama);
    REQUIRE(kama.size() == ramp.size());
    REQUIRE_NEAR(kama[90], 88.75, 1e-5); // efficiency 1, smoothing 4/9, steady lag 5/4
}

void test_macd()
{
    const auto m = TALIB_MACD({10, 10, 10, 10, 10, 10, 16, 10, 10}, 3, 5, 2);
    REQUIRE(m.warmup == 5 && m.macd.size() == 9);
    // After the jump to 16, fast/slow EMAs are 13 and 12. The signal's
    // smoothing weight is 2/3. One candle later the MACD is 1/6.
    REQUIRE_NEAR(m.macd[6], 1, 1e-6);
    REQUIRE_NEAR(m.signal[6], 2.0 / 3, 1e-6);
    REQUIRE_NEAR(m.histogram[6], 1.0 / 3, 1e-6);
    REQUIRE_NEAR(m.macd[7], 1.0 / 6, 1e-6);
    REQUIRE_NEAR(m.signal[7], 1.0 / 3, 1e-6);
    REQUIRE_NEAR(m.histogram[7], -1.0 / 6, 1e-6);
}

void test_stoch_and_aroon()
{
    std::vector<float> high, low, close;
    for (int i = 0; i < 60; ++i)
    {
        high.push_back(102 + i);
        low.push_back(98 + i);
        close.push_back(100 + i);
    }
    const auto st = TALIB_STOCH(high, low, close, 14, 3, 3);
    // Every mature 14-bar window has close-lowest=15 and highest-lowest=17.
    REQUIRE_NEAR(st.k[50], 1500.0 / 17, 1e-5);
    REQUIRE_NEAR(st.d[50], 1500.0 / 17, 1e-5);
    const auto ar = TALIB_AROON(high, low, 14);
    REQUIRE_NEAR(ar.up[50], 100, 1e-6);
    REQUIRE_NEAR(ar.down[50], 0, 1e-6);
}

void test_momentum_family()
{
    const std::vector<float> ramp = make_ramp(100);

    // MOM(n) on a unit ramp is exactly n.
    size_t warm_mom = 0;
    const std::vector<float> mom = TALIB_MOM(ramp, 10, &warm_mom);
    REQUIRE(warm_mom == 10);
    REQUIRE_NEAR(mom[50], 10.0f, 1e-3);

    // ROC(n) at i = (v[i] - v[i-n]) / v[i-n] * 100 = 10/40*100 at i=50.
    const std::vector<float> roc = TALIB_ROC(ramp, 10);
    REQUIRE_NEAR(roc[50], 25.0f, 1e-3);

    size_t warm_cci = 0;
    const auto cci = TALIB_CCI(ramp, ramp, ramp, 20, &warm_cci);
    REQUIRE(warm_cci == 19);
    REQUIRE_NEAR(cci[50], 126.6666666667, 1e-4); // 9.5 / (.015 * mean deviation 5)
    const auto ult = TALIB_ULTOSC(ramp, ramp, ramp, 7, 14, 28);
    REQUIRE_NEAR(ult[50], 100, 1e-6); // buying pressure equals true range in each window
}

void test_trend_strength_family()
{
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
    // Each upward move is 2 and true range is 3; DX/ADX therefore converge to 100.
    REQUIRE_NEAR(up.plus_di[100], 200.0 / 3, 1e-4);
    REQUIRE_NEAR(up.minus_di[100], 0, 1e-6);
    REQUIRE_NEAR(up.adx[100], 100, 1e-4);

    // Parabolic SAR trails an uptrend from below.
    size_t warm_sar = 0;
    const std::vector<float> sar = TALIB_SAR(up_h, up_l, 0.02, 0.2, &warm_sar);
    REQUIRE(sar.size() == up_c.size());
    REQUIRE(sar[100] < up_c[100]);
}

void test_volatility_family()
{
    const std::vector<float> high(60, 52), low(60, 48), close(60, 50);
    const auto natr = TALIB_NATR(high, low, close, 14);
    REQUIRE_NEAR(natr[40], 8, 1e-6); // 4 / 50 * 100
    const auto tr = TALIB_TRANGE({52, 52, 72}, {48, 48, 68}, {50, 50, 70});
    REQUIRE_NEAR(tr[1], 4, 1e-6);
    REQUIRE_NEAR(tr[2], 22, 1e-6); // opening gap dominates high-low
    const auto flat_sd = TALIB_STDDEV(close, 20, 1);
    REQUIRE_NEAR(flat_sd[40], 0, 1e-6);
    const auto ramp_sd = TALIB_STDDEV({1, 2, 3, 4, 5}, 3, 1);
    REQUIRE_NEAR(ramp_sd[4], 0.816496580927726, 1e-6);
    const auto kc = KELTNER_CHANNELS(high, low, close, 20, 10, 2);
    REQUIRE_NEAR(kc.middle[40], 50, 1e-6);
    REQUIRE_NEAR(kc.upper[40], 58, 1e-6);
    REQUIRE_NEAR(kc.lower[40], 42, 1e-6);
    const auto dc = DONCHIAN_CHANNELS({11, 15, 12, 13, 14}, {9, 7, 10, 11, 12}, 3);
    REQUIRE_NEAR(dc.upper[3], 15, 1e-6);
    REQUIRE_NEAR(dc.lower[3], 7, 1e-6);
    REQUIRE_NEAR(dc.middle[3], 11, 1e-6);
    REQUIRE_NEAR(dc.upper[4], 14, 1e-6); // old extrema have left the window
    REQUIRE_NEAR(dc.lower[4], 10, 1e-6);
    REQUIRE_NEAR(dc.middle[4], 12, 1e-6);
}

void test_volume_family()
{

    // OBV has an exact hand-checkable definition: add volume on an up close, subtract
    // on a down close, carry on unchanged.
    const std::vector<float> close{10.0f, 11.0f, 10.5f, 10.5f, 12.0f};
    const std::vector<float> vol{100.0f, 200.0f, 300.0f, 400.0f, 500.0f};
    const std::vector<float> obv = TALIB_OBV(close, vol);
    REQUIRE(obv.size() == close.size());
    REQUIRE_NEAR(obv[0], 100.0f, 1e-3); // seed
    REQUIRE_NEAR(obv[1], 300.0f, 1e-3); // up   -> +200
    REQUIRE_NEAR(obv[2], 0.0f, 1e-3);   // down -> -300
    REQUIRE_NEAR(obv[3], 0.0f, 1e-3);   // flat -> unchanged
    REQUIRE_NEAR(obv[4], 500.0f, 1e-3); // up   -> +500

    const std::vector<float> rising{10, 11, 12, 13, 14}, volume(5, 100);
    REQUIRE_NEAR(TALIB_MFI(rising, rising, rising, volume, 3)[4], 100, 1e-6);
    const std::vector<float> falling{14, 13, 12, 11, 10};
    REQUIRE_NEAR(TALIB_MFI(falling, falling, falling, volume, 3)[4], 0, 1e-6);
    const std::vector<float> high(5, 12), low(5, 8), last{11, 9, 10, 12, 8};
    const auto ad = TALIB_AD(high, low, last, vol);
    // Money-flow multipliers .5,-.5,0,1,-1 give cumulative flows 50,-50,-50,350,-150.
    REQUIRE_NEAR(ad[0], 50, 1e-6);
    REQUIRE_NEAR(ad[3], 350, 1e-6);
    REQUIRE_NEAR(ad[4], -150, 1e-6);
    const auto osc = TALIB_ADOSC(high, low, last, vol, 2, 3);
    // EMAs seeded at flow 50 have values -350/9 and -25 after the third candle.
    REQUIRE_NEAR(osc[2], -125.0 / 9, 1e-5);
    size_t warm_vwap = 0;
    const auto vwap = VWAP_ROLLING({10, 20, 30}, {10, 20, 30}, {10, 20, 30}, {1, 2, 3}, 2, &warm_vwap);
    REQUIRE(warm_vwap == 1);
    REQUIRE_NEAR(vwap[1], 50.0 / 3, 1e-5);
    REQUIRE_NEAR(vwap[2], 26, 1e-6);
    const auto rv = RELATIVE_VOLUME(std::vector<float>(60, 10), 20);
    REQUIRE_NEAR(rv[40], 1, 1e-6);
    std::vector<float> spiky(60, 10);
    spiky[50] = 100;
    const auto spike = RELATIVE_VOLUME(spiky, 20);
    REQUIRE_NEAR(spike[50], 200.0 / 29, 1e-6); // 100 against a window mean of 14.5
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
    REQUIRE_NEAR(ha.open[1], 11, 1e-5);
    REQUIRE_NEAR(ha.high[1], 16, 1e-6);
    REQUIRE_NEAR(ha.low[1], 9, 1e-6);
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
    // A partial first hour must be skipped, with its offset preserved in projection.
    KLINEf in{};
    const int64_t hour = 1687104000;  // an exact hour
    const int64_t start = hour - 300; // 5 minutes earlier: 15:55
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
    // Bar 0 is 15:55; bar 1 is 16:00 and starts the first whole hour.
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
    // Bars 0..11 close before the first higher candle completes.
    REQUIRE(proj[0] == -1.0f);
    REQUIRE(proj[11] == -1.0f);
    // The first higher value is visible at bar 12's close, the second at bar 24's.
    REQUIRE(proj[12] == 7.0f);
    REQUIRE(proj[23] == 7.0f);
    REQUIRE(proj[24] == 8.0f);
}

void test_project_htf_to_ltf_has_no_lookahead()
{
    // Three higher-timeframe values projected onto 36 lower-timeframe bars.
    const std::vector<float> htf{10.0f, 20.0f, 30.0f};
    const std::vector<float> ltf = PROJECT_HTF_TO_LTF(htf, 12, 36, 0, -777.0f);

    REQUIRE(ltf.size() == 36);

    // Bars 0..10 close before the higher candle is complete; bar 11 completes it.
    bool head_is_fill = true;
    for (size_t i = 0; i < 11; ++i)
    {
        head_is_fill = head_is_fill && ltf[i] == -777.0f;
    }
    REQUIRE(head_is_fill);

    // Each published close remains available until the next higher close.
    bool block1 = true, block2 = true;
    for (size_t i = 11; i < 23; ++i)
    {
        block1 = block1 && ltf[i] == 10.0f;
    }
    for (size_t i = 23; i < 35; ++i)
    {
        block2 = block2 && ltf[i] == 20.0f;
    }
    REQUIRE(block1);
    REQUIRE(block2);
    REQUIRE(ltf[35] == 30.0f); // available at this bar close, for next-open execution
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
    {"supertrend_known_reversal", test_supertrend_known_reversal},
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
    {"talib_bbands_known_values", test_talib_bbands_known_values},
    {"talib_ao_ramp", test_talib_ao_ramp},
    {"talib_stochrsi_known_values", test_talib_stochrsi_known_values},
    {"apply_funding_fee", test_apply_funding_fee},
    {"calculate_calmar_ratio", test_calculate_calmar_ratio},
    {"realign_timestamps_noop_when_aligned", test_realign_timestamps_noop_when_aligned},
};
} // namespace

int main()
{
    if (TA_Initialize() != TA_SUCCESS)
    {
        std::cerr << "FAIL: TA-Lib initialization\n";
        return 1;
    }
    for (const NamedTest &t : ALL_TESTS)
    {
        const int failures_before = g_fail_count;
        t.fn();
        const bool ok = g_fail_count == failures_before;
        std::cout << (ok ? "PASS " : "FAIL ") << t.name << "\n";
    }

    TA_Shutdown();
    std::cout << "\n" << (g_run_count - g_fail_count) << " / " << g_run_count << " checks passed\n";
    if (g_fail_count != 0)
    {
        std::cout << "TESTS FAILED: " << g_fail_count << "\n";
        return 1;
    }
    std::cout << "ALL TESTS PASSED\n";
    return 0;
}
