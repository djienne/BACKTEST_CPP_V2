// Independent small ledgers and causal properties through the production executor.
// No expected output is calculated by calling the implementation under test.
#include "engine/strategy_runner.hh"
#include "strategies/bbtrend.hh"
#include <stdexcept>
#include <sys/wait.h>
#include <csignal>
using trade_core::Intent;
static unsigned checks = 0;
void require(bool ok, const char *what)
{
    ++checks;
    if (!ok)
        throw std::runtime_error(what);
}
void near(double got, double expected, double tol = 1e-8)
{
    ++checks;
    if (!std::isfinite(got) || std::abs(got - expected) > tol)
        throw std::runtime_error("Expected " + std::to_string(expected) + ", got " + std::to_string(got));
}
MarketData market(std::vector<std::array<float, 4>> bars, int step = 60)
{
    MarketData d;
    KLINEf k{};
    k.name = "hand candles";
    for (size_t i = 0; i < bars.size(); ++i)
    {
        k.timestamp.push_back(1704067200 + i * step);
        k.open.push_back(bars[i][0]);
        k.high.push_back(bars[i][1]);
        k.low.push_back(bars[i][2]);
        k.close.push_back(bars[i][3]);
        k.volume.push_back(1);
    }
    k.nb = bars.size();
    d.signal = {k};
    d.execution = {k};
    d.funding.resize(1);
    d.signal_seconds = d.execution_seconds = step;
    d.start = k.timestamp[0];
    d.end = k.timestamp.back() + step;
    return d;
}
void fees()
{
    trade_core::PortfolioState<1> s(1001);
    trade_core::TradeStats stats;
    trade_core::open_futures_short(s, stats, 0, 100, 0.1, 1);
    near(s.coin_amounts[0], -10);
    near(s.usdt_amount, 0);
    trade_core::close_futures_short(s, stats, 0, 50, 0.1);
    near(s.usdt_amount, 1499.5);
    near(s.total_fees_paid_usdt, 1.5);
    trade_core::PortfolioState<1> spot(1000);
    stats = {};
    trade_core::open_spot_long(spot, stats, 0, 100, 0.1, 1);
    trade_core::close_spot_long(spot, stats, 0, 100, 0.1);
    near(spot.usdt_amount, 998.001);
    require(stats.nb_profit == 0 && stats.nb_loss == 1, "Fee-losing trade must lose");
}
void fills()
{
    auto d = market({{100, 101, 99, 100}, {110, 116, 109, 115}, {110, 125, 95, 115}, {115, 116, 114, 115}});
    auto signal = [](uint, size_t i)
    {
        Intent s;
        if (i == 0)
        {
            s.entry = 1;
            s.stop_distance = 10;
            s.target_distance = 10;
        }
        return s;
    };
    auto r = trade_core::simulate(d, {0, 4}, false, 1, signal, true, 1000, 0);
    require(r.fills.size() == 2, "One round trip");
    near(r.fills[0].price, 110);
    near(r.fills[1].price, 100);
    require(r.fills[0].timestamp == d.start + 60, "Next open, not signal close");
    require(r.fills[1].action == "stop" && r.ambiguous_bars == 1, "Stop first on ambiguous candle");
    // Exact cash ledger: buy 1000/110 units, sell at 100.
    near(r.WALLET_VAL_USDT, 909.090909090909);
    for (const auto &scenario : std::vector<std::pair<float, std::string>>{{80, "gap_stop"}, {130, "gap_target"}})
    {
        auto gap = d;
        gap.execution[0].open[2] = scenario.first;
        gap.execution[0].high[2] = scenario.first + 1;
        gap.execution[0].low[2] = scenario.first - 1;
        gap.execution[0].close[2] = scenario.first;
        auto out = trade_core::simulate(gap, {0, 4}, false, 1, signal, true, 1000, 0);
        near(out.fills[1].price, scenario.first);
        require(out.fills[1].action == scenario.second, "Gap opening execution");
    }
}
void trailing_and_timeouts()
{
    auto d = market({{100, 101, 99, 100}, {100, 115, 99, 110}, {105, 111, 95, 105}, {105, 106, 104, 105}});
    auto r = trade_core::simulate(
        d, {0, 4}, false, 1,
        [](uint, size_t i)
        {
        Intent s;
        if (i == 0)
        {
            s.entry = 1;
            s.stop_distance = 10;
            s.trailing = true;
        }
        return s;
        },
        true, 1000, 0);
    require(r.fills.size() == 2 && r.fills[1].timestamp == d.start + 180, "Trailing update cannot act retroactively");
    near(r.fills[1].price, 100);
    auto timeout = trade_core::simulate(
        d, {0, 4}, false, 1,
        [](uint, size_t i)
        {
        Intent s;
        if (i == 0)
        {
            s.entry = 1;
            s.max_hold_seconds = 120;
        }
        return s;
        },
        true, 1000, 0);
    require(timeout.fills[1].action == "timeout" && timeout.fills[1].timestamp == d.start + 180,
            "Timeout measured from fill");
    near(timeout.fills[1].price, 105);
}
void funding_and_drawdown()
{
    auto d = market({{100, 101, 99, 100}, {100, 101, 99, 100}, {100, 101, 99, 100}, {100, 101, 99, 100}});
    d.funding[0] = {{d.start + 60, 0.01, 110}, {d.start + 120, 0.01, 110}};
    for (int side : {1, -1})
    {
        auto r = trade_core::simulate(
            d, {0, 4}, true, 1,
            [=](uint, size_t i)
            {
            Intent s;
            if (i == 0)
                s.entry = side;
            return s;
            },
            true, 1000, 0);
        // At the entry timestamp the account was flat immediately before settlement.
        // Only the second funding event applies: ten units * 110 * 1% = 11.
        near(r.net_funding, side * 11);
        near(r.WALLET_VAL_USDT, 1000 - side * 11);
        near(r.total_fees_paid, 0);
        near(r.win_rate, side == 1 ? 0 : 100);
    }
    auto dip = market({{100, 101, 99, 100}, {100, 101, 99, 100}, {50, 51, 49, 50}, {100, 101, 99, 100}});
    auto r = trade_core::simulate(
        dip, {0, 4}, false, 1,
        [](uint, size_t i)
        {
        Intent s;
        if (i == 0)
            s.entry = 1;
        return s;
        },
        true);
    near(r.max_DD, -50.05);
    near(r.WALLET_VAL_USDT, 998.001);
    auto peak = market(
        {{100, 101, 99, 100}, {100, 101, 99, 100}, {200, 201, 199, 200}, {100, 101, 99, 100}, {180, 181, 179, 180}});
    auto recovered = trade_core::simulate(
        peak, {0, 5}, false, 1,
        [](uint, size_t i)
        {
        Intent s;
        if (i == 0)
            s.entry = 1;
        return s;
        },
        true, 1000, 0);
    near(recovered.max_DD, -50); // peak equity 2000 falls to 1000, then recovers
    near(recovered.WALLET_VAL_USDT, 1800);
}
void ordering_and_insolvency()
{
    auto d = market({{100, 101, 99, 100}, {100, 101, 99, 100}, {120, 121, 119, 120}, {120, 121, 119, 120}});
    d.execution.push_back(d.execution[0]);
    d.signal.push_back(d.signal[0]);
    d.funding.resize(2);
    auto r = trade_core::simulate(
        d, {0, 4}, false, 1,
        [](uint p, size_t i)
        {
        Intent s;
        if (p == 1 && i == 0)
            s.entry = 1;
        if (p == 1 && i == 1)
            s.exit_long = true;
        if (p == 0 && i == 1)
            s.entry = 1;
        return s;
        },
        true, 1000, 0);
    require(r.fills.size() == 4, "All pair exits free cash before the first pair's entry");
    require(r.fills[1].pair == 1 && r.fills[2].pair == 0, "Exit precedes allocation in coin order");
    near(r.WALLET_VAL_USDT, 1200);
    auto flat = market(std::vector<std::array<float, 4>>(4, {100, 101, 99, 100}));
    flat.funding[0] = {{flat.start + 60, 0.01, 110, 1}};
    r = trade_core::simulate(
        flat, {0, 4}, true, 1,
        [](uint, size_t i)
        {
        Intent s;
        if (i == 0)
            s.entry = 1;
        return s;
        },
        true, 1000, 0);
    near(r.net_funding, 11);
    near(r.WALLET_VAL_USDT, 989);
    auto boom = market({{100, 101, 99, 100}, {100, 101, 99, 100}, {210, 211, 209, 210}, {210, 211, 209, 210}});
    r = trade_core::simulate(
        boom, {0, 4}, true, 1,
        [](uint, size_t i)
        {
        Intent s;
        if (i == 0)
            s.entry = -1;
        return s;
        },
        true, 1000, 0);
    require(!r.valid, "Insolvent short candidates are invalid");
    r = trade_core::simulate(
        flat, {0, 4}, false, 1,
        [](uint, size_t i)
        {
        Intent s;
        if (i == 2)
            s.entry = 1;
        return s;
        },
        true, 1000, 0);
    require(r.fills.empty(), "No opening position on the terminal candle");
    for (float op : {80, 120})
    {
        auto gap = market({{100, 101, 99, 100}, {100, 101, 99, 100}, {op, op + 1, op - 1, op}, {100, 101, 99, 100}});
        r = trade_core::simulate(
            gap, {0, 4}, true, 1,
            [](uint, size_t i)
            {
            Intent s;
            if (i <= 1)
            {
                s.entry = -1;
                s.stop_distance = 10;
                s.target_distance = 10;
            }
            return s;
            },
            true, 1000, 0);
        require(r.fills.size() == 2, "Protective gap exit cannot reenter the same candle");
        require(r.fills[1].action == (op < 100 ? "gap_target" : "gap_stop"), "Short gap direction");
        near(r.fills[1].price, op);
        near(r.WALLET_VAL_USDT, op < 100 ? 1200 : 800);
    }
}
void shared_strategy_models()
{
    // With EMA(2) and 3-candle Bollinger bands, close 110 breaks the upper
    // band at index 3; close 90 breaks the lower at index 5. Orders follow at 4/6.
    auto d = market({{100, 101, 99, 100},
                     {100, 101, 98, 99},
                     {99, 101, 98, 100},
                     {100, 111, 99, 110},
                     {112, 113, 99, 100},
                     {100, 101, 89, 90},
                     {88, 89, 84, 85},
                     {85, 86, 84, 85}});
    strategy_runner::init_talib();
    for (bool futures : {false, true})
    {
        std::vector<IndicatorCache> cache(1);
        auto r = strategies::BBTrend{futures}(d, cache, {2, 3, 1, 1}, {0, 8}, true);
        require(r.fills.size() == (futures ? 4 : 2), "Shared model keeps spot long-only and enables futures shorts");
        require(r.fills[0].timestamp == d.start + 4 * 60, "Upper-band signal uses the next opening");
        near(r.fills[0].price, 112);
        if (futures)
        {
            require(r.fills[2].quantity < 0 && r.fills[2].timestamp == d.start + 6 * 60,
                    "Lower-band break opens a short");
            near(r.fills[2].price, 88);
        }
    }
    TA_Shutdown();
}
void causality_and_holdout()
{
    std::vector<std::array<float, 4>> candles;
    // Training rises from entry price 100 to 120. A fixed long beats a fixed
    // short. The holdout later crashes, so using it for selection reverses that choice.
    for (int i = 0; i < 30; ++i)
    {
        const float price = i < 2 ? 100 : i < 24 ? 100 + 20.0f * (i - 1) / 22 : 120 + 5.0f * (i - 23);
        candles.push_back({price, price + 1, price - 1, price});
    }
    auto d = market(candles);
    auto evaluate =
        [](const MarketData &data, auto &, const strategy_runner::Params &p, trade_core::Window w, bool trace)
    {
        return trade_core::simulate(
            data, w, true, 1,
            [=](uint, size_t)
            {
            Intent s;
            s.entry = static_cast<int>(p[0]);
            return s;
            },
            trace);
    };
    backtest_config::StrategyConfig cfg;
    cfg.workers = 1;
    cfg.holdout_fraction = 0.2;
    const std::vector<strategy_runner::Params> params{{1}, {-1}};
    strategy_runner::Selection selection{0, -100, -1e6, 1e6};
    auto first = strategy_runner::evaluate_search(d, cfg, params, selection, evaluate);
    require(first.at("parameters") == nlohmann::json::array({1.0}), "Rising training prices select the long");
    near(first.at("training").at("wallet"), 1198800.0 / 1001);
    cfg.workers = 2;
    auto threaded = strategy_runner::evaluate_search(d, cfg, params, selection, evaluate);
    require(first == threaded, "Worker count must not change results or winner");
    auto future = d;
    for (size_t i = 24; i < 30; ++i)
    {
        const float price = 120 - 20.0f * (i - 24);
        future.execution[0].open[i] = future.execution[0].close[i] = price;
        future.execution[0].high[i] = price + 1;
        future.execution[0].low[i] = price - 1;
    }
    future.signal[0] = future.execution[0];
    auto changed = strategy_runner::evaluate_search(future, cfg, params, selection, evaluate);
    require(first.at("parameters") == changed.at("parameters"), "Holdout must not select parameters");
    require(first.at("training") == changed.at("training"), "Future changes must not change training trades");
    require(first.at("holdout").at("wallet") != changed.at("holdout").at("wallet"),
            "Holdout perturbation must change actual returns");
    near(first.at("holdout").at("equity")[0], 1000);
    auto parameters =
        strategy_runner::sample(std::vector<strategy_runner::Axis>{{"period", integer_range(2, 200)}}, 100, 42,
                                [](const strategy_runner::Params &)
                                {
        return true;
        });
    std::set<strategy_runner::Params> distinct(parameters.begin(), parameters.end());
    require(distinct.size() == 100, "Distinct budgeted candidates");
}
void input_and_readiness()
{
    auto d = market(std::vector<std::array<float, 4>>(1000, {100, 101, 99, 100}), 3600);
    auto pairs = std::vector<KLINEf>{d.signal[0], d.signal[0]};
    auto starts = INITIALIZE_DATA(pairs);
    require(pairs[1].nb == 1000 && starts[1] == 0, "No arbitrary deletion of hourly history");
    pairs[1].timestamp[20] += 1;
    bool rejected = false;
    try
    {
        INITIALIZE_DATA(pairs);
    }
    catch (const std::exception &)
    {
        rejected = true;
    }
    require(rejected, "Reject timestamp mismatch before indexing");
    size_t ready = 999;
    TA_Initialize();
    std::vector<IndicatorCache> cache(1);
    strategy_runner::Indicators ind{d, cache};
    ind.ema(800);
    auto guarded = strategy_runner::evaluate(
        d, {0, 1000}, false, 1, ind.ready,
        [](uint, size_t)
        {
        Intent s;
        s.entry = 1;
        return s;
        },
        true);
    require(guarded.fills.size() == 2, "Readiness production path opens one position");
    require(guarded.fills[0].timestamp == d.start + 800 * 3600, "EMA warmup gates the first actual fill");
    auto e = TALIB_EMA(std::vector<float>(5, 100), 20, &ready);
    require(ready == e.size(), "Insufficient EMA history has no valid values");
    auto prefix = TALIB_EMA(std::vector<float>(50, 100), 20);
    auto extended = std::vector<float>(50, 100);
    extended.insert(extended.end(), 50, 10000);
    auto full = TALIB_EMA(extended, 20);
    for (size_t i = 0; i < prefix.size(); ++i)
        near(prefix[i], full[i]);
    auto histogram = TALIB_TRIX(std::vector<float>(80, 100), 5, 4);
    for (float value : histogram)
        near(value, 0);
    // Valid-stage EMA seeding gives a zero histogram for a constant price from
    // its first usable value; padding must not create a false initial trend.
    auto trend = std::vector<float>(80, 100);
    trend.insert(trend.end(), 20, 200);
    auto future_trix = TALIB_TRIX(trend, 5, 4);
    for (size_t i = 0; i < histogram.size(); ++i)
        near(future_trix[i], histogram[i]);
    near(TALIB_TRIX({100, 100, 100, 100, 100, 100, 200}, 2, 2).back(), 400.0 / 27, 1e-4);
    TA_Shutdown();
    auto gaps = market(std::vector<std::array<float, 4>>(24, {100, 101, 99, 100}), 300);
    for (size_t i = 1; i < 24; ++i)
        gaps.signal[0].timestamp[i] += 300;
    std::cout.flush();
    pid_t pid = fork();
    if (pid == 0)
    {
        freopen("/dev/null", "w", stderr);
        RESAMPLE_TIMEFRAME(gaps.signal[0], 12, 5, 60);
        _Exit(0);
    }
    int status = 0;
    waitpid(pid, &status, 0);
    require(WIFSIGNALED(status) && WTERMSIG(status) == SIGABRT, "Resampler must reject missing bars");
}
int main()
{
    try
    {
        fees();
        fills();
        trailing_and_timeouts();
        funding_and_drawdown();
        ordering_and_insolvency();
        shared_strategy_models();
        causality_and_holdout();
        input_and_readiness();
        std::cout << "Execution checks passed: " << checks << "\n";
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "FAIL: " << e.what() << "\n";
        return 1;
    }
}
