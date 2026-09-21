#pragma once
#include "trade_core.hh"
#include "indicators.hh"
#include <functional>
#include <future>
#include <numeric>
#include <unordered_set>
#ifndef BACKTEST_REVISION
#define BACKTEST_REVISION "unknown"
#endif

namespace strategy_runner
{
using Params = std::vector<double>;
using Column = std::array<const std::vector<float> *, backtest_config::MAX_PAIRS>;
struct Axis
{
    std::string name;
    std::vector<double> values;
    template <class T> Axis(std::string n, const std::vector<T> &v) : name(std::move(n)), values(v.begin(), v.end()) {}
};
struct Selection
{
    uint min_trades = 100;
    double min_dd = -40, min_gain = 0, max_gain = 1e6;
};
inline void init_talib()
{
    if (TA_Initialize() != TA_SUCCESS)
        throw std::runtime_error("TA-Lib initialization failed");
}
// The cache owns series; Columns keep stable references outside the bar loop.
struct Indicators
{
    const MarketData &data;
    std::vector<IndicatorCache> &cache;
    size_t ready = 0;
    template <class Compute> Column get(const std::string &key, Compute compute, bool higher = false)
    {
        Column out{};
        for (size_t p = 0; p < data.signal.size(); ++p)
        {
            auto &c = cache[p];
            if (!c.has(key))
            {
                size_t warm = 0;
                auto values = compute(higher ? data.higher.at(p) : data.signal[p], warm);
                if (higher)
                {
                    warm = data.higher_offset[p] + (warm + 1) * data.higher_ratio - 1;
                    values = PROJECT_HTF_TO_LTF(values, data.higher_ratio, data.signal[p].nb, data.higher_offset[p], 0);
                }
                c.put(key, std::move(values), warm);
            }
            out[p] = &c.get(key);
            ready = std::max(ready, c.first_valid(key));
        }
        return out;
    }
    Column ema(int n, bool higher = false)
    {
        return get(
            IndicatorCache::key(higher ? "EMA_HTF" : "EMA", n),
            [=](const KLINEf &k, size_t &w)
            {
            return TALIB_EMA(k.close, n, &w);
            },
            higher);
    }
    Column rsi(int n)
    {
        return get(IndicatorCache::key("RSI", n),
                   [=](const KLINEf &k, size_t &w)
                   {
            return TALIB_RSI(k.close, n, &w);
        });
    }
    Column atr(int n = 14)
    {
        return get(IndicatorCache::key("ATR", n),
                   [=](const KLINEf &k, size_t &w)
                   {
            return TALIB_ATR(k.high, k.low, k.close, n, &w);
        });
    }
    Column will()
    {
        return get("WILLR:14",
                   [](const KLINEf &k, size_t &w)
                   {
            return TALIB_WILLR(k.high, k.low, k.close, 14, &w);
        });
    }
    Column srsi(int smooth = 0)
    {
        return get(IndicatorCache::key("SRSI", smooth),
                   [=](const KLINEf &k, size_t &w)
                   {
            w = 27 + (smooth == 0 ? 0 : smooth == 1 ? 2 : 4);
            return smooth == 0   ? TALIB_STOCHRSI_not_averaged(k.close, 14, 14)
                   : smooth == 1 ? TALIB_STOCHRSI_K(k.close, 14, 14, 3, 3)
                                 : TALIB_STOCHRSI_D(k.close, 14, 14, 3, 3);
        });
    }
    Column ao(int fast, int slow)
    {
        return get(IndicatorCache::key("AO", fast, slow),
                   [=](const KLINEf &k, size_t &w)
                   {
            w = static_cast<size_t>(std::max(fast, slow) - 1);
            return TALIB_AO(k.high, k.low, fast, slow);
        });
    }
    Column trix(int length, int smooth)
    {
        return get(IndicatorCache::key("TRIX", length, smooth),
                   [=](const KLINEf &k, size_t &w)
                   {
            w = 3 * (length - 1) + smooth;
            return TALIB_TRIX(k.close, length, smooth);
        });
    }
    Column supertrend(int period, float mult, bool higher = false)
    {
        return get(
            IndicatorCache::key(higher ? "ST_HTF" : "ST", period, mult),
            [=](const KLINEf &k, size_t &w)
            {
            auto s = TALIB_SuperTrend(k.high, k.low, k.close, period, mult);
            w = s.warmup;
            return std::vector<float>(s.supertrend.begin(), s.supertrend.end());
            },
            higher);
    }
    std::array<Column, 3> bands(int length, float dev)
    {
        std::array<Column, 3> out{};
        const std::string base = IndicatorCache::key("BB", length, dev);
        for (size_t p = 0; p < data.signal.size(); ++p)
        {
            auto &c = cache[p];
            if (!c.has(base + ":0"))
            {
                auto b = TALIB_BBANDS_R(data.signal[p].close, dev, dev, length);
                c.put(base + ":0", std::move(b.upper), b.warmup);
                c.put(base + ":1", std::move(b.middle), b.warmup);
                c.put(base + ":2", std::move(b.lower), b.warmup);
            }
            for (int k = 0; k < 3; ++k)
            {
                const auto key = base + ":" + std::to_string(k);
                out[k][p] = &c.get(key);
                ready = std::max(ready, c.first_valid(key));
            }
        }
        return out;
    }
};
template <class Valid>
std::vector<Params> sample(const std::vector<Axis> &axes, uint64_t budget, uint64_t seed, Valid valid)
{
    uint64_t count = 1;
    for (const auto &a : axes)
    {
        if (a.values.empty() || count > std::numeric_limits<uint64_t>::max() / a.values.size())
            throw std::runtime_error("Invalid parameter grid");
        count *= a.values.size();
    }
    // Explicit exhaustive mode is bounded by memory. Above this ceiling use a trial
    // budget; a streaming grid would be the next step if a real study needs it.
    if (!budget && count > 10000000)
        throw std::runtime_error("Exhaustive grid exceeds 10 million combinations; set max_trials");
    std::vector<Params> result;
    std::mt19937_64 rng(seed);
    std::unordered_set<uint64_t> seen;
    const uint64_t wanted = budget ? std::min(budget, count) : count;
    for (uint64_t visits = 0; seen.size() < count && result.size() < wanted;)
    {
        uint64_t index = budget ? std::uniform_int_distribution<uint64_t>(0, count - 1)(rng) : visits++;
        if (!seen.insert(index).second)
            continue;
        Params p;
        for (const auto &a : axes)
        {
            p.push_back(a.values[index % a.values.size()]);
            index /= a.values.size();
        }
        if (valid(p))
            result.push_back(std::move(p));
    }
    if (result.empty())
        throw std::runtime_error("No valid parameter combinations");
    return result;
}
inline nlohmann::json result_json(const RUN_RESULTf &r)
{
    nlohmann::json j = {{"valid", r.valid},
                        {"invalid_reason", r.invalid_reason},
                        {"wallet", r.WALLET_VAL_USDT},
                        {"gain_percent", r.gain_pc},
                        {"win_rate_percent", r.win_rate},
                        {"max_drawdown_percent", r.max_DD},
                        {"score", r.score},
                        {"trades", r.nb_posi_entered},
                        {"commissions", r.total_fees_paid},
                        {"net_funding_paid", r.net_funding},
                        {"ambiguous_bars", r.ambiguous_bars}};
    j["calmar"] = r.calmar_ratio ? nlohmann::json(*r.calmar_ratio) : nlohmann::json(nullptr);
    j["fills"] = nlohmann::json::array();
    for (const auto &f : r.fills)
        j["fills"].push_back({{"timestamp", f.timestamp},
                              {"pair", f.pair},
                              {"action", f.action},
                              {"price", f.price},
                              {"quantity", f.quantity},
                              {"commission", f.commission},
                              {"net_pnl", f.net_pnl}});
    j["equity"] = r.equity;
    j["equity_times"] = r.equity_times;
    return j;
}
inline bool qualifies(const RUN_RESULTf &r, const Selection &s)
{
    return r.valid && std::isfinite(r.score) && r.gain_pc > s.min_gain && r.gain_pc < s.max_gain &&
           r.nb_posi_entered >= static_cast<int>(s.min_trades) && r.max_DD > s.min_dd;
}
template <class Evaluate>
nlohmann::json evaluate_search(const MarketData &data, const backtest_config::StrategyConfig &cfg,
                               const std::vector<Params> &params, const Selection &selection, Evaluate evaluate)
{
    const size_t n = data.execution[0].nb;
    const size_t split = static_cast<size_t>(std::floor(n * (1 - cfg.holdout_fraction)));
    if (!split || split >= n)
        throw std::runtime_error("Insufficient evaluation bars for holdout");
    std::vector<RUN_RESULTf> results(params.size());
    std::vector<std::future<void>> workers;
    for (unsigned worker = 0; worker < cfg.workers; ++worker)
        workers.push_back(std::async(std::launch::async,
                                     [&, worker]
                                     {
            std::vector<IndicatorCache> cache(data.signal.size());
            for (size_t k = worker; k < params.size(); k += cfg.workers)
            {
                for (auto &c : cache)
                    c.begin_trial();
                results[k] = evaluate(data, cache, params[k], trade_core::Window{0, split}, false);
                for (auto &c : cache)
                    c.discard_unused();
            }
        }));
    for (auto &worker : workers)
        worker.get();
    size_t winner = params.size();
    for (size_t i = 0; i < results.size(); ++i)
        if (qualifies(results[i], selection) && (winner == params.size() || results[i].score > results[winner].score))
            winner = i;
    nlohmann::json report = {{"status", winner == params.size() ? "no_eligible_candidate" : "complete"},
                             {"trials", params.size()},
                             {"seed", cfg.seed},
                             {"split_timestamp", data.execution[0].timestamp[split]},
                             {"data", data.provenance},
                             {"selection",
                              {{"min_trades", selection.min_trades},
                               {"min_drawdown", selection.min_dd},
                               {"min_gain", selection.min_gain},
                               {"max_gain", selection.max_gain}}},
                             {"invalid_candidates", std::count_if(results.begin(), results.end(),
                                                                  [](const auto &r)
                                                                  {
        return !r.valid;
                              })},
                             {"max_training_trades", std::max_element(results.begin(), results.end(),
                                                                      [](const auto &a, const auto &b)
                                                                      {
        return a.nb_posi_entered < b.nb_posi_entered;
                                                     })->nb_posi_entered}};
    if (winner != params.size())
    {
        std::vector<IndicatorCache> cache(data.signal.size());
        report["parameters"] = params[winner];
        report["training"] = result_json(evaluate(data, cache, params[winner], trade_core::Window{0, split}, true));
        report["holdout"] = result_json(evaluate(data, cache, params[winner], trade_core::Window{split, n}, true));
    }
    return report;
}
template <class Valid, class Evaluate>
int run(backtest_config::StrategyConfig cfg, bool futures, size_t warmup, const std::vector<Axis> &axes, Valid valid,
        Evaluate evaluate, Selection selection = {})
{
    try
    {
        const auto &name = cfg.name;
        if (cfg.is_futures() != futures)
            throw std::runtime_error("Strategy market does not match configuration");
        // The single-asset EMA benchmark intentionally uses the first configured coin.
        if (name == "2EMA_crossover")
            cfg.coins.resize(1);
        backtest_config::print_summary(cfg);
        init_talib();
        auto data = load_market(cfg, warmup);
        auto parameters = sample(axes, cfg.max_trials, cfg.seed, valid);
        const double begin = get_wall_time();
        auto report = evaluate_search(data, cfg, parameters, selection, evaluate);
        report["code_revision"] = BACKTEST_REVISION;
        report["workers"] = cfg.workers;
        report["holdout_fraction"] = cfg.holdout_fraction;
        report["strategy"] = name;
        report["model"] = "next_open_stop_first_v2";
        report["parameter_names"] = nlohmann::json::array();
        for (const auto &a : axes)
            report["parameter_names"].push_back(a.name);
        report["elapsed_seconds"] = get_wall_time() - begin;
        report["resident_mb"] = process_mem_usage();
        report["assumptions"] = {
            {"fee_percent", 0.1},
            {"slippage_bps", 0},
            {"leverage", 1},
            {"liquidation_model", false},
            {"funding_order",
             "exact opening: before orders; fractional milliseconds: after orders, before intrabar stops"},
            {"drawdown", "execution_bar_close"},
            {"intrabar_fill_timestamp", "bar close upper bound"}};
        std::filesystem::create_directories("results");
        const auto stamp = std::chrono::system_clock::now().time_since_epoch().count();
        const std::string path =
            "results/" + name + "-" + std::to_string(stamp) + "-" + std::to_string(getpid()) + ".json";
        std::ofstream f(path + ".tmp");
        f << report.dump(2) << "\n";
        f.close();
        if (!f)
            throw std::runtime_error("Could not write result");
        std::filesystem::rename(path + ".tmp", path);
        std::cout << report["status"] << " | " << parameters.size() << " trials | result " << path << "\n";
        if (report.contains("holdout"))
            std::cout << "Training gain " << report["training"]["gain_percent"] << "% | holdout gain "
                      << report["holdout"]["gain_percent"] << "%\n";
        TA_Shutdown();
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << "\n";
        TA_Shutdown();
        return 1;
    }
}
inline std::vector<int> slots(const backtest_config::StrategyConfig &c)
{
    return integer_range(1, static_cast<int>(c.nb_pairs()));
}
template <class Body> int configure(const std::string &name, Body body)
{
    try
    {
        return body(backtest_config::load(name));
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
template <class Signal>
RUN_RESULTf evaluate(const MarketData &d, trade_core::Window w, bool futures, uint limit, size_t ready, Signal signal,
                     bool trace)
{
    return trade_core::simulate(
        d, w, futures, limit,
        [=](uint p, size_t i)
        {
        return i < ready ? trade_core::Intent{} : signal(p, i);
        },
        trace);
}
} // namespace strategy_runner
