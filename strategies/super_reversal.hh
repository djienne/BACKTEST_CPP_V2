#pragma once
#include "engine/strategy_runner.hh"

namespace strategies
{
using namespace strategy_runner;
using trade_core::Intent;

struct SuperReversal
{
    bool futures;
    RUN_RESULTf operator()(const MarketData &d, std::vector<IndicatorCache> &cache, const Params &p,
        trade_core::Window w, bool trace) const
    {
        Indicators ind{d, cache};
        const auto f = ind.ema(p[0], true), s = ind.ema(p[1], true),
                   st = ind.supertrend(futures ? 10 : 15, futures ? 5.5f : 5.0f, true);
        return evaluate(
            d, w, futures, p[2], ind.ready,
            [&d, f, s, st, this](uint k, size_t i)
            {
                const bool cross = d.signal[k].high[i] > (*f[k])[i] && d.signal[k].low[i] < (*f[k])[i];
                Intent v;
                v.entry = (*f[k])[i] > (*s[k])[i] && (*st[k])[i] == 1 && cross ? 1 : 0;
                if (futures && (*f[k])[i] < (*s[k])[i] && (*st[k])[i] == -1 && cross)
                    v.entry = -1;
                v.exit_long = ((*f[k])[i] < (*s[k])[i] || (*st[k])[i] == -1) && cross;
                v.exit_short = ((*f[k])[i] > (*s[k])[i] || (*st[k])[i] == 1) && cross;
                return v;
            },
            trace);
    }
};

inline int run_super_reversal(bool futures)
{
    return configure(futures ? "F_SuperReversal_mtf" : "SuperReversal_mtf",
        [=](auto cfg)
        {
            const auto fast = futures ? generateRange_int(3, 600, 300) : integer_range(3, 204, 2);
            const auto slow = futures ? generateRange_int(3, 600, 300) : integer_range(70, 590, 5);
            if (cfg.htf.empty())
                throw std::runtime_error("This strategy requires htf");
            const int ratio =
                backtest_config::timeframe_seconds(cfg.htf) / backtest_config::timeframe_seconds(cfg.timeframe);
            const size_t warm = (std::max(fast.back(), slow.back()) + 1) * ratio + 1;
            return run(
                cfg, futures, warm, {{"ema_fast", fast}, {"ema_slow", slow}, {"max_open", slots(cfg)}},
                [](const Params &)
                {
                    return true;
                },
                SuperReversal{futures}, futures ? Selection{20, -40, -1e6, 1e6} : Selection{200, -40, 0, 1e6});
        });
}
} // namespace strategies
