#pragma once
#include "engine/strategy_runner.hh"

namespace strategies
{
using namespace strategy_runner;
using trade_core::Intent;

struct BBTrend
{
    bool futures;
    RUN_RESULTf operator()(const MarketData &d, std::vector<IndicatorCache> &cache, const Params &p,
        trade_core::Window w, bool trace) const
    {
        Indicators ind{d, cache};
        const auto e = ind.ema(p[0]);
        const auto bands = ind.bands(p[1], p[2]);
        return evaluate(
            d, w, futures, p[3], ind.ready + 1,
            [&d, e, bands, this](uint k, size_t i)
            {
                const auto &c = d.signal[k].close, &up = *bands[0][k], &mid = *bands[1][k], &low = *bands[2][k];
                Intent s;
                s.entry = c[i - 1] < up[i - 1] && c[i] > up[i] && c[i] > (*e[k])[i] ? 1 : 0;
                if (futures && c[i - 1] > low[i - 1] && c[i] < low[i] && c[i] < (*e[k])[i])
                    s.entry = -1;
                s.exit_long = c[i] < mid[i];
                s.exit_short = c[i] > mid[i];
                return s;
            },
            trace);
    }
};

inline int run_bbtrend(bool futures)
{
    return configure(futures ? "F_BBTREND" : "BBTREND",
        [=](auto cfg)
        {
            const auto ema = futures ? generateRange_int(5, 602, 120) : integer_range(5, 602, 5);
            const auto length = futures ? generateRange_int(5, 400, 120) : integer_range(5, 250, 5);
            return run(
                cfg, futures, std::max(ema.back(), length.back()) + 2,
                {{"ema", ema}, {"length", length},
                    {"deviation", float_Nvalues_range(futures ? 1.0f : .25f, 4.5f, futures ? 40 : 10)},
                    {"max_open", slots(cfg)}},
                [](const Params &)
                {
                    return true;
                },
                BBTrend{futures}, futures ? Selection{100, -50, -1e6, 1e6} : Selection{100, -40, 0, 1e6});
        });
}
} // namespace strategies
