#pragma once
#include "engine/strategy_runner.hh"

namespace strategies
{
using namespace strategy_runner;
using trade_core::Intent;

// Spot and futures share indicators/exits; their equality rules and grids differ.
struct BigWill
{
    bool futures;
    RUN_RESULTf operator()(const MarketData &d, std::vector<IndicatorCache> &cache, const Params &p,
        trade_core::Window w, bool trace) const
    {
        Indicators ind{d, cache};
        const auto ao = ind.ao(p[0], p[1]), f = ind.ema(p[2]), s = ind.ema(p[3]), r = ind.srsi(), will = ind.will();
        return evaluate(
            d, w, futures, p[9], ind.ready + 1,
            [=](uint k, size_t i)
            {
                const bool trend = futures ? (*f[k])[i] > (*s[k])[i] : (*f[k])[i] >= (*s[k])[i];
                const bool momentum = futures ? (*ao[k])[i] >= 0 : (*ao[k])[i] > 0;
                Intent v;
                v.entry = trend && (*will[k])[i] < p[6] && momentum && (*ao[k])[i - 1] > (*ao[k])[i] ? 1 : 0;
                if (futures && (*f[k])[i] < (*s[k])[i] && (*will[k])[i] > p[7] && (*ao[k])[i] <= 0 &&
                    (*ao[k])[i - 1] < (*ao[k])[i])
                    v.entry = -1;
                v.exit_long = ((*ao[k])[i] < 0 && (*r[k])[i] > p[4]) || (*will[k])[i] > p[7];
                v.exit_short = ((*ao[k])[i] > 0 && (*r[k])[i] < p[5]) || (*will[k])[i] < p[6];
                v.target_fraction = p[8] / 100;
                return v;
            },
            trace);
    }
};

inline int run_big_will(bool futures)
{
    return configure(futures ? "F_BigWill" : "BigWill",
        [=](auto cfg)
        {
            const auto ef = integer_range(2, futures ? 305 : 105, 5), es = integer_range(50, futures ? 610 : 310, 10);
            return run(
                cfg, futures, es.back() + 2,
                {{"ao_fast", integer_range(2, 102, 2)}, {"ao_slow", integer_range(2, 105, 5)}, {"ema_fast", ef},
                    {"ema_slow", es},
                    {"srsi_lower", futures ? float_Nvalues_range(.08f, .92f, 40) : std::vector<float>{.2f}},
                    {"srsi_upper", futures ? float_Nvalues_range(.08f, .92f, 40) : std::vector<float>{.8f}},
                    {"will_lower", futures ? float_Nvalues_range(-1, -100, 50) : std::vector<float>{-85}},
                    {"will_upper", futures ? float_Nvalues_range(-1, -100, 50) : std::vector<float>{-10}},
                    {"target_percent", futures ? float_Nvalues_range(1, 20, 10) : std::vector<float>{15}},
                    {"max_open", slots(cfg)}},
                [=](const Params &p)
                {
                    return futures ? p[1] - p[0] >= 5 && p[3] - p[2] >= 5 && p[5] - p[4] >= .2 && p[7] - p[6] >= 20
                                   : std::abs(p[0] - p[1]) >= 7;
                },
                BigWill{futures}, futures ? Selection{100, -80, 15, 1e6} : Selection{100, -42, 0, 1e6});
        });
}
} // namespace strategies
