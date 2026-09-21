#pragma once
#include "engine/strategy_runner.hh"

namespace strategies
{
using namespace strategy_runner;
using trade_core::Intent;

struct EMA3
{
    bool futures;
    RUN_RESULTf operator()(const MarketData &d, std::vector<IndicatorCache> &cache, const Params &p,
        trade_core::Window w, bool trace) const
    {
        Indicators ind{d, cache};
        const auto e1 = ind.ema(p[0]), e2 = ind.ema(p[1]), e3 = ind.ema(p[2]), a = ind.atr(), k = ind.srsi(1),
                   v = ind.srsi(2);
        return evaluate(
            d, w, futures, p[7], ind.ready + 1,
            [&d, p, e1, e2, e3, a, k, v, this](uint pair, size_t i)
            {
                const auto &f = *e1[pair], &m = *e2[pair], &s = *e3[pair], &sk = *k[pair], &sd = *v[pair];
                Intent intent;
                // Preserve the model-specific equality rules rather than silently unifying them.
                if (futures)
                {
                    intent.entry = f[i] > m[i] && m[i] > s[i] && d.signal[pair].close[i] > f[i] && sk[i] <= p[5] &&
                                           sd[i] <= p[5] && sk[i - 1] >= sd[i - 1] && sk[i] <= sd[i]
                                       ? 1
                                       : 0;
                    if (f[i] < m[i] && m[i] < s[i] && d.signal[pair].close[i] < f[i] && sk[i] >= p[6] &&
                        sd[i] >= p[6] && sk[i - 1] <= sd[i - 1] && sk[i] >= sd[i])
                        intent.entry = -1;
                }
                else
                {
                    intent.entry = f[i] >= m[i] && m[i] >= s[i] && d.signal[pair].close[i] >= f[i] && sk[i] < p[5] &&
                                   sd[i] < p[5] && sk[i - 1] > sd[i - 1] && sk[i] <= sd[i];
                    intent.target_fraction = .15;
                }
                intent.stop_distance = p[4] * (*a[pair])[i];
                intent.target_distance = p[3] * (*a[pair])[i];
                intent.max_hold_seconds = 2 * 86400;
                return intent;
            },
            trace);
    }
};

inline int run_ema3(bool futures)
{
    return configure(futures ? "F_EMA3_SRSI_ATR" : "EMA3_SRSI_ATR",
        [=](auto cfg)
        {
            const auto e1 = futures ? generateRange_int(3, 100, 30) : integer_range(3, 100, 5);
            const auto e2 = futures ? generateRange_int(5, 400, 120) : integer_range(5, 400, 10);
            const auto e3 = futures ? generateRange_int(5, 580, 250) : integer_range(5, 580, 10);
            const auto multipliers = futures ? float_Nvalues_range(2, 20, 10) : float_Nvalues_range(3, 12, 10);
            return run(
                cfg, futures, e3.back() + 2,
                {{"ema1", e1}, {"ema2", e2}, {"ema3", e3}, {"target_atr", multipliers}, {"stop_atr", multipliers},
                    {"srsi_lower", futures ? float_Nvalues_range(.09f, .91f, 12) : float_Nvalues_range(.1f, .9f, 9)},
                    {"srsi_upper", futures ? float_Nvalues_range(.09f, .91f, 12) : std::vector<float>{0}},
                    {"max_open", slots(cfg)}},
                [](const Params &p)
                {
                    return p[1] - p[0] >= 5 && p[2] - p[1] >= 5;
                },
                EMA3{futures}, futures ? Selection{100, -50, -1e6, 1e6} : Selection{300, -80, 500, 1e6});
        });
}
} // namespace strategies
