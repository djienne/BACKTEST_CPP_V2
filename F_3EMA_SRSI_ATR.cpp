#include "strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

int main()
{
    return configure("F_EMA3_SRSI_ATR",
                     [](auto cfg)
                     {
        const auto range_e1 = generateRange_int(3, 100, 30);
        const auto range_e2 = generateRange_int(5, 400, 120);
        const auto range_e3 = generateRange_int(5, 580, 250);
        return run(
            cfg, true, range_e3.back() + 2,
            {{"ema1", range_e1},
             {"ema2", range_e2},
             {"ema3", range_e3},
             {"target_atr", float_Nvalues_range(2, 20, 10)},
             {"stop_atr", float_Nvalues_range(2, 20, 10)},
             {"srsi_lower", float_Nvalues_range(0.09f, 0.91f, 12)},
             {"srsi_upper", float_Nvalues_range(0.09f, 0.91f, 12)},
             {"max_open", slots(cfg)}},
            [](const Params &p)
            {
            return p[1] - p[0] >= 5 && p[2] - p[1] >= 5;
            },
            [](const MarketData &d, auto &cache, const Params &p, auto w, bool trace)
            {
            Indicators ind{d, cache};
            const auto e1 = ind.ema(p[0]), e2 = ind.ema(p[1]), e3 = ind.ema(p[2]), a = ind.atr(), k = ind.srsi(1),
                       v = ind.srsi(2);
            return evaluate(
                d, w, true, p[7], ind.ready + 1,
                [&d, p, e1, e2, e3, a, k, v](uint pair, size_t i)
                {
                const auto &f = *e1[pair], &m = *e2[pair], &s = *e3[pair], &sk = *k[pair], &sd = *v[pair];
                Intent intent;
                intent.entry = f[i] > m[i] && m[i] > s[i] && d.signal[pair].close[i] > f[i] && sk[i] <= p[5] &&
                                       sd[i] <= p[5] && sk[i - 1] >= sd[i - 1] && sk[i] <= sd[i]
                                   ? 1
                                   : 0;
                if (f[i] < m[i] && m[i] < s[i] && d.signal[pair].close[i] < f[i] && sk[i] >= p[6] && sd[i] >= p[6] &&
                    sk[i - 1] <= sd[i - 1] && sk[i] >= sd[i])
                    intent.entry = -1;
                intent.stop_distance = p[4] * (*a[pair])[i];
                intent.target_distance = p[3] * (*a[pair])[i];
                intent.max_hold_seconds = 2 * 86400;
                return intent;
                },
                trace);
        },
            {100, -50, -1e6, 1e6});
    });
}
