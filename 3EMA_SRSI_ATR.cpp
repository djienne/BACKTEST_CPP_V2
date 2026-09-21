#include "strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

int main()
{
    return configure("EMA3_SRSI_ATR",
                     [](auto cfg)
                     {
        const auto range_e1 = integer_range(3, 100, 5);
        const auto range_e2 = integer_range(5, 400, 10);
        const auto range_e3 = integer_range(5, 580, 10);
        return run(
            cfg, false, range_e3.back() + 2,
            {{"ema1", range_e1},
             {"ema2", range_e2},
             {"ema3", range_e3},
             {"target_atr", integer_range(3, 12)},
             {"stop_atr", integer_range(3, 12)},
             {"srsi_lower", float_Nvalues_range(0.1f, 0.9f, 9)},
             {"srsi_upper", std::vector<int>{0}},
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
                d, w, false, p[7], ind.ready + 1,
                [&d, p, e1, e2, e3, a, k, v](uint pair, size_t i)
                {
                const auto &f = *e1[pair], &m = *e2[pair], &s = *e3[pair], &sk = *k[pair], &sd = *v[pair];
                Intent intent;
                intent.entry = f[i] >= m[i] && m[i] >= s[i] && d.signal[pair].close[i] >= f[i] && sk[i] < p[5] &&
                               sd[i] < p[5] && sk[i - 1] > sd[i - 1] && sk[i] <= sd[i];
                intent.target_fraction = 0.15;
                intent.stop_distance = p[4] * (*a[pair])[i];
                intent.target_distance = p[3] * (*a[pair])[i];
                intent.max_hold_seconds = 2 * 86400;
                return intent;
                },
                trace);
        },
            {300, -80, 500, 1e6});
    });
}
