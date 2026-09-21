#include "strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

int main()
{
    return configure("F_BigWill",
                     [](auto cfg)
                     {
        const auto ef = integer_range(2, 305, 5), es = integer_range(50, 610, 10);
        return run(
            cfg, true, es.back() + 2,
            {{"ao_fast", integer_range(2, 102, 2)},
             {"ao_slow", integer_range(2, 105, 5)},
             {"ema_fast", ef},
             {"ema_slow", es},
             {"srsi_lower", float_Nvalues_range(0.08f, 0.92f, 40)},
             {"srsi_upper", float_Nvalues_range(0.08f, 0.92f, 40)},
             {"will_lower", float_Nvalues_range(-1, -100, 50)},
             {"will_upper", float_Nvalues_range(-1, -100, 50)},
             {"target_percent", float_Nvalues_range(1, 20, 10)},
             {"max_open", slots(cfg)}},
            [](const Params &p)
            {
            return p[1] - p[0] >= 5 && p[3] - p[2] >= 5 && p[5] - p[4] >= 0.2 && p[7] - p[6] >= 20;
            },
            [](const MarketData &d, auto &cache, const Params &p, auto w, bool trace)
            {
            Indicators ind{d, cache};
            const auto ao = ind.ao(p[0], p[1]), f = ind.ema(p[2]), s = ind.ema(p[3]), r = ind.srsi(), will = ind.will();
            return evaluate(
                d, w, true, p[9], ind.ready + 1,
                [=](uint k, size_t i)
                {
                Intent v;
                v.entry =
                    (*f[k])[i] > (*s[k])[i] && (*will[k])[i] < p[6] && (*ao[k])[i] >= 0 && (*ao[k])[i - 1] > (*ao[k])[i]
                        ? 1
                        : 0;
                if ((*f[k])[i] < (*s[k])[i] && (*will[k])[i] > p[7] && (*ao[k])[i] <= 0 &&
                    (*ao[k])[i - 1] < (*ao[k])[i])
                    v.entry = -1;
                v.exit_long = ((*ao[k])[i] < 0 && (*r[k])[i] > p[4]) || (*will[k])[i] > p[7];
                v.exit_short = ((*ao[k])[i] > 0 && (*r[k])[i] < p[5]) || (*will[k])[i] < p[6];
                v.target_fraction = p[8] / 100;
                return v;
                },
                trace);
        },
            {100, -80, 15, 1e6});
    });
}
