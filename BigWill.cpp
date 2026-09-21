#include "strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

int main()
{
    return configure("BigWill",
                     [](auto cfg)
                     {
        const auto ef = integer_range(2, 105, 5), es = integer_range(50, 310, 10);
        return run(
            cfg, false, es.back() + 2,
            {{"ao_fast", integer_range(2, 102, 2)},
             {"ao_slow", integer_range(2, 105, 5)},
             {"ema_fast", ef},
             {"ema_slow", es},
             {"srsi_lower", std::vector<float>{0.2f}},
             {"srsi_upper", std::vector<float>{0.8f}},
             {"will_lower", std::vector<int>{-85}},
             {"will_upper", std::vector<int>{-10}},
             {"target_percent", std::vector<int>{15}},
             {"max_open", slots(cfg)}},
            [](const Params &p)
            {
            return std::abs(p[0] - p[1]) >= 7;
            },
            [](const MarketData &d, auto &cache, const Params &p, auto w, bool trace)
            {
            Indicators ind{d, cache};
            const auto ao = ind.ao(p[0], p[1]), f = ind.ema(p[2]), s = ind.ema(p[3]), r = ind.srsi(), will = ind.will();
            return evaluate(
                d, w, false, p[9], ind.ready + 1,
                [=](uint k, size_t i)
                {
                Intent v;
                v.entry =
                    (*f[k])[i] >= (*s[k])[i] && (*will[k])[i] < p[6] && (*ao[k])[i] > 0 && (*ao[k])[i - 1] > (*ao[k])[i]
                        ? 1
                        : 0;

                v.exit_long = ((*ao[k])[i] < 0 && (*r[k])[i] > p[4]) || (*will[k])[i] > p[7];
                v.exit_short = ((*ao[k])[i] > 0 && (*r[k])[i] < p[5]) || (*will[k])[i] < p[6];
                v.target_fraction = p[8] / 100;
                return v;
                },
                trace);
        },
            {100, -42, 0, 1e6});
    });
}
