#include "strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

int main()
{
    return configure("F_BBTREND",
                     [](auto cfg)
                     {
        const auto ema = generateRange_int(5, 602, 120);
        const auto length = generateRange_int(5, 400, 120);
        return run(
            cfg, true, std::max(ema.back(), length.back()) + 2,
            {{"ema", ema},
             {"length", length},
             {"deviation", float_Nvalues_range(1.0, 4.5, 40)},
             {"max_open", slots(cfg)}},
            [](const Params &)
            {
            return true;
            },
            [](const MarketData &d, auto &cache, const Params &p, auto w, bool trace)
            {
            Indicators ind{d, cache};
            const auto e = ind.ema(p[0]);
            const auto bands = ind.bands(p[1], p[2]);
            return evaluate(
                d, w, true, p[3], ind.ready + 1,
                [&d, e, bands](uint k, size_t i)
                {
                const auto &c = d.signal[k].close, &up = *bands[0][k], &mid = *bands[1][k], &low = *bands[2][k];
                Intent s;
                s.entry = c[i - 1] < up[i - 1] && c[i] > up[i] && c[i] > (*e[k])[i] ? 1 : 0;
                if (c[i - 1] > low[i - 1] && c[i] < low[i] && c[i] < (*e[k])[i])
                    s.entry = -1;
                s.exit_long = c[i] < mid[i];
                s.exit_short = c[i] > mid[i];
                return s;
                },
                trace);
        },
            {100, -50, -1e6, 1e6});
    });
}
