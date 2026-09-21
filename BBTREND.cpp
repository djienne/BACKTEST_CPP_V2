#include "strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

int main()
{
    return configure("BBTREND",
                     [](auto cfg)
                     {
        const auto ema = integer_range(5, 602, 5);
        const auto length = integer_range(5, 250, 5);
        return run(
            cfg, false, std::max(ema.back(), length.back()) + 2,
            {{"ema", ema},
             {"length", length},
             {"deviation", float_Nvalues_range(0.25, 4.5, 10)},
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
                d, w, false, p[3], ind.ready + 1,
                [&d, e, bands](uint k, size_t i)
                {
                const auto &c = d.signal[k].close, &up = *bands[0][k], &mid = *bands[1][k], &low = *bands[2][k];
                Intent s;
                s.entry = c[i - 1] < up[i - 1] && c[i] > up[i] && c[i] > (*e[k])[i] ? 1 : 0;
                (void)low;
                s.exit_long = c[i] < mid[i];
                s.exit_short = c[i] > mid[i];
                return s;
                },
                trace);
        },
            {100, -40, 0, 1e6});
    });
}
