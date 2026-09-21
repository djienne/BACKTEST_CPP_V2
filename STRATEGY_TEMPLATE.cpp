#include "strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

// A minimal strategy: define axes, indicator preparation and completed-bar signals.
// The shared runner owns data repair, sampling, execution, holdout and reporting.
int main()
{
    return configure("TEMPLATE_RSI_EMA",
                     [](auto cfg)
                     {
        const auto ema = integer_range(50, 350, 50);
        return run(
            cfg, false, ema.back() + 1,
            {{"rsi", integer_range(7, 29, 7)},
             {"ema", ema},
             {"buy_below", float_Nvalues_range(20, 40, 5)},
             {"sell_above", float_Nvalues_range(60, 80, 5)},
             {"max_open", slots(cfg)}},
            [](const Params &)
            {
            return true;
            },
            [](const MarketData &d, auto &cache, const Params &p, auto w, bool trace)
            {
            Indicators ind{d, cache};
            const auto r = ind.rsi(p[0]), e = ind.ema(p[1]);
            return evaluate(
                d, w, false, p[4], ind.ready,
                [&d, p, r, e](uint k, size_t i)
                {
                Intent s;
                s.entry = (*r[k])[i] < p[2] && d.signal[k].close[i] > (*e[k])[i];
                s.exit_long = (*r[k])[i] > p[3];
                return s;
                },
                trace);
        },
            {100, -50, 0, 1e6});
    });
}
