#include "engine/strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

int main()
{
    return configure("SuperTrend_EMA_ATR",
        [](auto cfg)
        {
            const auto ema = integer_range(5, 602, 1), multipliers = integer_range(1, 11);
            return run(
                cfg, false, ema.back() + 1,
                {{"ema", ema}, {"target_atr", multipliers}, {"stop_atr", multipliers}, {"max_open", slots(cfg)}},
                [](const Params &)
                {
                    return true;
                },
                [](const MarketData &d, auto &cache, const Params &p, auto w, bool trace)
                {
                    Indicators ind{d, cache};
                    const auto e = ind.ema(p[0]), a = ind.atr(17), st = ind.supertrend(10, 3);
                    return evaluate(
                        d, w, false, p[3], ind.ready,
                        [&d, p, e, a, st](uint k, size_t i)
                        {
                            Intent s;
                            s.entry = d.signal[k].close[i] > (*e[k])[i] && (*st[k])[i] == 1;
                            s.exit_long = d.signal[k].low[i] < (*e[k])[i] || (*st[k])[i] == -1;
                            s.stop_distance = p[2] * (*a[k])[i];
                            s.target_distance = p[1] * (*a[k])[i];
                            s.trailing = true;
                            return s;
                        },
                        trace);
                });
        });
}
