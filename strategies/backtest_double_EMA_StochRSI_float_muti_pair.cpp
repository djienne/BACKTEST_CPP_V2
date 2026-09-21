#include "engine/strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

int main()
{
    return configure("2EMA_crossover_StochRSI",
        [](auto cfg)
        {
            const auto fast = integer_range(2, 304, 1), slow = integer_range(70, 604, 1);
            // Preserve the original comparator orientation: the longer-range EMA is
            // on the left. This strategy has always allocated one slot per pair.
            return run(
                cfg, false, slow.back() + 1, {{"ema_left", slow}, {"ema_right", fast}},
                [](const Params &)
                {
                    return true;
                },
                [](const MarketData &d, auto &cache, const Params &p, auto w, bool trace)
                {
                    Indicators ind{d, cache};
                    const auto f = ind.ema(p[0]), s = ind.ema(p[1]), r = ind.srsi();
                    return evaluate(
                        d, w, false, static_cast<uint>(d.signal.size()), ind.ready,
                        [=](uint k, size_t i)
                        {
                            Intent v;
                            v.entry = (*f[k])[i] >= (*s[k])[i] && (*r[k])[i] < 0.8f;
                            v.exit_long = (*f[k])[i] <= (*s[k])[i] && (*r[k])[i] > 0.2f;
                            return v;
                        },
                        trace);
                },
                {100, -50, 0, 1e6});
        });
}
