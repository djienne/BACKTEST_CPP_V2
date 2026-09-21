#include "strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

int main()
{
    return configure("TRIX",
                     [](auto cfg)
                     {
        const auto ema_periods = integer_range(40, 602, 3), length = integer_range(2, 100, 2),
                   smooth = integer_range(10, 100, 2);
        const size_t warm = std::max(ema_periods.back(), 3 * (length.back() - 1) + smooth.back()) + 1;
        return run(
            cfg, false, warm,
            {{"ema", ema_periods}, {"trix_length", length}, {"trix_signal", smooth}, {"max_open", slots(cfg)}},
            [](const Params &)
            {
            return true;
            },
            [](const MarketData &d, auto &cache, const Params &p, auto w, bool trace)
            {
            Indicators ind{d, cache};
            const auto ema = ind.ema(p[0]), t = ind.trix(p[1], p[2]), r = ind.srsi();
            return evaluate(
                d, w, false, p[3], ind.ready,
                [&d, ema, t, r](uint k, size_t i)
                {
                Intent s;
                s.entry = d.signal[k].close[i] > (*ema[k])[i] && (*t[k])[i] > 0 && (*r[k])[i] < 0.8f;
                s.exit_long = (*t[k])[i] < 0 && (*r[k])[i] > 0.2f;
                return s;
                },
                trace);
        },
            {200, -40, 0, 1e6});
    });
}
