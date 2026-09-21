#include "strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

int main()
{
    return configure("F_SuperReversal_mtf",
                     [](auto cfg)
                     {
        const auto fast = generateRange_int(3, 600, 300);
        const auto slow = generateRange_int(3, 600, 300);
        if (cfg.htf.empty())
            throw std::runtime_error("This strategy requires htf");
        const int ratio =
            backtest_config::timeframe_seconds(cfg.htf) / backtest_config::timeframe_seconds(cfg.timeframe);
        const size_t warm = (std::max(fast.back(), slow.back()) + 1) * ratio + 1;
        return run(
            cfg, true, warm, {{"ema_fast", fast}, {"ema_slow", slow}, {"max_open", slots(cfg)}},
            [](const Params &)
            {
            return true;
            },
            [](const MarketData &d, auto &cache, const Params &p, auto w, bool trace)
            {
            Indicators ind{d, cache};
            const auto f = ind.ema(p[0], true), s = ind.ema(p[1], true), st = ind.supertrend(10, 5.5f, true);
            return evaluate(
                d, w, true, p[2], ind.ready,
                [&d, f, s, st](uint k, size_t i)
                {
                const bool cross = d.signal[k].high[i] > (*f[k])[i] && d.signal[k].low[i] < (*f[k])[i];
                Intent v;
                v.entry = (*f[k])[i] > (*s[k])[i] && (*st[k])[i] == 1 && cross ? 1 : 0;
                if ((*f[k])[i] < (*s[k])[i] && (*st[k])[i] == -1 && cross)
                    v.entry = -1;
                v.exit_long = ((*f[k])[i] < (*s[k])[i] || (*st[k])[i] == -1) && cross;
                v.exit_short = ((*f[k])[i] > (*s[k])[i] || (*st[k])[i] == 1) && cross;
                return v;
                },
                trace);
        },
            {20, -40, -1e6, 1e6});
    });
}
