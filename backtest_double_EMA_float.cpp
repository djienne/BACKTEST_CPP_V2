#include "strategy_runner.hh"
using namespace strategy_runner;
using trade_core::Intent;

RUN_RESULTf evaluate_double_ema(const MarketData &d, std::vector<IndicatorCache> &cache, const Params &p,
                                trade_core::Window w, bool trace)
{
    Indicators ind{d, cache};
    const auto slow = ind.ema(p[0]), fast = ind.ema(p[1]);
    return evaluate(
        d, w, false, 1, ind.ready + 1,
        [=](uint k, size_t i)
        {
        Intent s;
        s.entry = (*fast[k])[i] >= (*slow[k])[i] && (*fast[k])[i - 1] <= (*slow[k])[i - 1];
        s.exit_long = (*fast[k])[i] <= (*slow[k])[i] && (*fast[k])[i - 1] >= (*slow[k])[i - 1];
        return s;
        },
        trace);
}

int main()
{
    return configure("2EMA_crossover",
                     [](auto cfg)
                     {
        cfg.coins.resize(1);
        const auto periods = integer_range(3, 600);
        return run(cfg, false, periods.back() + 1, {{"ema_slow", periods}, {"ema_fast", periods}},
                   [](const Params &p)
                   {
            return p[0] - p[1] >= 7;
        },
                   evaluate_double_ema, {100, -40, -1e6, 10000});
    });
}
