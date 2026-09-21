#include "double_ema.hh"
using namespace strategy_runner;

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
                strategies::evaluate_double_ema, {100, -40, -1e6, 10000});
        });
}
