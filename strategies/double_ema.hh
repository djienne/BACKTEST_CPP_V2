#pragma once
#include "engine/strategy_runner.hh"

namespace strategies
{
using namespace strategy_runner;
using trade_core::Intent;
inline RUN_RESULTf evaluate_double_ema(
    const MarketData &d, std::vector<IndicatorCache> &cache, const Params &p, trade_core::Window w, bool trace)
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
} // namespace strategies
