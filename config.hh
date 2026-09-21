#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace backtest_config
{
constexpr unsigned MAX_PAIRS = 16;
struct StrategyConfig
{
    std::string name, market, timeframe, htf, data_dir, path;
    std::string start, end; // UTC ISO dates/times; end is exclusive.
    uint64_t max_trials = 1000, seed = 42;
    unsigned workers = 1;
    double holdout_fraction = 0.20;
    bool offline = false;
    unsigned nb_pairs() const
    {
        return static_cast<unsigned>(coins.size());
    }
    bool is_futures() const
    {
        return market == "futures";
    }
    std::vector<std::string> coins;
};
int timeframe_seconds(const std::string &value);
StrategyConfig load(const std::string &name);
void print_summary(const StrategyConfig &);
} // namespace backtest_config
