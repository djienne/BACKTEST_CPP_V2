#pragma once

#include <string>
#include <vector>

#include "Klinef.hh" // uint

// ---------------------------------------------------------------------------
//  Strategy configuration, read from backtest_config.json at startup.
//
//  The same file drives tools/download_data.py, so the data a strategy opens and
//  the data the downloader fetches cannot drift apart -- which is how the repo
//  ended up shipping strategies that referenced symbols with no data, and symbols
//  (MATIC, XMR, EOS) that Binance no longer lists at all.
//
//  Path resolution: $BACKTEST_CONFIG if set, else ./backtest_config.json relative
//  to the working directory. Strategies are already run from the repo root because
//  their data paths are relative, so the default works unchanged.
// ---------------------------------------------------------------------------

namespace backtest_config
{

// Compile-time ceiling for the per-pair arrays in the hot loop.
//
// The traded pair count is a runtime value (see StrategyConfig::nb_pairs), but
// trade_core::PortfolioState<N> and the std::array<float, N> price buffers are
// fixed-size on purpose: PROCESS runs millions of times per sweep, and a
// std::vector there would put a heap allocation in that path. Sizing to a ceiling
// and iterating only nb_pairs of it costs a few unused slots and nothing else.
constexpr uint MAX_PAIRS = 16;

struct StrategyConfig
{
    std::string name;
    std::string market;    // "spot" or "futures"
    std::string timeframe; // the timeframe the strategy trades on
    std::string htf;       // higher timeframe, empty unless multi-timeframe
    std::vector<std::string> coins;
    std::string data_dir;

    uint nb_pairs() const { return static_cast<uint>(coins.size()); }
    bool is_futures() const { return market == "futures"; }
};

// Returns the entry for `strategy_name`. Aborts with BACKTEST_FATAL if the config
// file is missing, malformed, names no coins, lists more than MAX_PAIRS of them, or
// has no entry for this strategy -- a strategy silently falling back to a guessed
// universe is exactly the drift this file exists to prevent.
StrategyConfig load(const std::string &strategy_name);

// Data file paths for the configured coins, in coin order.
//   spot     <data_dir>/binance/<tf>/<COIN>-USDT.csv
//   futures  <data_dir>/futures/<COIN>_USDT-<tf>-futures.json
//   funding  <data_dir>/futures/<COIN>_USDT-8h-funding_rate.json
std::vector<std::string> spot_paths(const StrategyConfig &cfg, const std::string &timeframe = "");
std::vector<std::string> futures_paths(const StrategyConfig &cfg, const std::string &timeframe = "");
std::vector<std::string> funding_paths(const StrategyConfig &cfg);

// Convenience: spot_paths or futures_paths according to cfg.market.
std::vector<std::string> data_paths(const StrategyConfig &cfg, const std::string &timeframe = "");

// Prints the resolved universe, so a run's log states which coins and timeframe it
// actually used rather than leaving it to be inferred from the file names.
void print_summary(const StrategyConfig &cfg);

} // namespace backtest_config
