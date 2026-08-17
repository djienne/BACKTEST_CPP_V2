#include "config.hh"

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "nlohmann/json.hpp"
#include "tools_fatal.hh"

using json = nlohmann::json;

namespace backtest_config
{
namespace
{

std::string config_path()
{
    if (const char *env = std::getenv("BACKTEST_CONFIG"))
    {
        return std::string(env);
    }
    return "./backtest_config.json";
}

std::string join_path(const std::string &dir, const std::string &leaf)
{
    if (dir.empty())
    {
        return leaf;
    }
    return dir.back() == '/' ? dir + leaf : dir + "/" + leaf;
}

} // namespace

StrategyConfig load(const std::string &strategy_name)
{
    const std::string path = config_path();

    std::ifstream file(path);
    if (!file.is_open())
    {
        BACKTEST_FATAL("Cannot open config '" + path +
                       "'. Run from the repo root, or set BACKTEST_CONFIG to its location.");
    }

    json root;
    try
    {
        root = json::parse(file, nullptr, true, /*ignore_comments=*/true);
    }
    catch (const json::parse_error &err)
    {
        BACKTEST_FATAL("Malformed config '" + path + "': " + err.what());
    }

    StrategyConfig cfg{};
    cfg.name = strategy_name;
    cfg.data_dir = root.value("data_dir", std::string("./data/data"));

    if (!root.contains("coins") || !root["coins"].is_array() || root["coins"].empty())
    {
        BACKTEST_FATAL("Config '" + path + "' has no non-empty \"coins\" array.");
    }
    cfg.coins = root["coins"].get<std::vector<std::string>>();

    if (cfg.coins.size() > MAX_PAIRS)
    {
        BACKTEST_FATAL("Config lists " + std::to_string(cfg.coins.size()) + " coins but MAX_PAIRS is " +
                       std::to_string(MAX_PAIRS) + ". Raise MAX_PAIRS in config.hh and rebuild.");
    }

    const auto strategies = root.find("strategies");
    if (strategies == root.end() || !strategies->is_object())
    {
        BACKTEST_FATAL("Config '" + path + "' has no \"strategies\" object.");
    }

    const auto entry = strategies->find(strategy_name);
    if (entry == strategies->end())
    {
        // Deliberately fatal rather than defaulting: a strategy quietly guessing its
        // market and timeframe is how a binary ends up reading files nobody downloads.
        std::string known;
        for (auto it = strategies->begin(); it != strategies->end(); ++it)
        {
            known += "\n    " + it.key();
        }
        BACKTEST_FATAL("Config '" + path + "' has no entry for strategy \"" + strategy_name +
                       "\". Known strategies:" + known);
    }

    cfg.market = entry->value("market", std::string("spot"));
    cfg.timeframe = entry->value("timeframe", std::string(""));
    cfg.htf = entry->value("htf", std::string(""));

    if (cfg.market != "spot" && cfg.market != "futures")
    {
        BACKTEST_FATAL("Strategy \"" + strategy_name + "\" has market \"" + cfg.market +
                       "\"; expected \"spot\" or \"futures\".");
    }
    if (cfg.timeframe.empty())
    {
        BACKTEST_FATAL("Strategy \"" + strategy_name + "\" has no \"timeframe\".");
    }

    return cfg;
}

std::vector<std::string> spot_paths(const StrategyConfig &cfg, const std::string &timeframe)
{
    const std::string tf = timeframe.empty() ? cfg.timeframe : timeframe;
    std::vector<std::string> out;
    out.reserve(cfg.coins.size());
    for (const std::string &coin : cfg.coins)
    {
        out.push_back(join_path(cfg.data_dir, "binance/" + tf + "/" + coin + "-USDT.csv"));
    }
    return out;
}

std::vector<std::string> futures_paths(const StrategyConfig &cfg, const std::string &timeframe)
{
    const std::string tf = timeframe.empty() ? cfg.timeframe : timeframe;
    std::vector<std::string> out;
    out.reserve(cfg.coins.size());
    for (const std::string &coin : cfg.coins)
    {
        out.push_back(join_path(cfg.data_dir, "futures/" + coin + "_USDT-" + tf + "-futures.json"));
    }
    return out;
}

std::vector<std::string> funding_paths(const StrategyConfig &cfg)
{
    std::vector<std::string> out;
    out.reserve(cfg.coins.size());
    for (const std::string &coin : cfg.coins)
    {
        out.push_back(join_path(cfg.data_dir, "futures/" + coin + "_USDT-8h-funding_rate.json"));
    }
    return out;
}

std::vector<std::string> data_paths(const StrategyConfig &cfg, const std::string &timeframe)
{
    return cfg.is_futures() ? futures_paths(cfg, timeframe) : spot_paths(cfg, timeframe);
}

void print_summary(const StrategyConfig &cfg)
{
    std::cout << "Config           : " << config_path() << std::endl;
    std::cout << "Market           : " << cfg.market << std::endl;
    std::cout << "Timeframe        : " << cfg.timeframe;
    if (!cfg.htf.empty())
    {
        std::cout << "  (higher: " << cfg.htf << ")";
    }
    std::cout << std::endl;
    std::cout << "Coins (" << cfg.nb_pairs() << ")       : ";
    for (size_t i = 0; i < cfg.coins.size(); ++i)
    {
        std::cout << cfg.coins[i] << (i + 1 < cfg.coins.size() ? ", " : "");
    }
    std::cout << std::endl;
}

} // namespace backtest_config
