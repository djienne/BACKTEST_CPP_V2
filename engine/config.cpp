#include "config.hh"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <nlohmann/json.hpp>

namespace backtest_config
{
int timeframe_seconds(const std::string &v)
{
    static const std::vector<std::pair<std::string, int>> supported = {
        {"1m", 60},   {"3m", 180},   {"5m", 300},   {"15m", 900},  {"30m", 1800},  {"1h", 3600},
        {"2h", 7200}, {"4h", 14400}, {"6h", 21600}, {"8h", 28800}, {"12h", 43200}, {"1d", 86400}};
    for (const auto &item : supported)
        if (item.first == v)
            return item.second;
    throw std::runtime_error("Unsupported timeframe: " + v);
}
StrategyConfig load(const std::string &name)
{
    StrategyConfig c;
    c.name = name;
    c.path = std::getenv("BACKTEST_CONFIG") ? std::getenv("BACKTEST_CONFIG") : "backtest_config.json";
    std::ifstream f(c.path);
    if (!f)
        throw std::runtime_error("Cannot open config: " + c.path);
    const auto j = nlohmann::json::parse(f);
    c.data_dir = j.value("data_dir", std::string("data/research"));
    c.coins = j.at("coins").get<std::vector<std::string>>();
    std::set<std::string> seen;
    for (const auto &coin : c.coins)
    {
        if (coin.empty() || coin.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") != std::string::npos ||
            !seen.insert(coin).second)
            throw std::runtime_error("Invalid or repeated coin: " + coin);
    }
    if (c.coins.empty() || c.coins.size() > MAX_PAIRS)
        throw std::runtime_error("Expected 1..16 unique coins");
    const auto &s = j.at("strategies").at(name);
    c.market = s.at("market").get<std::string>();
    c.timeframe = s.at("timeframe").get<std::string>();
    c.htf = s.value("htf", std::string());
    if (c.market != "spot" && c.market != "futures")
        throw std::runtime_error("Invalid market");
    const int tf = timeframe_seconds(c.timeframe);
    if (!c.htf.empty() && (timeframe_seconds(c.htf) <= tf || timeframe_seconds(c.htf) % tf))
        throw std::runtime_error("Higher timeframe must be a larger integral multiple");
    const auto r = j.value("run", nlohmann::json::object());
    c.start = r.value("start", std::string());
    c.end = r.value("end", std::string());
    const int64_t trials = r.value("max_trials", int64_t(1000));
    const int64_t seed = r.value("seed", int64_t(42));
    const int workers = r.value("workers", 1);
    if (trials < 0 || seed < 0 || workers < 1 || workers > 16)
        throw std::runtime_error("Invalid search settings");
    c.max_trials = static_cast<uint64_t>(trials);
    c.seed = static_cast<uint64_t>(seed);
    c.workers = static_cast<unsigned>(workers);
    c.holdout_fraction = r.value("holdout_fraction", 0.20);
    if (!(c.holdout_fraction > 0 && c.holdout_fraction < 1))
        throw std::runtime_error("Holdout fraction must be between 0 and 1");
    c.offline = std::getenv("BACKTEST_OFFLINE") && std::string(std::getenv("BACKTEST_OFFLINE")) == "1";
    return c;
}
void print_summary(const StrategyConfig &c)
{
    std::cout << c.name << " | " << c.market << " " << c.timeframe << " | seed " << c.seed << " | trials "
              << c.max_trials << "\n";
}
} // namespace backtest_config
