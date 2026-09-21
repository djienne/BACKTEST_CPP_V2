#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "tools_fatal.hh"

// Worker-local indicator series and their first usable indices. Keep references
// outside the bar loop; discard unreferenced parameter series after each trial.
class IndicatorCache
{
public:
    // "EMA:200", "STOCHRSI:14:14", "SUPERTREND:10:3.000000" -- readable in a debugger
    // Floating parameters use six decimals. Current grids are coarser; increase
    // key precision before introducing steps smaller than 1e-6.
    static std::string key(const std::string &name)
    {
        return name;
    }

    template <typename T, typename... Rest>
    static std::string key(const std::string &name, const T &first, const Rest &...rest)
    {
        return key(name + ":" + std::to_string(first), rest...);
    }

    void put(const std::string &k, std::vector<float> series, size_t first_valid = 0)
    {
        series_by_key_[k] = std::move(series);
        first_valid_[k] = first_valid;
    }

    size_t first_valid(const std::string &k) const
    {
        return first_valid_.at(k);
    }

    bool has(const std::string &k) const
    {
        return series_by_key_.find(k) != series_by_key_.end();
    }

    // Aborts on a missing key rather than returning an empty series: an empty vector
    // would make every comparison against it read out of bounds or silently compare
    // against nothing, which is exactly the class of bug this cache exists to prevent.
    const std::vector<float> &get(const std::string &k) const
    {
        const auto it = series_by_key_.find(k);
        if (it == series_by_key_.end())
        {
            BACKTEST_FATAL("IndicatorCache: no series for key '" + k + "'. Compute it before evaluation.");
        }
        used_.insert(k);
        return it->second;
    }

    void begin_trial()
    {
        used_.clear();
    }
    // Retain only the current parameter set. Common series survive across trials;
    // varying series use at most previous + current trial memory during preparation.
    void discard_unused()
    {
        for (auto it = series_by_key_.begin(); it != series_by_key_.end();)
        {
            if (!used_.count(it->first))
            {
                first_valid_.erase(it->first);
                it = series_by_key_.erase(it);
            }
            else
                ++it;
        }
    }

    // Targeted eviction preserves references to other entries. Normal searches
    // use discard_unused() after each trial.
    void erase(const std::string &k)
    {
        series_by_key_.erase(k);
        first_valid_.erase(k);
    }

    void clear()
    {
        series_by_key_.clear();
        first_valid_.clear();
    }

    size_t size() const
    {
        return series_by_key_.size();
    }

private:
    std::unordered_map<std::string, std::vector<float>> series_by_key_;
    std::unordered_map<std::string, size_t> first_valid_;
    mutable std::unordered_set<std::string> used_;
};
