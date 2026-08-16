#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "tools_fatal.hh"

// Per-pair store of precomputed indicator series, keyed by indicator name plus its
// parameters.
//
// This replaces two things in KLINEf:
//
//   * std::array<std::vector<float>, 1000> EMA (and EMA_1h), indexed directly by
//     period. Any period >= 1000 was silent out-of-bounds, and the two arrays cost
//     ~48 KB of mostly-empty vector headers per pair.
//   * a fixed list of named members (ATR, StochRSI, StochRSI_K/D, AO, WILLR,
//     BollB_U/M/L, SuperTrend_1h) that every strategy shared whether it used them or
//     not -- so adding an indicator to one strategy meant editing the struct that all
//     of them depend on.
//
// PERFORMANCE: get() hashes a string. Hoist the reference out of the per-bar loop:
//
//     const std::vector<float> &ema = pair.indicators.get(IndicatorCache::key("EMA", n));
//     for (uint ii = ...) { ... ema[ii] ... }
//
// not `pair.indicators.get(...)[ii]` inside the loop, which is the mistake the
// string-keyed EMA_LISTS maps used to make once or twice per bar per pair.
class IndicatorCache
{
public:
    // "EMA:200", "STOCHRSI:14:14", "SUPERTREND:10:3.000000" -- readable in a debugger
    // and collision-free as long as the name is unique per indicator.
    static std::string key(const std::string &name)
    {
        return name;
    }

    template <typename T, typename... Rest>
    static std::string key(const std::string &name, const T &first, const Rest &...rest)
    {
        return key(name + ":" + std::to_string(first), rest...);
    }

    void put(const std::string &k, std::vector<float> series)
    {
        series_by_key_[k] = std::move(series);
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
            BACKTEST_FATAL("IndicatorCache: no series for key '" + k + "'. Compute it in CALCULATE_INDICATORS first.");
        }
        return it->second;
    }

    // Drop one series. Used by strategies that cache a sweep-expensive indicator
    // lazily: they evict the previous parameter's series when the sweep moves on, so
    // the cache stays a fixed size instead of accumulating every combination.
    //
    // Safe to call while holding references to *other* entries: unordered_map keeps
    // references to surviving elements valid across erase and rehash.
    void erase(const std::string &k)
    {
        series_by_key_.erase(k);
    }

    void clear()
    {
        series_by_key_.clear();
    }

    size_t size() const
    {
        return series_by_key_.size();
    }

private:
    std::unordered_map<std::string, std::vector<float>> series_by_key_;
};
