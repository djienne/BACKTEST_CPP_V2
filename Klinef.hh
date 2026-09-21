#pragma once
#include <cstdint>
#include <vector>
#include <array>
#include <string>
#include <unordered_map>

#include "indicator_cache.hh"

using uint = unsigned int;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct KLINEf
{
    // Seconds since the Unix epoch, always UTC. int64_t rather than int: the calendar
    // helpers used to narrow to int, which overflows in January 2038.
    std::vector<int64_t> timestamp;
    std::vector<float> open;
    std::vector<float> high;
    std::vector<float> low;
    std::vector<float> close;
    // Base-asset volume, retained for volume-based indicators.
    std::vector<float> volume;
    uint nb = 0;
    std::string name;
    uint start_idx = 0;

    // Optional container cache. Searches keep their caches separately, per worker.
    IndicatorCache indicators;
};

struct fundings
{
    std::vector<int64_t> timestamp;
    std::vector<float> funding;
    // Exact lookup for legacy rate-only fixture consumers. Research uses FundingEvent.
    std::unordered_map<int64_t, float> funding_by_timestamp;
    uint nb = 0;
    std::string name;
};
