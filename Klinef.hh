#pragma once
#include <cstdint>
#include <vector>
#include <array>
#include <string>
#include <unordered_map>
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
    std::array<std::vector<float>, 1000> EMA;
    std::array<std::vector<float>, 1000> EMA_1h;
    uint nb;
    std::string name;
    uint start_idx;
    std::vector<float> ATR;
    std::vector<float> StochRSI_K;
    std::vector<float> StochRSI_D;
    std::vector<float> StochRSI;
    std::vector<float> AO;
    std::vector<float> WILLR;
    std::vector<float> BollB_U;
    std::vector<float> BollB_M;
    std::vector<float> BollB_L;
    std::vector<float> SuperTrend_1h;
};

struct fundings
{
    std::vector<int64_t> timestamp;
    std::vector<float> funding;
    // Mirrors `timestamp` for fast exact funding lookups during futures backtests.
    std::unordered_map<int64_t, float> funding_by_timestamp;
    uint nb;
    std::string name;
};
