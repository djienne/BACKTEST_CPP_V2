#pragma once
#include <stdio.h>
#include "tools_fatal.hh" // BACKTEST_FATAL
#include <vector>
#include <array>
#include <time.h>
#include <chrono>
#include <limits>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <algorithm>
#include <random>
#include <ctime>
#include <sstream>
#include <regex>
#include <filesystem>
#include <set>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <optional>
#include "nlohmann/json.hpp"
#include "Klinef.hh"
#include "data_io.hh"
constexpr const char *RESET = "\033[0m";
constexpr const char *RED = "\033[31m";
constexpr const char *GREEN = "\033[32m";
constexpr const char *YELLOW = "\033[33m";
constexpr const char *BLUE = "\033[34m";
constexpr const char *MAGENTA = "\033[35m";
constexpr const char *CYAN = "\033[36m";
constexpr const char *WHITE = "\033[37m";
constexpr const char *GREY = "\033[90m";

struct FillRecord
{
    int64_t timestamp;
    uint pair;
    std::string action;
    double price, quantity, commission, net_pnl;
};
struct RUN_RESULTf
{
    double WALLET_VAL_USDT = 0.0;
    double gain_pc = 0.0;
    double win_rate = 0.0;
    double max_DD = 0.0;
    double gain_over_DDC = 0.0;
    double score = 0.0;
    int nb_posi_entered = 0;
    double total_fees_paid = 0.0;
    std::optional<double> calmar_ratio;
    double net_funding = 0.0;
    uint ambiguous_bars = 0;
    bool valid = true;
    std::string invalid_reason;
    std::vector<FillRecord> fills;
    std::vector<double> equity;
    std::vector<int64_t> equity_times;
    uint max_open_trades = 0;
};

float find_min(const std::vector<float> &vec);
float find_max(const std::vector<float> &vec);

// Calendar fields, always interpreted in UTC (see tools.cpp for why).
int get_hour_from_timestamp(const int64_t timestamp);

int get_year_from_timestamp(const int64_t timestamp);

int get_month_from_timestamp(const int64_t timestamp);

int get_day_from_timestamp(const int64_t timestamp);

double get_wall_time();

double process_mem_usage();

std::vector<int> integer_range(const int min, const int max, const int step);
std::vector<int> integer_range(const int min, const int max);

std::vector<float> float_Nvalues_range(const float &vmin, const float &vmax, const int &N);

float get_funding_fee_if_any(const fundings &FUND, const int64_t current_timestamp);

// Mark-to-market value of a futures book. A short is valued as its mirrored long,
// abs(quantity) * (2 * entry - price), representing collateral plus short P&L.
// This 1x accounting model does not model exchange liquidation.
template <size_t N>
double calculate_wallet_val_usdt(const double USDT_amount, const std::array<double, N> &COIN_AMOUNTS,
                                 const std::array<float, N> &current_prices,
                                 const std::array<double, N> &prices_position_open, const size_t count = N)
{
    double VAL = USDT_amount;

    for (size_t ic = 0; ic < count && ic < COIN_AMOUNTS.size(); ic++)
    {
        if (COIN_AMOUNTS[ic] > 0.0)
        {
            VAL += COIN_AMOUNTS[ic] * static_cast<double>(current_prices[ic]);
        }
        else if (COIN_AMOUNTS[ic] < 0.0)
        {
            VAL +=
                std::abs(COIN_AMOUNTS[ic]) * (2.0 * prices_position_open[ic] - static_cast<double>(current_prices[ic]));
        }
    }

    return VAL;
}

std::vector<int> generateRange_int(const int &vmin, const int &vmax, const int &N);

time_t convertToUnixTimestamp(const std::string &dateString);
std::string getCurrentDateMinusTwoDays();
