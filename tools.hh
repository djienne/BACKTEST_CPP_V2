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
#include <algorithm> // std::shuffle
#include <random>    // std::default_random_engine
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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
constexpr const char *RESET = "\033[0m";
constexpr const char *RED = "\033[31m";
constexpr const char *GREEN = "\033[32m";
constexpr const char *YELLOW = "\033[33m";
constexpr const char *BLUE = "\033[34m";
constexpr const char *MAGENTA = "\033[35m";
constexpr const char *CYAN = "\033[36m";
constexpr const char *WHITE = "\033[37m";
constexpr const char *GREY = "\033[90m";
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Strategy-agnostic run result. New strategies should stash their parameters in
// `param_str` instead of adding numeric fields here. `ema1`/`ema2` remain because the
// 2-EMA strategy reads them to display a compact legacy banner.
//
// Monetary and score fields are double. Price and indicator arrays stay float -- that
// is where the memory win lives -- but the wallet is an accumulator: thousands of
// fee subtractions and position round-trips ran through float's ~7 significant digits,
// and every ranking decision the sweep makes depends on those totals.
struct RUN_RESULTf
{
    double WALLET_VAL_USDT = 0.0;
    double gain_pc = 0.0;
    double win_rate = 0.0;
    double max_DD = 0.0;
    double gain_over_DDC = 0.0;
    double score = 0.0;
    int nb_posi_entered = 0;
    int ema1 = 0;
    int ema2 = 0;
    double total_fees_paid = 0.0;
    double calmar_ratio = 0.0;
    double calmar_ratio_monthly = 0.0;
    uint max_open_trades = 0;
    std::string param_str;
};

struct doubleEMA_params
{
    int ema1;
    int ema2;
};

struct trix_params
{
    int ema1;
    int trixLength;
    int trixSignal;
    uint max_open_trades;
};

struct BigWill_params
{
    int AO_fast;
    int AO_slow;
    int ema_f;
    int ema_s;
    uint max_open_trades;
};

struct EMA3_params
{
    int ema1;
    int ema2;
    int ema3;
    float up;
    float down;
    float SRSIL;
    float SRSIU;
    uint max_open_trades;
};

struct SR_params
{
    int ema_fast;
    int ema_slow;
    uint max_open_trades;
};

struct ST_EMA_ATR_params
{
    int ema;
    float up;
    float down;
    uint max_open_trades;
};

struct BBTREND_params
{
    int ema;
    int BBlength;
    float BBstd;
    uint max_open_trades;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Position sizes accumulate in double while prices stay float, so the dot product is
// accumulated in double.
template <size_t N>
double vector_product(const std::array<double, N> &vec, const std::array<float, N> &vec2)
{
    double out = 0.0;
    for (size_t i = 0; i < vec.size(); i++)
    {
        out += vec[i] * static_cast<double>(vec2[i]);
    }
    return out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float find_average(const std::vector<float> &vec);
double find_average(const std::vector<double> &vec);
float find_min(const std::vector<float> &vec);
float find_max(const std::vector<float> &vec);
int find_max(const std::vector<int> &vec);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Calendar fields, always interpreted in UTC (see tools.cpp for why).
int get_hour_from_timestamp(const int64_t timestamp);

int get_year_from_timestamp(const int64_t timestamp);

int get_month_from_timestamp(const int64_t timestamp);

int get_day_from_timestamp(const int64_t timestamp);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double get_wall_time();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double process_mem_usage();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<int> integer_range(const int min, const int max, const int step);
std::vector<int> integer_range(const int min, const int max);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> float_Nvalues_range(const float &vmin, const float &vmax, const int &N);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string ReplaceAll(std::string str, const std::string &from, const std::string &to);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// `drawdown_normalizer` is the value the averaged periodic return is divided by. Every
// strategy passes ResultMetrics::ddc (the drawdown-corrected gain, a positive number),
// not the raw negative max drawdown -- the parameter used to be named `max_DD`, which
// described neither what callers pass nor what the function does with it.
double calculate_calmar_ratio(const std::vector<int64_t> &times, const std::vector<double> &wallet_vals, const double drawdown_normalizer);
double calculate_calmar_ratio_monthly(const std::vector<int64_t> &times, const std::vector<double> &wallet_vals, const double drawdown_normalizer);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// `verbose` logs the seed, which is useful once per run but not when shuffling many
// small groups in a loop -- pass false there.
template <typename T>
void random_shuffle_vector(std::vector<T> &vec_in, const bool verbose = true)
{
    std::random_device rd;
    std::mt19937::result_type seed = rd() + std::hash<std::thread::id>{}(std::this_thread::get_id());
    std::mt19937 rng(seed);
    std::shuffle(vec_in.begin(), vec_in.end(), rng);

    if (verbose)
    {
        std::cout << "Thread " << std::this_thread::get_id() << " Random number seed: " << seed << std::endl;
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void WRITE_OR_UPDATE_BEST_SCORE_FILE(const std::string &STRAT_NAME, const std::string &out_filename, const RUN_RESULTf &result);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string GET_CURRENT_TIME_STR();

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool check_timestamp_consistencies(const std::vector<KLINEf> &PAIRS);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void realign_timestamps(const KLINEf &klines_btc, KLINEf &klines2);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
std::vector<T> add_zeros(const std::vector<T> &vec_in, const int &nb_to_add)
{
    std::vector<T> vec_to_add(nb_to_add, 0);
    vec_to_add.insert(vec_to_add.end(), vec_in.begin(), vec_in.end());
    return vec_to_add;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

KLINEf read_input_data(const std::string &input_file_path);
KLINEf read_input_data_f(const std::string &input_file_path, const std::string &max_time = "2099-06-20");
fundings read_funding_rates_data(const std::string &input_file_path);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<uint> INITIALIZE_DATA(std::vector<KLINEf> &PAIRS);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template <typename T>
std::vector<T> combineAndRemoveDuplicates(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
    std::set<T> uniqueElements;
    std::vector<T> combinedVector;

    // Insert elements from vec1
    uniqueElements.insert(vec1.begin(), vec1.end());

    // Insert elements from vec2
    uniqueElements.insert(vec2.begin(), vec2.end());

    // Copy unique elements to combinedVector
    combinedVector.insert(combinedVector.end(), uniqueElements.begin(), uniqueElements.end());

    return combinedVector;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class RandomNumberGenerator
{
public:
    // The thread id is mixed into the seed so that two workers constructing a generator
    // in the same tick still diverge. Previously the constructor body built a *local*
    // mt19937 that shadowed the member and was then discarded, so this mixing never
    // took effect and the seed printed below was not the one actually in use.
    RandomNumberGenerator()
    {
        std::random_device rd;
        const std::mt19937::result_type seed = rd() + std::hash<std::thread::id>{}(std::this_thread::get_id());
        rng.seed(seed);

        std::cout << "Thread " << std::this_thread::get_id() << ": Random number seed = " << seed << std::endl;
    }

    int getRandomNumber(const int upperLimit)
    {
        std::uniform_int_distribution<int> distribution(0, upperLimit);
        return distribution(rng);
    }

private:
    std::mt19937 rng{};
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float get_funding_fee_if_any(const fundings &FUND, const int64_t current_timestamp);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mark-to-market value of a futures book. A short is valued as its mirrored long,
// (2 * entry - price), which is exact at leverage 1 and has no liquidation model --
// see the note on the short helpers in trade_core.hh.
template <size_t N>
double calculate_wallet_val_usdt(const double USDT_amount, const std::array<double, N> &COIN_AMOUNTS, const std::array<float, N> &current_prices, const std::array<double, N> &prices_position_open)
{
    double VAL = USDT_amount;

    for (size_t ic = 0; ic < COIN_AMOUNTS.size(); ic++)
    {
        if (COIN_AMOUNTS[ic] > 0.0)
        {
            VAL += COIN_AMOUNTS[ic] * static_cast<double>(current_prices[ic]);
        }
        else if (COIN_AMOUNTS[ic] < 0.0)
        {
            VAL += std::abs(COIN_AMOUNTS[ic]) * (2.0 * prices_position_open[ic] - static_cast<double>(current_prices[ic]));
        }
    }

    return VAL;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fill_datafile_paths_f(const std::vector<std::string> &COINS, const std::string &timeframe, std::vector<std::string> &DATAFILES, std::vector<std::string> &DATAFILES_fundings);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<int> generateRange_int(const int &vmin, const int &vmax, const int &N);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

time_t convertToUnixTimestamp(const std::string &dateString);
std::string getCurrentDateMinusTwoDays();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
