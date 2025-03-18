#pragma once
#include <stdio.h>
#include <thread>
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
using json = nlohmann::json;
using namespace std::chrono;
using namespace std;
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

struct RUN_RESULTf
{
    float WALLET_VAL_USDT;
    float gain_pc;
    float win_rate;
    float max_DD;
    float gain_over_DDC;
    float score;
    int nb_posi_entered;
    int ema1;
    int ema2;
    int ema3;
    int ema4;
    int AO_fast;
    int AO_slow;
    int trixLength;
    int trixSignal;
    float UP;
    float DOWN;
    float min_yearly_gain;
    float max_yearly_gain;
    std::vector<float> yearly_gains;
    std::vector<float> years_yearly_gains;
    float RSI_limit;
    float RSI_limit2;
    float gain_limit;
    float up;
    float down;
    float SRSIL;
    float total_fees_paid;
    int max_delta_t_new_ATH;
    float calmar_ratio;
    float calmar_ratio_monthly;
    uint max_open_trades;
    float AMOUNT_USDT;
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

template <int N>
float vector_product(const std::array<float, N> &vec, const std::array<float, N> &vec2)
{
    float out = 0.0;
    for (uint i = 0; i < vec.size(); i++)
    {
        out += vec[i] * vec2[i];
    }
    return out;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float find_average(const std::vector<float> &vec);
float find_min(const std::vector<float> &vec);
float find_max(const std::vector<float> &vec);
int find_max(const std::vector<int> &vec);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int get_hour_from_timestamp(const int timestamp);

int get_year_from_timestamp(const int timestamp);

int get_month_from_timestamp(const int &timestamp);

int get_day_from_timestamp(const int timestamp);

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

float calculate_calmar_ratio(const std::vector<int> &times, const std::vector<float> &wallet_vals, const float &max_DD);
float calculate_calmar_ratio_monthly(const std::vector<int> &times, const std::vector<float> &wallet_vals, const float &max_DD);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
void random_shuffle_vector(std::vector<T> &vec_in)
{
    std::random_device rd;
    std::mt19937::result_type seed = rd() + std::hash<std::thread::id>{}(std::this_thread::get_id());
    std::mt19937 rng(seed);
    std::shuffle(vec_in.begin(), vec_in.end(), rng);

    std::cout << "Thread " << std::this_thread::get_id() << " Random number seed: " << seed << std::endl;
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

bool key_exists(std::unordered_map<std::string, std::vector<float>> m, const std::string &ch);

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

std::vector<std::vector<float>> splitVector(const std::vector<float> &inputVector, const int splitSize);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> duplicateElements(const std::vector<float> &inputVector, const int n);

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
    RandomNumberGenerator() : rng(std::random_device()())
    {
        std::random_device rd;
        std::mt19937::result_type seed = rd() + std::hash<std::thread::id>{}(std::this_thread::get_id());
        std::mt19937 rng(seed);

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

float get_funding_fee_if_any(const fundings &FUND, const int &current_timestamp);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <int N>
float calculate_wallet_val_usdt(const float &USDT_amount, const std::array<float, N> &COIN_AMOUNTS, const std::array<float, N> &current_prices, const std::array<float, N> &prices_position_open)
{
    float VAL = USDT_amount;

    for (int ic = 0; ic < COIN_AMOUNTS.size(); ic++)
    {
        if (COIN_AMOUNTS[ic] > 0.0f)
        {
            VAL += COIN_AMOUNTS[ic] * current_prices[ic];
        }
        else if (COIN_AMOUNTS[ic] < 0.0f)
        {
            VAL += std::abs(COIN_AMOUNTS[ic]) * (2.0f * prices_position_open[ic] - current_prices[ic]);
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
