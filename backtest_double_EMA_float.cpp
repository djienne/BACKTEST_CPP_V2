#include <iostream>
#include <ctime>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <unordered_map>
#include "tools.hh"
#include "trade_core.hh"
#include "custom_talib_wrapper.hh"
#include <ta-lib/ta_libc.h>

using namespace std;

const string STRAT_NAME = "2EMA_crossover";
const string out_filename = STRAT_NAME + "_best.txt";

const float FEE = 0.1f;                   // FEES in %
const float USDT_amount_initial = 1000.0f;
const double MAX_ALLOWED_DD = -40.0;
const std::string PAIR = "BTC-USDT";
const std::string DATAFILE = "./data/data/binance/1h/" + PAIR + ".csv";

// RANGE OF EMA PERIDOS TO TEST
const int MIN_EMA = 3;
std::vector<int> range1 = integer_range(MIN_EMA, 600);
std::vector<int> range2 = integer_range(MIN_EMA, 600);
//////////////////////////

uint i_print = 0;

std::unordered_map<std::string, std::vector<float>> EMA_LISTS{};
std::vector<float> year{};
std::vector<float> hour{};
std::vector<float> month{};
std::vector<float> day{};

KLINEf kline{};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_best_res(const RUN_RESULTf best)
{
    const int year = get_year_from_timestamp(kline.timestamp[0]);
    const int month = get_month_from_timestamp(kline.timestamp[0]);
    const int day = get_day_from_timestamp(kline.timestamp[0]);
    const uint last_idx = kline.nb - 1;
    const int last_year = get_year_from_timestamp(kline.timestamp[last_idx]);
    const int last_month = get_month_from_timestamp(kline.timestamp[last_idx]);
    const int last_day = get_day_from_timestamp(kline.timestamp[last_idx]);

    std::time_t difference = std::abs(int(kline.timestamp[last_idx]) - int(kline.timestamp[0]));
    const int days = difference / (24 * 60 * 60);

    // Display info
    std::cout << "\n-------------------------------------" << std::endl;
    std::cout << "TIME RANGE: " << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << " Begin day : " << year << "/" << month << "/" << day << endl;
    std::cout << " End day   : " << last_year << "/" << last_month << "/" << last_day << endl;
    std::cout << " Duration  : " << days << " days" << endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "BEST PARAMETER SET FOUND: " << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << " Time             : " << GREY << GET_CURRENT_TIME_STR() << RESET << std::endl;
    std::cout << " Strategy         : " << BLUE << STRAT_NAME << RESET << std::endl;
    std::cout << " PAIR             : " << PAIR << std::endl;
    std::cout << " EMAs             : " << best.ema1 << " " << best.ema2 << std::endl;
    std::cout << " Gain             : " << best.gain_pc << "%" << std::endl;
    std::cout << " Win rate         : " << best.win_rate << "%" << std::endl;
    std::cout << " max DD           : " << best.max_DD << "%" << std::endl;
    std::cout << " Gain over DDC    : " << best.gain_over_DDC << std::endl;
    std::cout << " Score            : " << GREEN << best.score << RESET << std::endl;
    std::cout << " Number of trades : " << best.nb_posi_entered << std::endl;
    std::cout << "-------------------------------------"
              << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void INITIALIZE_DATA0(const KLINEf &kline)
{
    std::cout << "Initializing data (calculating indicators)...\n"
              << std::endl;

    std::vector<int> list_ema = {};
    EMA_LISTS.clear();
    year.clear();
    hour.clear();
    month.clear();
    day.clear();
    year.reserve(kline.nb);
    hour.reserve(kline.nb);
    month.reserve(kline.nb);
    day.reserve(kline.nb);

    for (int i = MIN_EMA; i <= range2.size() + 5; i++)
    {
        list_ema.push_back(i);
    }

    EMA_LISTS.reserve(list_ema.size());

    for (const int i : list_ema)
    {
        EMA_LISTS["EMA" + std::to_string(i)] = TALIB_EMA(kline.close, i);
    }
    cout << "Done calculating EMAs." << endl;

    for (uint ii = 0; ii < kline.nb; ii++)
    {
        year.push_back(get_year_from_timestamp(kline.timestamp[ii]));
        hour.push_back(get_hour_from_timestamp(kline.timestamp[ii]));
        month.push_back(get_month_from_timestamp(kline.timestamp[ii]));
        day.push_back(get_day_from_timestamp(kline.timestamp[ii]));
    }

    cout << "Initialized calculations." << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(const KLINEf &kline_data, const int ema1_v, const int ema2_v)
{
    const std::vector<float> &EMA1 = EMA_LISTS["EMA" + std::to_string(ema1_v)];
    const std::vector<float> &EMA2 = EMA_LISTS["EMA" + std::to_string(ema2_v)];
    const std::vector<float> &close = kline_data.close;
    const int nb_max = kline_data.nb;
    const int ii_begin = find_max(range1) + 2;

    trade_core::TradeStats stats{};
    trade_core::PortfolioState<1> portfolio(USDT_amount_initial);

    // Signal-only strategy on a single asset: use trade_core's spot long helpers.
    // Drawdown is tracked off realized USDT (updated only at position transitions),
    // which matches the previous inline-PnL behavior exactly since USDT is constant
    // while a position is open.
    for (int ii = ii_begin; ii < nb_max; ii++)
    {
        const bool LAST_ITERATION = (ii == nb_max - 1);
        const bool OPEN_LONG_CONDI = (EMA2[ii] >= EMA1[ii]) && (EMA2[ii - 1] <= EMA1[ii - 1]);
        const bool CLOSE_LONG_CONDI = (EMA2[ii] <= EMA1[ii]) && (EMA2[ii - 1] >= EMA1[ii - 1]);

        // IMPORTANT: close before open.
        if (portfolio.coin_amounts[0] > 0.0f && (CLOSE_LONG_CONDI || LAST_ITERATION))
        {
            trade_core::close_spot_long(portfolio, stats, 0, close[ii], FEE);

            const float wallet_val_usdt = portfolio.usdt_amount;
            if (wallet_val_usdt > portfolio.max_wallet_val_usdt)
            {
                portfolio.max_wallet_val_usdt = wallet_val_usdt;
            }
            const float pc_change_with_max = (wallet_val_usdt - portfolio.max_wallet_val_usdt) / portfolio.max_wallet_val_usdt * 100.0f;
            if (pc_change_with_max < portfolio.max_drawdown)
            {
                portfolio.max_drawdown = pc_change_with_max;
            }
        }

        if (portfolio.coin_amounts[0] == 0.0f && OPEN_LONG_CONDI && !LAST_ITERATION)
        {
            trade_core::open_spot_long(portfolio, stats, 0, close[ii], FEE, 1);
        }
    }

    const float wallet_val_usdt = portfolio.usdt_amount;
    const trade_core::ResultMetrics metrics = trade_core::calculate_result_metrics(wallet_val_usdt, USDT_amount_initial, portfolio.max_drawdown, stats);

    i_print++;

    RUN_RESULTf result{};
    trade_core::populate_common_result(result, metrics, wallet_val_usdt, portfolio.max_drawdown, portfolio.total_fees_paid_usdt, stats, 1);
    result.ema1 = ema1_v;
    result.ema2 = ema2_v;
    result.param_str = "PAIR: " + PAIR + "; EMA1: " + std::to_string(ema1_v) + "; EMA2: " + std::to_string(ema2_v);
    return result;
}

int main()
{
    auto startTime = std::chrono::high_resolution_clock::now();

    std::cout << "\n-------------------------------------" << std::endl;
    std::cout << "Strategy to test: EMA crossover simple" << std::endl;
    std::cout << "DATA FILE TO PROCESS: " << DATAFILE << std::endl;

    TA_RetCode retCode;
    retCode = TA_Initialize();
    if (retCode != TA_SUCCESS)
    {
        std::cout << "Cannot initialize TA-Lib !\n"
                  << retCode << "\n";
    }
    else
    {
        std::cout << "Initialized TA-Lib !\n";
    }

    kline = read_input_data(DATAFILE);

    INITIALIZE_DATA0(kline);

    RUN_RESULTf best{};
    best.gain_over_DDC = -100.0;

    int i_print2 = 0;
    int i_print3 = 0;

    std::cout << "Begining day : " << year[0] << "/" << month[0] << "/" << day[0] << std::endl;
    std::cout << "End day      : " << year.back() << "/" << month.back() << "/" << day.back() << std::endl;
    std::cout << "CAN LONG     : 1" << std::endl;
    std::cout << "CAN SHORT    : 0" << std::endl;

    std::vector<doubleEMA_params> param_list{};
    param_list.reserve(range1.size() * range2.size());

    for (const int ema1 : range1)
    {
        for (const int ema2 : range2)
        {
            if (std::abs(ema1 - ema2) < 7 || ema1 < ema2)
                continue;

            const doubleEMA_params to_add{ema1, ema2};

            param_list.push_back(to_add);
        }
    }

    random_shuffle_vector(param_list);

    // MAIN LOOP

    uint nb_done = 0;

    for (const auto para : param_list)
    {
        RUN_RESULTf res = PROCESS(kline, para.ema1, para.ema2);
        nb_done++;

        if (res.score > best.score && res.gain_pc < 10000.0 && res.nb_posi_entered >= 100 && res.max_DD > MAX_ALLOWED_DD)
        {
            best = res;
        }

        i_print3++;
        if (i_print3 == 10000)
        {
            print_best_res(best);
            WRITE_OR_UPDATE_BEST_SCORE_FILE(STRAT_NAME, out_filename, best);
            i_print3 = 0;
            const double pc_done = std::round(double(nb_done) / double(param_list.size()) * 100.0 * 100.0) / 100.0;
            std::cout << "DONE " << nb_done << " / " << param_list.size() << "   = " << pc_done << "%" << std::endl;
        }
    }

    print_best_res(best);
    WRITE_OR_UPDATE_BEST_SCORE_FILE(STRAT_NAME, out_filename, best);

    // End the timer
    auto endTime = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration<double>(endTime - startTime).count();

    std::cout << "CPU wall time: " << duration << " seconds" << std::endl;
    std::cout << "-------------------------------------\n"
              << std::endl;
}
