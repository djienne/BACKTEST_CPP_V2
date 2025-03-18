#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <math.h>
#include <unordered_map>
#include "tools.hh"
#include "custom_talib_wrapper.hh"
#include <ta-lib/ta_libc.h>
using namespace std;
using uint = unsigned int;

const string STRAT_NAME = "2EMA_crossover_StochRSI";
const string out_filename = STRAT_NAME + "_best.txt";

static const uint NB_PAIRS = 2;
const uint NB_POSITION_MAX = 2;
const vector<string> COINS = {"BTC", "ETH"};
const string timeframe = "1h";

const vector<string> DATAFILES = {"./data/data/binance/" + timeframe + "/" + COINS[0] + "-USDT.csv",
                                  "./data/data/binance/" + timeframe + "/" + COINS[1] + "-USDT.csv"};

const float start_year = 2017; // forced year to start (applies if data below is available)
const float FEE = 0.1f;        // FEES in %
const float USDT_amount_initial = 1000.0f;
const uint MIN_NUMBER_OF_TRADES = 100;         // minimum number of trades required (to avoid some noise / lucky circunstances)
const float MIN_ALLOWED_MAX_DRAWBACK = -50.0f; // %
const float STOCH_RSI_UPPER = 0.800;
const float STOCH_RSI_LOWER = 0.200;

std::vector<uint> start_indexes{};

// RANGE OF EMA PERIDOS TO TESTs
const int period_max_EMA = 600;
// const int range_step = 2;
// vector<int> range_EMA = {180};
vector<int> range_EMA_short = integer_range(2, 300 + 4, 1);
vector<int> range_EMA_long = integer_range(70, period_max_EMA + 4, 1);
//////////////////////////
array<std::unordered_map<string, vector<float>>, NB_PAIRS> EMA_LISTS{};
array<vector<float>, NB_PAIRS> StochRSI{};

array<long int, NB_PAIRS> last_times{};

uint i_print = 0;
uint nb_tested = 0;

RUN_RESULTf best{};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_best_res(const RUN_RESULTf &bestt)
{
    std::cout << "\n--------------------------------------------------------------------------" << std::endl;
    std::cout << "BEST PARAMETER SET FOUND: " << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "Time             : " << GREY << GET_CURRENT_TIME_STR() << RESET << std::endl;
    std::cout << "Strategy         : " << BLUE << STRAT_NAME << RESET << std::endl;
    std::cout << "Parameters       : " << YELLOW << bestt.param_str << RESET << std::endl;
    std::cout << "Gain             : " << bestt.gain_pc << "%" << std::endl;
    std::cout << "Porfolio         : " << bestt.WALLET_VAL_USDT << "$ (started with 1000$)" << std::endl;
    std::cout << "Win rate         : " << bestt.win_rate << "%" << std::endl;
    std::cout << "max DD           : " << bestt.max_DD << "%" << std::endl;
    std::cout << "Gain/DDC         : " << bestt.gain_over_DDC << std::endl;
    std::cout << "Score            : " << GREEN << bestt.score << RESET << std::endl;
    std::cout << "Calmar ratio     : " << bestt.calmar_ratio << std::endl;
    std::cout << "Number of trades : " << bestt.nb_posi_entered << std::endl;
    std::cout << "Total fees paid  : " << round(bestt.total_fees_paid * 100.0f) / 100.0f << "$ (started with 1000$)" << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(const vector<KLINEf> &PAIRS, const int ema_s, const int ema_l)
{
    nb_tested++;

    RUN_RESULTf result{};

    vector<float> USDT_tracking{};
    vector<int> USDT_tracking_ts{};

    const uint nb_max = PAIRS[0].nb;

    bool LAST_ITERATION = false;
    bool OPEN_LONG_CONDI = false;
    bool CLOSE_LONG_CONDI = false;
    uint nb_profit = 0;
    uint nb_loss = 0;
    uint NB_POSI_ENTERED = 0;
    float pc_change_with_max = 0, max_drawdown = 0;
    float USDT_amount = USDT_amount_initial;
    float MAX_WALLET_VAL_USDT = USDT_amount_initial;
    float total_fees_paid_USDT = 0.0f;
    float WALLET_VAL_USDT = USDT_amount_initial;
    array<float, NB_PAIRS> price_position_open{};
    array<float, NB_PAIRS> COIN_AMOUNTS{};
    array<float, NB_PAIRS> take_profit{};
    array<float, NB_PAIRS> stop_loss{};
    array<float, NB_PAIRS> stop_loss_at_open{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        COIN_AMOUNTS[ic] = 0.0f;
    }
    uint ACTIVE_POSITIONS = 0;

    const uint ii_begin = start_indexes[0];

    for (uint ii = ii_begin; ii < nb_max; ii++)
    {
        if (ii == nb_max - 1)
            LAST_ITERATION = true;

        bool closed = false;
        // For all pairs, check to close / open positions
        for (uint ic = 0; ic < NB_PAIRS; ic++)
        {
            if (ii < start_indexes[ic])
                continue;

            // conditions for open / close position

            OPEN_LONG_CONDI = EMA_LISTS[ic]["EMA_" + std::to_string(ema_s)][ii] >= EMA_LISTS[ic]["EMA_" + std::to_string(ema_l)][ii] && StochRSI[ic][ii] < STOCH_RSI_UPPER;
            CLOSE_LONG_CONDI = EMA_LISTS[ic]["EMA_" + std::to_string(ema_s)][ii] <= EMA_LISTS[ic]["EMA_" + std::to_string(ema_l)][ii] && StochRSI[ic][ii] > STOCH_RSI_LOWER;

            // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING POSITION AND ONLY THEN FOR OPENING POSITION

            // CLOSE LONG
            if (COIN_AMOUNTS[ic] > 0.0f && (CLOSE_LONG_CONDI || LAST_ITERATION))
            {
                const float to_add = COIN_AMOUNTS[ic] * PAIRS[ic].close[ii];
                USDT_amount += to_add;
                COIN_AMOUNTS[ic] = 0.0f;

                ACTIVE_POSITIONS--;

                // apply FEEs
                const float fe = to_add * FEE / 100.0f;
                USDT_amount -= fe;
                total_fees_paid_USDT += fe;
                //
                if (PAIRS[ic].close[ii] >= price_position_open[ic])
                {
                    nb_profit++;
                }
                else
                {
                    nb_loss++;
                }
                closed = true;
            }

            // OPEN LONG
            if (COIN_AMOUNTS[ic] == 0.0f && OPEN_LONG_CONDI && LAST_ITERATION == false && ACTIVE_POSITIONS < NB_POSITION_MAX)
            {
                price_position_open[ic] = PAIRS[ic].close[ii];

                const float usdMultiplier = 1.0f / float(NB_POSITION_MAX - ACTIVE_POSITIONS);

                COIN_AMOUNTS[ic] = USDT_amount * usdMultiplier / PAIRS[ic].close[ii];
                USDT_amount -= USDT_amount * usdMultiplier;

                // apply FEEs
                const float fe = COIN_AMOUNTS[ic] * FEE / 100.0f;
                COIN_AMOUNTS[ic] -= fe;
                total_fees_paid_USDT += fe * PAIRS[ic].close[ii];
                //

                ACTIVE_POSITIONS++;
                NB_POSI_ENTERED++;
            }
        }

        // check wallet status
        if (closed || LAST_ITERATION)
        {
            array<float, NB_PAIRS> closes{};
            for (uint ic = 0; ic < NB_PAIRS; ic++)
                closes[ic] = PAIRS[ic].close[ii];

            WALLET_VAL_USDT = USDT_amount + vector_product<NB_PAIRS>(COIN_AMOUNTS, closes);
            if (WALLET_VAL_USDT > MAX_WALLET_VAL_USDT)
                MAX_WALLET_VAL_USDT = WALLET_VAL_USDT;

            pc_change_with_max = (WALLET_VAL_USDT - MAX_WALLET_VAL_USDT) / MAX_WALLET_VAL_USDT * 100.0f;
            if (pc_change_with_max < max_drawdown)
                max_drawdown = pc_change_with_max;

            USDT_tracking.push_back(WALLET_VAL_USDT);
            USDT_tracking_ts.push_back(PAIRS[0].timestamp[ii]);
        }
    }

    array<float, NB_PAIRS> last_closes{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
        last_closes[ic] = PAIRS[ic].close[nb_max - 1];

    WALLET_VAL_USDT = USDT_amount + vector_product<NB_PAIRS>(COIN_AMOUNTS, last_closes);

    const float gain = (WALLET_VAL_USDT - USDT_amount_initial) / USDT_amount_initial * 100.0f;
    const float WR = float(nb_profit) / float(NB_POSI_ENTERED) * 100.0f;
    const float DDC = (1.0f / (1.0f + max_drawdown / 100.0f) - 1.0f) * 100.0f;
    const float score = gain / DDC * WR;

    i_print++;

    result.param_str = "\n  EMA_short: " + std::to_string(ema_s) + " ; EMA_long: " + std::to_string(ema_l) +
                       "\n  STOCH_RSI_LOWER: " + std::to_string(STOCH_RSI_LOWER) + " ; STOCH_RSI_UPPER: " + std::to_string(STOCH_RSI_UPPER);
    result.WALLET_VAL_USDT = USDT_amount;
    result.gain_over_DDC = gain / DDC;
    result.gain_pc = gain;
    result.max_DD = max_drawdown;
    result.nb_posi_entered = NB_POSI_ENTERED;
    result.win_rate = WR;
    result.score = score;
    result.calmar_ratio = calculate_calmar_ratio(USDT_tracking_ts, USDT_tracking, DDC);
    result.ema1 = ema_s;
    result.ema2 = ema_l;
    result.total_fees_paid = total_fees_paid_USDT;
    result.max_open_trades = NB_POSITION_MAX;

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CALCULATE_INDICATORS(std::vector<KLINEf> &PAIRS)
{
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        std::cout << "Calculating for " << COINS[ic] << std::endl;

        // std::cout << "Calculated STOCHRSI." << std::endl;

        std::vector<int> ema_values = combineAndRemoveDuplicates(range_EMA_long, range_EMA_short);

        for (const int ema_per : ema_values)
        {
            EMA_LISTS[ic]["EMA_" + std::to_string(ema_per)] = TALIB_EMA(PAIRS[ic].close, ema_per);
        }
        // std::cout << "Calculated EMAs." << std::endl;
    }
    cout << "Calculated EMAs." << std::endl;

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        StochRSI[ic] = TALIB_STOCHRSI_not_averaged(PAIRS[ic].close, 14, 14);
    }
    cout << "Calculated STOCHRSI." << std::endl;

    std::cout << "Initialized calculations." << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    const double t_begin = get_wall_time();
    std::cout << "\n-------------------------------------------------------------------------" << std::endl;
    std::cout << "Strategy to test : " << BLUE << STRAT_NAME << RESET << std::endl;
    std::cout << "DATA FILES TO PROCESS: " << std::endl;
    for (const string &dataf : DATAFILES)
    {
        std::cout << " " << YELLOW << "  " << dataf << RESET << std::endl;
    }

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

    vector<KLINEf> PAIRS;
    PAIRS.reserve(NB_PAIRS);
    for (const string &dataf : DATAFILES)
    {
        PAIRS.push_back(read_input_data(dataf));
    }

    random_shuffle_vector(range_EMA_long);

    start_indexes = INITIALIZE_DATA(PAIRS); // this function modifies PAIRS
    CALCULATE_INDICATORS(PAIRS);

    best.gain_over_DDC = -100.0f;
    best.calmar_ratio = -100.0f;
    best.score = -100.0f;

    const uint last_idx = PAIRS[0].nb - 1;

    const int year = get_year_from_timestamp(PAIRS[0].timestamp[0]);
    const int month = get_month_from_timestamp(PAIRS[0].timestamp[0]);
    const int day = get_day_from_timestamp(PAIRS[0].timestamp[0]);

    const int last_year = get_year_from_timestamp(PAIRS[0].timestamp[last_idx]);
    const int last_month = get_month_from_timestamp(PAIRS[0].timestamp[last_idx]);
    const int last_day = get_day_from_timestamp(PAIRS[0].timestamp[last_idx]);

    std::time_t difference = std::abs(int(PAIRS[0].timestamp[last_idx]) - int(PAIRS[0].timestamp[0]));
    const int days = difference / (24 * 60 * 60);

    // Display info
    std::cout << "Begin day                         : " << year << "/" << month << "/" << day << std::endl;
    std::cout << "End day                           : " << last_year << "/" << last_month << "/" << last_day << std::endl;
    std::cout << "Number of days                    : " << YELLOW << days << RESET << std::endl;
    std::cout << "OPEN/CLOSE FEE                    : " << FEE << " %" << std::endl;
    std::cout << "Maximum number open trades        : " << NB_POSITION_MAX << std::endl;
    std::cout << "Minimum number of trades required : " << MIN_NUMBER_OF_TRADES << std::endl;
    std::cout << "Maximum drawdown allowed          : " << MIN_ALLOWED_MAX_DRAWBACK << " %" << std::endl;
    std::cout << "StochRSI Upper Band               : " << STOCH_RSI_UPPER << std::endl;
    std::cout << "StochRSI Lower Band               : " << STOCH_RSI_LOWER << std::endl;
    std::cout << "EMA short period max tested       : " << find_max(range_EMA_short) << std::endl;
    std::cout << "EMA long period max tested        : " << find_max(range_EMA_long) << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;

    std::vector<doubleEMA_params> param_list{};
    param_list.reserve(range_EMA_short.size() * range_EMA_long.size());

    for (const int ema1 : range_EMA_long)
    {
        for (const int ema2 : range_EMA_short)
        {
            if (std::abs(ema1 - ema2) < 5)
                continue;

            const doubleEMA_params to_add{ema1, ema2};

            param_list.push_back(to_add);
        }
    }

    random_shuffle_vector(param_list);

    // MAIN LOOP

    uint nb_done = 0;
    uint i_print3 = 0;

    for (const auto para : param_list)
    {
        const RUN_RESULTf res = PROCESS(PAIRS, para.ema1, para.ema2);
        nb_done++;
        i_print3++;

        if (res.score > best.score && res.gain_pc < 1000000.0f && res.nb_posi_entered >= MIN_NUMBER_OF_TRADES && res.max_DD > MIN_ALLOWED_MAX_DRAWBACK)
        {
            best = res;
        }

        if (i_print3 == 1000)
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

    const double t_end = get_wall_time();

    std::cout << "Number of backtests performed : " << nb_tested << std::endl;
    std::cout << "Time taken                    : " << t_end - t_begin << " seconds " << std::endl;
    const double ram_usage = process_mem_usage();
    std::cout << "RAM usage                     : " << std::round(ram_usage * 10.0) / 10.0 << " MB" << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;

    TA_Shutdown();

    return 0;
}
