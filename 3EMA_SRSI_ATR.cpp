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

const string STRAT_NAME = "EMA3_SRSI_ATR";
const string out_filename = STRAT_NAME + "_best.txt";

const std::vector<uint> MAX_OPEN_TRADES_TO_TEST{4, 5, 6, 7, 8, 9, 10, 11, 12};
const std::vector<string> COINS = {"BTC",
                                   "ETH",
                                   "BNB",
                                   "XRP",
                                   "TRX",
                                   "MATIC",
                                   "LTC",
                                   "XMR",
                                   "XLM",
                                   "EOS",
                                   "ETC"};

static const uint NB_PAIRS = 11;

string timeframe = "5m";

vector<string> DATAFILES{};

const float start_year = 2017; // forced year to start (applies if data below is available)
const float FEE = 0.1f;        // FEES in %
const float USDT_amount_initial = 1000.0f;

const uint MIN_NUMBER_OF_TRADES = 300; // minimum number of trades required (to avoid some noise / lucky circunstances)

const float MIN_ALLOWED_MAX_DRAWBACK = -80.0f; // %

std::vector<uint> start_indexes{};

// RANGE OF PARAMETERS TO TEST
const std::vector<int> range_EMA1 = integer_range(3, 100, 5);
const std::vector<int> range_EMA2 = integer_range(5, 400, 10);
const std::vector<int> range_EMA3 = integer_range(5, 580, 10);
const std::vector<float> range_UP{3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
const std::vector<float> range_DOWN{3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
const std::vector<float> range_STOCH_RSI_LOWER = {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f};
//////////////////////////

array<std::unordered_map<string, std::vector<float>>, NB_PAIRS> EMA_LISTS{};
std::array<std::vector<float>, NB_PAIRS> StochRSI_K{};
std::array<std::vector<float>, NB_PAIRS> StochRSI_D{};
std::array<std::vector<float>, NB_PAIRS> ATR{};

uint i_print = 0;
uint nb_tested = 0;

RUN_RESULTf best{};

uint end_timestamp_datasets = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fill_datafile_paths()
{
    for (uint i = 0; i < COINS.size(); i++)
    {
        DATAFILES.push_back("./data/data/binance/" + timeframe + "/" + COINS[i] + "-USDT.csv");
    }
}

void print_best_res(const RUN_RESULTf &bestt)
{
    std::cout << "\n-------------------------------------------------------------------------" << std::endl;
    std::cout << "BEST PARAMETER SET FOUND: " << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "Time                 : " << GREY << GET_CURRENT_TIME_STR() << RESET << std::endl;
    std::cout << "Strategy             : " << BLUE << STRAT_NAME << RESET << std::endl;
    std::cout << "Parameters           : " << YELLOW << bestt.param_str << RESET << std::endl;
    std::cout << "Max Open Trades      : " << YELLOW << bestt.max_open_trades << RESET << std::endl;
    std::cout << "Gain                 : " << bestt.gain_pc << "%" << std::endl;
    std::cout << "Porfolio             : " << bestt.WALLET_VAL_USDT << "$ (started with 1000$)" << std::endl;
    std::cout << "Win rate             : " << bestt.win_rate << "%" << std::endl;
    std::cout << "max DD               : " << bestt.max_DD << "%" << std::endl;
    std::cout << "Gain/DDC             : " << bestt.gain_over_DDC << std::endl;
    std::cout << "Score                : " << GREEN << bestt.score << RESET << std::endl;
    std::cout << "Calmar ratio monthly : " << bestt.calmar_ratio_monthly << std::endl;
    std::cout << "Calmar ratio         : " << bestt.calmar_ratio << std::endl;
    std::cout << "Number of trades     : " << bestt.nb_posi_entered << std::endl;
    std::cout << "Total fees paid      : " << round(bestt.total_fees_paid * 100.0f) / 100.0f << "$ (started with 1000$)" << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(const std::vector<KLINEf> &PAIRS, const int &ema1, const int &ema2, const int &ema3, const float &up, const float &down, const float &STOCH_RSI_LOWER, const uint &MAX_OPEN_TRADES)
{
    nb_tested++;

    RUN_RESULTf result{};

    std::vector<float> USDT_tracking{};
    std::vector<int> USDT_tracking_ts{};

    const uint nb_max = PAIRS[0].nb;

    bool LAST_ITERATION = false;
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
    array<float, NB_PAIRS> ATR_AT_OPEN{};
    array<uint, NB_PAIRS> OPEN_TS{};

    uint ACTIVE_POSITIONS = 0;

    const uint ii_begin = start_indexes[0];

    bool NEW_MONTH = false;

    for (uint ii = ii_begin + 1; ii < nb_max; ii++)
    {
        if (ii == nb_max - 1)
        {
            LAST_ITERATION = true;
        }

        const int month_b = get_month_from_timestamp(PAIRS[0].timestamp[ii - 1]);
        const int month = get_month_from_timestamp(PAIRS[0].timestamp[ii]);

        if (month_b != month)
        {
            NEW_MONTH = true;
        }
        else
        {
            NEW_MONTH = false;
        }

        bool closed = false;
        // For all pairs, check to close / open positions
        for (uint ic = 0; ic < NB_PAIRS; ic++)
        {
            if (ii < start_indexes[ic])
            {
                continue;
            }

            // conditions for open / close position

            const bool OPEN_LONG_CONDI = (EMA_LISTS[ic]["EMA_" + std::to_string(ema1)][ii] >= EMA_LISTS[ic]["EMA_" + std::to_string(ema2)][ii] && EMA_LISTS[ic]["EMA_" + std::to_string(ema2)][ii] >= EMA_LISTS[ic]["EMA_" + std::to_string(ema3)][ii] && PAIRS[ic].close[ii] >= EMA_LISTS[ic]["EMA_" + std::to_string(ema1)][ii] && StochRSI_K[ic][ii] < STOCH_RSI_LOWER && StochRSI_D[ic][ii] < STOCH_RSI_LOWER && StochRSI_K[ic][ii - 1] > StochRSI_D[ic][ii - 1] && StochRSI_K[ic][ii] <= StochRSI_D[ic][ii]);

            bool timeout = false;
            bool hard_TP_condition = false;
            if (COIN_AMOUNTS[ic] != 0.0f)
            {
                timeout = (PAIRS[ic].timestamp[ii] - OPEN_TS[ic]) >= 2 * 24 * 3600;
                const float pc_gain = (PAIRS[ic].close[ii] - price_position_open[ic]) / price_position_open[ic] * 100.0f;
                hard_TP_condition = pc_gain > 15.0f;
            }
            else
            {
                timeout = false;
            }

            bool CLOSE_LONG_CONDI = PAIRS[ic].close[ii] > price_position_open[ic] + up * ATR_AT_OPEN[ic] || PAIRS[ic].close[ii] < price_position_open[ic] - down * ATR_AT_OPEN[ic] || timeout || hard_TP_condition;

            // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING POSITION AND ONLY THEN FOR OPENING POSITION

            // CLOSE LONG
            if (COIN_AMOUNTS[ic] > 0.0f && (CLOSE_LONG_CONDI || LAST_ITERATION))
            {
                const float to_add = COIN_AMOUNTS[ic] * PAIRS[ic].close[ii];
                USDT_amount += to_add;
                COIN_AMOUNTS[ic] = 0.0f;
                ATR_AT_OPEN[ic] = 0.0f;
                OPEN_TS[ic] = -10;

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
            if (COIN_AMOUNTS[ic] == 0.0f && OPEN_LONG_CONDI && LAST_ITERATION == false && ACTIVE_POSITIONS < MAX_OPEN_TRADES)
            {
                price_position_open[ic] = PAIRS[ic].close[ii];

                const float usdMultiplier = 1.0f / float(MAX_OPEN_TRADES - ACTIVE_POSITIONS);

                COIN_AMOUNTS[ic] = USDT_amount * usdMultiplier / PAIRS[ic].close[ii];
                USDT_amount -= USDT_amount * usdMultiplier;
                ATR_AT_OPEN[ic] = ATR[ic][ii];
                OPEN_TS[ic] = PAIRS[ic].timestamp[ii];

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
        if (closed || LAST_ITERATION || NEW_MONTH)
        {
            std::array<float, NB_PAIRS> closes{};
            for (uint ic = 0; ic < NB_PAIRS; ic++)
            {
                closes[ic] = PAIRS[ic].close[ii];
            }

            WALLET_VAL_USDT = USDT_amount + vector_product<NB_PAIRS>(COIN_AMOUNTS, closes);
            if (WALLET_VAL_USDT > MAX_WALLET_VAL_USDT)
            {
                MAX_WALLET_VAL_USDT = WALLET_VAL_USDT;
            }

            pc_change_with_max = (WALLET_VAL_USDT - MAX_WALLET_VAL_USDT) / MAX_WALLET_VAL_USDT * 100.0f;
            if (pc_change_with_max < max_drawdown)
            {
                max_drawdown = pc_change_with_max;
            }

            USDT_tracking.push_back(WALLET_VAL_USDT);
            USDT_tracking_ts.push_back(PAIRS[0].timestamp[ii]);
        }
    }

    std::array<float, NB_PAIRS> last_closes{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        last_closes[ic] = PAIRS[ic].close[nb_max - 1];
    }

    WALLET_VAL_USDT = USDT_amount + vector_product<NB_PAIRS>(COIN_AMOUNTS, last_closes);

    const float gain = (WALLET_VAL_USDT - USDT_amount_initial) / USDT_amount_initial * 100.0f;
    const float WR = float(nb_profit) / float(NB_POSI_ENTERED) * 100.0f;
    const float DDC = (1.0f / (1.0f + max_drawdown / 100.0f) - 1.0f) * 100.0f;
    const float score = gain / DDC * WR;

    i_print++;
    if (i_print == 11)
    {
        i_print = 0;
        print_best_res(best);
    }

    result.WALLET_VAL_USDT = USDT_amount;
    result.gain_over_DDC = gain / DDC;
    result.gain_pc = gain;
    result.max_DD = max_drawdown;
    result.nb_posi_entered = NB_POSI_ENTERED;
    result.win_rate = WR;
    result.score = score;
    result.calmar_ratio_monthly = calculate_calmar_ratio_monthly(USDT_tracking_ts, USDT_tracking, DDC);
    result.calmar_ratio = calculate_calmar_ratio(USDT_tracking_ts, USDT_tracking, DDC);
    result.ema1 = ema1;
    result.ema2 = ema2;
    result.ema3 = ema3;
    result.up = up;
    result.down = down;
    result.SRSIL = STOCH_RSI_LOWER;
    result.total_fees_paid = total_fees_paid_USDT;
    result.max_open_trades = MAX_OPEN_TRADES;
    result.param_str = "\n  EMA1: " + std::to_string(ema1) + " ; EMA2: " + std::to_string(ema2) + " ; EMA3: " + std::to_string(ema3) +
                       "\n  up: " + std::to_string(up) + " ; down: " + std::to_string(down) + " ; STOCH_RSI_LOWER: " + std::to_string(STOCH_RSI_LOWER);

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CALCULATE_SOME_INDICATORS(const std::vector<KLINEf> &PAIRS)
{
    std::vector<int> ema_values = combineAndRemoveDuplicates(range_EMA1, range_EMA2);
    ema_values = combineAndRemoveDuplicates(ema_values, range_EMA3);

    for (uint ic = 0; ic < COINS.size(); ic++)
    {
        std::cout << "Calculating some indicators for " << COINS[ic] << std::endl;

        for (const int ema_per : ema_values)
        {
            EMA_LISTS[ic]["EMA_" + std::to_string(ema_per)] = TALIB_EMA(PAIRS[ic].close, ema_per);
        }

        ATR[ic] = TALIB_ATR(PAIRS[ic].high, PAIRS[ic].low, PAIRS[ic].close, 14);
        StochRSI_K[ic] = TALIB_STOCHRSI_K(PAIRS[ic].close, 14, 14, 3, 3);
        StochRSI_D[ic] = TALIB_STOCHRSI_D(PAIRS[ic].close, 14, 14, 3, 3);

        std::cout << "Done." << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    const double t_begin = get_wall_time();
    std::cout << "\n-------------------------------------------------------------------------" << endl;
    std::cout << "Strategy to test: " << BLUE << STRAT_NAME << RESET << endl;
    std::cout << "DATA FILES TO PROCESS: " << endl;

    fill_datafile_paths();

    for (const string &dataf : DATAFILES)
    {
        std::cout << YELLOW << "  " << dataf << RESET << endl;
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

    std::vector<KLINEf> PAIRS;
    PAIRS.reserve(NB_PAIRS);
    for (const string &dataf : DATAFILES)
    {
        PAIRS.push_back(read_input_data(dataf));
    }

    start_indexes = INITIALIZE_DATA(PAIRS); // this function modifies PAIRS
    CALCULATE_SOME_INDICATORS(PAIRS);

    best.gain_over_DDC = -100.0f;
    best.calmar_ratio = -100.0f;
    best.calmar_ratio_monthly = -100.0f;
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
    std::cout << "Begin day                         : " << year << "/" << month << "/" << day << endl;
    std::cout << "End day                           : " << last_year << "/" << last_month << "/" << last_day << endl;
    std::cout << "Number of days                    : " << YELLOW << days << RESET << std::endl;
    std::cout << "OPEN/CLOSE FEE                    : " << FEE << " %" << endl;
    std::cout << "Minimum number of trades required : " << MIN_NUMBER_OF_TRADES << endl;
    std::cout << "Maximum drawback allowed          : " << MIN_ALLOWED_MAX_DRAWBACK << " %" << endl;
    std::cout << "EMA1 max tested                   : " << find_max(range_EMA1) << endl;
    std::cout << "EMA2 max tested                   : " << find_max(range_EMA2) << endl;
    std::cout << "EMA3 max tested                   : " << find_max(range_EMA3) << endl;
    std::cout << "StochRSI max tested               : " << find_max(range_STOCH_RSI_LOWER) << endl;
    std::cout << "-------------------------------------------------------------------------" << endl;

    // MAIN LOOP
    const long int nb_total = range_EMA1.size() * range_EMA2.size() * range_EMA3.size() * range_UP.size() * range_DOWN.size() * range_STOCH_RSI_LOWER.size() * MAX_OPEN_TRADES_TO_TEST.size();

    std::cout << "Running all backtests..." << std::endl;

    uint nb_done = 0;
    uint i_print3 = 0;

    RandomNumberGenerator rng1;

    for (long int ii = 0; ii < nb_total; ii++)
    {
        const int i1 = rng1.getRandomNumber(range_EMA1.size() - 1);
        const int i2 = rng1.getRandomNumber(range_EMA2.size() - 1);
        const int i3 = rng1.getRandomNumber(range_EMA3.size() - 1);
        const int i4 = rng1.getRandomNumber(range_UP.size() - 1);
        const int i5 = rng1.getRandomNumber(range_DOWN.size() - 1);
        const int i6 = rng1.getRandomNumber(range_STOCH_RSI_LOWER.size() - 1);
        const int i7 = rng1.getRandomNumber(MAX_OPEN_TRADES_TO_TEST.size() - 1);

        const EMA3_params para{range_EMA1[i1],
                               range_EMA2[i2],
                               range_EMA3[i3],
                               range_UP[i4],
                               range_DOWN[i5],
                               range_STOCH_RSI_LOWER[i6],
                               MAX_OPEN_TRADES_TO_TEST[i7]};

        if (para.ema3 < para.ema2 || para.ema3 < para.ema1 || para.ema2 < para.ema1 || std::abs(para.ema2 - para.ema1) < 5 || std::abs(para.ema3 - para.ema1) < 5 || std::abs(para.ema2 - para.ema3) < 5)
        {
            nb_done++;
            continue;
        }

        const RUN_RESULTf res = PROCESS(PAIRS, para.ema1, para.ema2, para.ema3, para.up, para.down, para.SRSIL, para.max_open_trades);

        if (res.score > best.score && res.gain_pc > 500.0f && res.gain_pc < 1000000.0f && res.nb_posi_entered >= MIN_NUMBER_OF_TRADES && res.max_DD > MIN_ALLOWED_MAX_DRAWBACK)
        {
            best = res;
        }

        nb_done++;
        i_print3++;

        if (i_print3 == 11)
        {

            WRITE_OR_UPDATE_BEST_SCORE_FILE(STRAT_NAME, out_filename, best);
            i_print3 = 0;
            const float pc_done = std::round(float(nb_done) / float(nb_total) * 100.0 * 100.0) / 100.0;
            std::cout << "DONE " << nb_done << " / " << nb_total << "   = " << pc_done << "%" << std::endl;
            std::cout << "Last done: " << res.param_str << "\n  " << i1 << " " << i2 << " " << i3 << " " << i4 << " " << i5 << " " << i6 << " " << i7 << std::endl;
        }
    }

    print_best_res(best);
    WRITE_OR_UPDATE_BEST_SCORE_FILE(STRAT_NAME, out_filename, best);

    const double t_end = get_wall_time();

    std::cout << "Number of backtests performed : " << nb_tested << endl;
    std::cout << "Time taken                    : " << t_end - t_begin << " seconds " << endl;
    const double ram_usage = process_mem_usage();
    std::cout << "RAM usage                     : " << std::round(ram_usage * 10.0) / 10.0 << " MB" << endl;
    std::cout << "-------------------------------------------------------------------------" << endl;

    TA_Shutdown();

    return 0;
}
