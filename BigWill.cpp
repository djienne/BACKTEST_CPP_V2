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

const string STRAT_NAME = "BigWill";
const string out_filename = STRAT_NAME + "_best.txt";

const vector<uint> MAX_OPEN_TRADES_TO_TEST{2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
const vector<string> COINS = {"BTC",
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
string timeframe = "1h";
vector<string> DATAFILES{};

const float start_year = 2017; // forced year to start (applies if data below is available)
const float FEE = 0.1f;        // FEES in %
const float USDT_amount_initial = 1000.0f;
const uint MIN_NUMBER_OF_TRADES = 100;         // minimum number of trades required (to avoid some noise / lucky circunstances)
const float MIN_ALLOWED_MAX_DRAWBACK = -42.0f; // %
const float stockOverBought = 0.800f;
const float stochOverSold = 0.200f;
const float WillOverSold = -85.0f;
const float WillOverBought = -10.0f;
const float HARD_TP_PC = 15.0f;
std::vector<uint> start_indexes{};

// RANGE OF PERIDOS TO TEST
// const int range_step = 2;
// vector<int> range_EMA = {180};
vector<int> range_AO_fast = integer_range(2, 102, 2);
vector<int> range_AO_slow = integer_range(2, 105, 5);
vector<int> range_EMA_fast = integer_range(2, 105, 5);
vector<int> range_EMA_slow = integer_range(50, 310, 10);
//////////////////////////
array<vector<float>, NB_PAIRS> AO{};
array<std::unordered_map<string, vector<float>>, NB_PAIRS> EMA_LISTS{};
array<vector<float>, NB_PAIRS> StochRSI{};
array<vector<float>, NB_PAIRS> WILLR{};

array<long int, NB_PAIRS> last_times{};

uint i_print = 0;
uint i_print3 = 0;
uint nb_tested = 0;
uint i_update_output = 0;

RUN_RESULTf best{};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fill_datafile_paths()
{
    for (uint i = 0; i < COINS.size(); i++)
    {
        DATAFILES.push_back("./data/data/binance/" + timeframe + "/" + COINS[i] + "-USDT.csv");
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_best_res(const RUN_RESULTf &bestt)
{
    std::cout << "\n--------------------------------------------------------------------------" << std::endl;
    std::cout << "BEST PARAMETER SET FOUND: " << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
    std::cout << "Time             : " << GREY << GET_CURRENT_TIME_STR() << RESET << std::endl;
    std::cout << "Strategy         : " << BLUE << STRAT_NAME << RESET << std::endl;
    std::cout << "Parameters       : " << YELLOW << bestt.param_str << RESET << std::endl;
    std::cout << "Max Open Trades  : " << YELLOW << bestt.max_open_trades << RESET << std::endl;
    std::cout << "Gain             : " << bestt.gain_pc << "%" << std::endl;
    std::cout << "Porfolio         : " << bestt.WALLET_VAL_USDT << "$ (started with 1000$)" << std::endl;
    std::cout << "Win rate         : " << bestt.win_rate << "%" << std::endl;
    std::cout << "max DD           : " << bestt.max_DD << "%" << std::endl;
    std::cout << "Gain/DDC         : " << bestt.gain_over_DDC << std::endl;
    std::cout << "Score            : " << GREEN << bestt.score << RESET << std::endl;
    std::cout << "Calmar ratio     : " << bestt.calmar_ratio << std::endl;
    std::cout << "Number of trades : " << bestt.nb_posi_entered << std::endl;
    std::cout << "Total fees paid  : " << round(bestt.total_fees_paid * 100.0f) / 100.0f << "$ (started with 1000$)" << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(const vector<KLINEf> &PAIRS, const int fast, const int slow, int ema_fast, int ema_slow, const uint MAX_OPEN_TRADES)
{
    nb_tested++;

    const std::string emaf_str = "EMA_" + std::to_string(ema_fast);
    const std::string emas_str = "EMA_" + std::to_string(ema_slow);

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
        AO[ic] = TALIB_AO(PAIRS[ic].high, PAIRS[ic].low, fast, slow);
    }
    uint ACTIVE_POSITIONS = 0;

    const uint ii_begin = start_indexes[0];

    for (uint ii = ii_begin + 1; ii < nb_max; ii++)
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

            bool TP_condition = false;
            if (COIN_AMOUNTS[ic] > 0)
            {
                const float pc_gain = (PAIRS[ic].high[ii] - price_position_open[ic]) / price_position_open[ic] * 100.0f;
                TP_condition = pc_gain > HARD_TP_PC;
            }

            OPEN_LONG_CONDI = EMA_LISTS[ic][emaf_str][ii] >= EMA_LISTS[ic][emas_str][ii] && WILLR[ic][ii] < WillOverSold && AO[ic][ii] > 0.0f && AO[ic][ii - 1] > AO[ic][ii];
            CLOSE_LONG_CONDI = (AO[ic][ii] < 0.0f && StochRSI[ic][ii] > stochOverSold) || WILLR[ic][ii] > WillOverBought || TP_condition;

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
            if (COIN_AMOUNTS[ic] == 0.0f && OPEN_LONG_CONDI && LAST_ITERATION == false && ACTIVE_POSITIONS < MAX_OPEN_TRADES)
            {
                price_position_open[ic] = PAIRS[ic].close[ii];

                const float usdMultiplier = 1.0f / float(MAX_OPEN_TRADES - ACTIVE_POSITIONS);

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
    if (i_print == 1000)
    {
        i_print = 0;
        std::cout << "DONE: Fast - Slow : " << fast << " - " << slow << std::endl;
        std::cout << "NB tested = " << nb_tested << "/" << range_AO_slow.size() * range_AO_fast.size() * MAX_OPEN_TRADES_TO_TEST.size() * range_EMA_fast.size() * range_EMA_slow.size() << std::endl;
        std::cout << "Done " << std::round(float(nb_tested) / float(range_AO_slow.size() * range_AO_fast.size() * MAX_OPEN_TRADES_TO_TEST.size() * range_EMA_fast.size() * range_EMA_slow.size()) * 100.0f * 100.0f) / 100.0f << " %" << std::endl;
        print_best_res(best);
    }

    result.WALLET_VAL_USDT = USDT_amount;
    result.gain_over_DDC = gain / DDC;
    result.gain_pc = gain;
    result.max_DD = max_drawdown;
    result.nb_posi_entered = NB_POSI_ENTERED;
    result.win_rate = WR;
    result.score = score;
    result.calmar_ratio = calculate_calmar_ratio(USDT_tracking_ts, USDT_tracking, DDC);
    result.ema1 = fast;
    result.ema2 = slow;
    result.ema3 = ema_fast;
    result.ema4 = ema_slow;
    result.total_fees_paid = total_fees_paid_USDT;
    result.max_open_trades = MAX_OPEN_TRADES;
    result.param_str = "\n  EMA_f: " + std::to_string(ema_fast) +
                       " ; EMA_s: " + std::to_string(ema_slow) +
                       " ; AO_f: " + std::to_string(fast) +
                       " ; AO_s: " + std::to_string(slow) +
                       "\n  stochOverSold: " + std::to_string(stochOverSold) +
                       " ; WillOverSold: " + std::to_string(WillOverSold) +
                       " ; WillOverBought: " + std::to_string(WillOverBought) +
                       "\n  HARD_TP: " + std::to_string(HARD_TP_PC);

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CALCULATE_INDICATORS(std::vector<KLINEf> &PAIRS)
{
    std::cout << "Calculating indicators..." << std::endl;

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        StochRSI[ic] = TALIB_STOCHRSI_not_averaged(PAIRS[ic].close, 14, 14);
    }
    cout << "Calculated STOCHRSI." << std::endl;

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        WILLR[ic] = TALIB_WILLR(PAIRS[ic].high, PAIRS[ic].low, PAIRS[ic].close, 14);
    }
    cout << "Calculated WILLR." << std::endl;

    std::vector<int> ema_values = combineAndRemoveDuplicates(range_EMA_fast, range_EMA_slow);
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        for (const int ema_per : ema_values)
        {
            EMA_LISTS[ic]["EMA_" + std::to_string(ema_per)] = TALIB_EMA(PAIRS[ic].close, ema_per);
        }
    }
    cout << "Calculated EMAs." << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    const double t_begin = get_wall_time();
    std::cout << "\n--------------------------------------------------------------------------" << std::endl;
    std::cout << "Strategy to test: " << STRAT_NAME << std::endl;
    std::cout << "DATA FILES TO PROCESS: " << std::endl;

    fill_datafile_paths();

    for (const string &dataf : DATAFILES)
    {
        std::cout << "  " << YELLOW << dataf << RESET << std::endl;
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
    std::cout << "Begin day      : " << year << "/" << month << "/" << day << std::endl;
    std::cout << "End day        : " << last_year << "/" << last_month << "/" << last_day << std::endl;
    std::cout << "Number of days : " << days << std::endl;
    std::cout << "OPEN/CLOSE FEE : " << FEE << " %" << std::endl;
    // std::cout << "Maximum number open trades: " << MAX_OPEN_TRADES << std::endl;
    std::cout << "Minimum number of trades required    : " << MIN_NUMBER_OF_TRADES << std::endl;
    std::cout << "Maximum drawback (=drawdown) allowed : " << MIN_ALLOWED_MAX_DRAWBACK << " %" << std::endl;
    std::cout << "fast period max tested : " << find_max(range_AO_fast) << std::endl;
    std::cout << "slow period max tested : " << find_max(range_AO_slow) << std::endl;
    std::cout << "StochRSI Upper Band   : " << stockOverBought << std::endl;
    std::cout << "StochRSI Lower Band   : " << stochOverSold << std::endl;
    std::cout << "WillR Over Bought     : " << WillOverBought << std::endl;
    std::cout << "WillR Over Sold       : " << WillOverSold << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;

    // MAIN LOOP
    std::vector<BigWill_params> param_list{};
    param_list.reserve(range_AO_slow.size() * range_AO_fast.size() * MAX_OPEN_TRADES_TO_TEST.size());

    for (const uint max_op_tr : MAX_OPEN_TRADES_TO_TEST)
    {
        for (const int fast : range_AO_fast)
        {
            for (const int slow : range_AO_slow)
            {
                for (const int ema_f : range_EMA_fast)
                {
                    for (const int ema_s : range_EMA_slow)
                    {
                        if (std::abs(fast - slow) < 7)
                            continue;

                        const BigWill_params to_add{fast, slow, ema_f, ema_s, max_op_tr};

                        param_list.push_back(to_add);
                    }
                }
            }
        }
    }

    std::cout << "Saved parameter list to test." << std::endl;
    std::cout << "Running all backtests..." << std::endl;

    random_shuffle_vector(param_list);

    uint nb_done = 0;

    for (const auto para : param_list)
    {
        const RUN_RESULTf res = PROCESS(PAIRS, para.AO_fast, para.AO_slow, para.ema_f, para.ema_s, para.max_open_trades);
        i_print3++;
        nb_done++;

        if (res.score > best.score && res.gain_pc > 50.0f && res.gain_pc < 1000000.0f && res.nb_posi_entered >= MIN_NUMBER_OF_TRADES && res.max_DD > MIN_ALLOWED_MAX_DRAWBACK)
        {
            best = res;
        }

        if (i_print3 == 1000)
        {
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
    std::cout << "--------------------------------------------------------------------------" << std::endl;

    TA_Shutdown();

    return 0;
}
