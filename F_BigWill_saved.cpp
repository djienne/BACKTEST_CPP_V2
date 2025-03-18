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

const std::string STRAT_NAME = "F_BigWill";
const std::string out_filename = STRAT_NAME + "_best.txt";

const bool CAN_SHORT = true;

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

const std::string timeframe = "1h";
std::vector<string> DATAFILES{};
std::vector<string> DATAFILES_fundings{};

const float start_year = 2017; // forced year to start (applies if data below is available)
const float FEE = 0.1f;        // FEES in %
const float USDT_amount_initial = 1000.0f;
const uint MIN_NUMBER_OF_TRADES = 100;         // minimum number of trades required (to avoid some noise / lucky circunstances)
const float MIN_ALLOWED_MAX_DRAWBACK = -80.0f; // %
const float stockOverBought = 0.700f;
const float stochOverSold = 0.300f;
const float WillOverSold = -85.0f;
const float WillOverBought = -10.0f;
const float HARD_TP_PC = 15.0f;
std::vector<uint> start_indexes{};

// RANGE OF PERIDOS TO TEST
// const int range_step = 2;
// vector<int> range_EMA = {180};
const std::vector<uint> MAX_OPEN_TRADES_TO_TEST{2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
const std::vector<int> range_AO_fast = integer_range(2, 102, 2);
const std::vector<int> range_AO_slow = integer_range(2, 105, 5);
const std::vector<int> range_EMA_fast = integer_range(2, 305, 5);
const std::vector<int> range_EMA_slow = integer_range(50, 610, 10);
//////////////////////////

std::array<long int, NB_PAIRS> last_times{};

uint i_print = 0;
uint i_print3 = 0;
uint nb_tested = 0;
uint i_update_output = 0;

RUN_RESULTf best{};

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

RUN_RESULTf PROCESS(std::vector<KLINEf> &df, const std::vector<fundings> &FUNDINGS, const int &AO_fast, const int &AO_slow, const int &ema_fast, const int &ema_slow, const uint &MAX_OPEN_TRADES)
{
    nb_tested++;

    RUN_RESULTf result{};

    std::vector<float> USDT_tracking{};
    std::vector<int> USDT_tracking_ts{};

    const uint nb_max = df[0].close.size();

    bool LAST_ITERATION = false;
    bool OPEN_LONG_CONDI = false;
    bool CLOSE_LONG_CONDI = false;
    bool OPEN_SHORT_CONDI = false;
    bool CLOSE_SHORT_CONDI = false;

    uint nb_profit = 0;
    uint nb_loss = 0;
    uint NB_POSI_ENTERED = 0;
    float pc_change_with_max = 0, max_drawdown = 0;
    float USDT_amount = USDT_amount_initial;
    float MAX_WALLET_VAL_USDT = USDT_amount_initial;
    float total_fees_paid_USDT = 0.0f;
    float WALLET_VAL_USDT = USDT_amount_initial;
    std::array<float, NB_PAIRS> price_position_open{};
    std::array<float, NB_PAIRS> COIN_AMOUNTS{};
    std::array<float, NB_PAIRS> take_profit{};
    std::array<float, NB_PAIRS> stop_loss{};
    std::array<float, NB_PAIRS> stop_loss_at_open{};

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        COIN_AMOUNTS[ic] = 0.0f;
        df[ic].AO = TALIB_AO(df[ic].high, df[ic].low, AO_fast, AO_slow);
    }

    uint ACTIVE_POSITIONS = 0;

    const uint ii_begin = start_indexes[0] + 2;

    for (uint ii = ii_begin + 1; ii < nb_max; ii++)
    {
        if (ii == nb_max - 1)
        {
            LAST_ITERATION = true;
        }

        bool closed = false;
        // For all pairs, check to close / open positions
        for (uint ic = 0; ic < NB_PAIRS; ic++)
        {
            if (ii < start_indexes[ic] + 2)
            {
                continue;
            }

            // APPLY FUNDING FEES
            const float funding_fee = get_funding_fee_if_any(FUNDINGS[ic], df[ic].timestamp[ii]);
            // std::cout << funding_fee << std::endl;

            if (funding_fee != 0.0f && COIN_AMOUNTS[ic] != 0.0f)
            {
                const float fe = COIN_AMOUNTS[ic] * df[ic].close[ii] * funding_fee; // positive funding fees means shorts get paid
                total_fees_paid_USDT += fe;
                USDT_amount -= fe;
            }

            bool TP_condition = false;
            if (COIN_AMOUNTS[ic] > 0.0f)
            {
                const float pc_gain = (df[ic].high[ii] - price_position_open[ic]) / price_position_open[ic] * 100.0f;
                TP_condition = pc_gain > HARD_TP_PC;
            }
            else if (COIN_AMOUNTS[ic] < 0.0f)
            {
                const float pc_gain = -1.0 * (df[ic].high[ii] - price_position_open[ic]) / price_position_open[ic] * 100.0f;
                TP_condition = pc_gain > HARD_TP_PC;
            }

            bool c01 = df[ic].AO[ii] >= 0.0f;
            bool c02 = df[ic].AO[ii] < df[ic].AO[ii - 1];
            bool c03 = df[ic].WILLR[ii] < WillOverSold;
            bool c04 = df[ic].EMA[ema_fast][ii] > df[ic].EMA[ema_slow][ii];
            OPEN_LONG_CONDI = c01 && c02 && c03 && c04;
            c01 = df[ic].AO[ii] <= 0.0f;
            c02 = df[ic].AO[ii] > df[ic].AO[ii - 1];
            c03 = df[ic].WILLR[ii] > WillOverBought;
            c04 = df[ic].EMA[ema_fast][ii] < df[ic].EMA[ema_slow][ii];
            OPEN_SHORT_CONDI = c01 && c02 && c03 && c04;

            bool c1 = (df[ic].AO[ii] < 0.0f || df[ic].WILLR[ii] > WillOverBought);
            bool c2 = (df[ic].StochRSI[ii] > stochOverSold || df[ic].WILLR[ii] > WillOverBought);
            CLOSE_LONG_CONDI = (c1 && c2) || TP_condition;
                 c1 = (df[ic].AO[ii] > 0.0f || df[ic].WILLR[ii] < WillOverSold);
                 c2 = (df[ic].StochRSI[ii] < stockOverBought || df[ic].WILLR[ii] < WillOverSold);
            CLOSE_SHORT_CONDI = (c1 && c2) || TP_condition;

            // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING POSITION AND ONLY THEN FOR OPENING POSITION

            // CLOSE LONG
            if (COIN_AMOUNTS[ic] > 0.0f && (CLOSE_LONG_CONDI || LAST_ITERATION))
            {
                const float to_add = COIN_AMOUNTS[ic] * df[ic].close[ii];
                USDT_amount += to_add;
                // apply FEEs
                const float fe = to_add * FEE / 100.0f;
                USDT_amount -= std::abs(fe);
                total_fees_paid_USDT += std::abs(fe);
                //
                if (df[ic].close[ii] >= price_position_open[ic])
                {
                    nb_profit++;
                }
                else
                {
                    nb_loss++;
                }
                closed = true;
                ACTIVE_POSITIONS--;
                COIN_AMOUNTS[ic] = 0.0f;
            }

            // CLOSE SHORT
            if (CAN_SHORT && COIN_AMOUNTS[ic] < 0.0f && (CLOSE_SHORT_CONDI || LAST_ITERATION))
            {
                const float to_add = std::abs(COIN_AMOUNTS[ic]) * (2.0f * price_position_open[ic] - df[ic].close[ii]);
                USDT_amount += to_add;
                // apply FEEs
                const float fe = to_add * FEE / 100.0f;
                USDT_amount -= std::abs(fe);
                total_fees_paid_USDT += std::abs(fe);
                //
                if (df[ic].close[ii] <= price_position_open[ic])
                {
                    nb_profit++;
                }
                else
                {
                    nb_loss++;
                }
                closed = true;
                ACTIVE_POSITIONS--;
                COIN_AMOUNTS[ic] = 0.0f;
            }

            // OPEN LONG
            if (COIN_AMOUNTS[ic] == 0.0f && OPEN_LONG_CONDI && !LAST_ITERATION && ACTIVE_POSITIONS < MAX_OPEN_TRADES && USDT_amount > 0.0f)
            {
                // order of lines is important
                const float usdMultiplier = 1.0f / float(MAX_OPEN_TRADES - ACTIVE_POSITIONS);
                COIN_AMOUNTS[ic] = USDT_amount * usdMultiplier / df[ic].close[ii];
                USDT_amount -= USDT_amount * usdMultiplier;
                // apply FEEs
                const float fe = std::abs(COIN_AMOUNTS[ic] * df[ic].close[ii] * FEE / 100.0f);
                USDT_amount -= fe;
                total_fees_paid_USDT += fe;
                //
                price_position_open[ic] = df[ic].close[ii];
                ACTIVE_POSITIONS++;
                NB_POSI_ENTERED++;
            }

            // OPEN SHORT
            if (CAN_SHORT && COIN_AMOUNTS[ic] == 0.0f && OPEN_SHORT_CONDI && !LAST_ITERATION && ACTIVE_POSITIONS < MAX_OPEN_TRADES && USDT_amount > 0.0f)
            {
                // order of lines is important
                const float usdMultiplier = 1.0f / float(MAX_OPEN_TRADES - ACTIVE_POSITIONS);
                COIN_AMOUNTS[ic] = -1.0f * USDT_amount * usdMultiplier / df[ic].close[ii];
                USDT_amount -= USDT_amount * usdMultiplier;
                // apply FEEs
                const float fe = std::abs(COIN_AMOUNTS[ic] * df[ic].close[ii] * FEE / 100.0f);
                USDT_amount -= fe;
                total_fees_paid_USDT += fe;
                //
                price_position_open[ic] = df[ic].close[ii];
                ACTIVE_POSITIONS++;
                NB_POSI_ENTERED++;
            }
        }

        // check wallet status
        if (closed || LAST_ITERATION)
        {
            std::array<float, NB_PAIRS> current_closes{};
            for (uint ic = 0; ic < NB_PAIRS; ic++)
            {
                current_closes[ic] = df[ic].close[ii];
            }
            WALLET_VAL_USDT = calculate_wallet_val_usdt<NB_PAIRS>(USDT_amount, COIN_AMOUNTS, current_closes, price_position_open);

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
            USDT_tracking_ts.push_back(df[0].timestamp[ii]);
        }
    }

    std::array<float, NB_PAIRS> last_closes{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        last_closes[ic] = df[ic].close[nb_max - 1];
    }
    WALLET_VAL_USDT = calculate_wallet_val_usdt<NB_PAIRS>(USDT_amount, COIN_AMOUNTS, last_closes, price_position_open);

    const float gain = (WALLET_VAL_USDT - USDT_amount_initial) / USDT_amount_initial * 100.0f;
    const float WR = float(nb_profit) / float(NB_POSI_ENTERED) * 100.0f;
    const float DDC = (1.0f / (1.0f + max_drawdown / 100.0f) - 1.0f) * 100.0f;
    const float score = gain / DDC * WR;

    result.WALLET_VAL_USDT = USDT_amount;
    result.gain_over_DDC = gain / DDC;
    result.gain_pc = gain;
    result.max_DD = max_drawdown;
    result.nb_posi_entered = NB_POSI_ENTERED;
    result.win_rate = WR;
    result.score = score;
    result.calmar_ratio = calculate_calmar_ratio(USDT_tracking_ts, USDT_tracking, DDC);
    result.ema1 = ema_fast;
    result.ema2 = ema_slow;
    result.AO_fast = ema_fast;
    result.AO_slow = ema_slow;
    result.total_fees_paid = total_fees_paid_USDT;
    result.max_open_trades = MAX_OPEN_TRADES;
    result.param_str = "\n  EMA_f: " + std::to_string(ema_fast) +
                       " ; EMA_s: " + std::to_string(ema_slow) +
                       " ; AO_f: " + std::to_string(AO_fast) +
                       " ; AO_s: " + std::to_string(AO_slow) +
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
        PAIRS[ic].StochRSI = TALIB_STOCHRSI_not_averaged(PAIRS[ic].close, 14, 14);
    }
    cout << "Calculated STOCHRSI." << std::endl;

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        PAIRS[ic].WILLR = TALIB_WILLR(PAIRS[ic].high, PAIRS[ic].low, PAIRS[ic].close, 14);
    }
    cout << "Calculated WILLR." << std::endl;

    std::vector<int> ema_values = combineAndRemoveDuplicates(range_EMA_fast, range_EMA_slow);
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        for (const int ema_per : ema_values)
        {
            PAIRS[ic].EMA[ema_per] = TALIB_EMA(PAIRS[ic].close, ema_per);
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

    fill_datafile_paths_f(COINS, timeframe, DATAFILES, DATAFILES_fundings);

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

    vector<KLINEf> PAIRS{};
    PAIRS.reserve(NB_PAIRS);
    for (const string &dataf : DATAFILES)
    {
        PAIRS.push_back(read_input_data_f(dataf));
    }
    std::vector<fundings> FUNDINGS{};
    FUNDINGS.reserve(NB_PAIRS);
    for (const string &dataf : DATAFILES_fundings)
    {
        FUNDINGS.push_back(read_funding_rates_data(dataf));
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

    const std::time_t difference = std::abs(int(PAIRS[0].timestamp[last_idx]) - int(PAIRS[0].timestamp[0]));
    const int days = difference / (24 * 60 * 60);

    // Display info
    std::cout << "Begin day      : " << year << "/" << month << "/" << day << std::endl;
    std::cout << "End day        : " << last_year << "/" << last_month << "/" << last_day << std::endl;
    std::cout << "Number of days : " << days << std::endl;
    std::cout << "OPEN/CLOSE FEE : " << FEE << " %" << std::endl;
    std::cout << "Minimum number of trades required    : " << MIN_NUMBER_OF_TRADES << std::endl;
    std::cout << "Maximum drawback (=drawdown) allowed : " << MIN_ALLOWED_MAX_DRAWBACK << " %" << std::endl;
    std::cout << "ema_fast max tested : " << find_max(range_EMA_fast) << std::endl;
    std::cout << "ema_slow max tested : " << find_max(range_EMA_slow) << std::endl;
    std::cout << "fast period max tested : " << find_max(range_AO_fast) << std::endl;
    std::cout << "slow period max tested : " << find_max(range_AO_slow) << std::endl;
    std::cout << "StochRSI Upper Band   : " << stockOverBought << std::endl;
    std::cout << "StochRSI Lower Band   : " << stochOverSold << std::endl;
    std::cout << "WillR Over Bought     : " << WillOverBought << std::endl;
    std::cout << "WillR Over Sold       : " << WillOverSold << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;

    // MAIN LOOP
    std::vector<BigWill_params> param_list{};

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
                        if (std::abs(fast - slow) < 7 || ema_f > ema_s || std::abs(ema_f - ema_s) < 10)
                        {
                            continue;
                        }

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

    for (const BigWill_params para : param_list)
    {
        const RUN_RESULTf res = PROCESS(PAIRS, FUNDINGS, para.AO_fast, para.AO_slow, para.ema_f, para.ema_s, para.max_open_trades);
        i_print3++;
        nb_done++;

        if (res.score > best.score && res.gain_pc > 15.0 && res.gain_pc < 1000000.0f && res.nb_posi_entered >= MIN_NUMBER_OF_TRADES && res.max_DD > MIN_ALLOWED_MAX_DRAWBACK)
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
    std::cout << "--------------------------------------------------------------------------" << std::endl;

    TA_Shutdown();

    return 0;
}
