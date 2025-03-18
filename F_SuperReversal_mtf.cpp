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

const std::string STRAT_NAME = "F_SuperReversal_mtf";
const std::string out_filename = STRAT_NAME + "_best.txt";

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
const bool CAN_SHORT = true;
const float ST_mult = 5.5;
const int ST_ATRper = 10;

const std::string timeframe_1 = "1h";
const std::string timeframe = "5m";
const int int_htf = 60; // minutes
const int int_ltf = 5;  // minutes

std::vector<std::string> DATAFILES{};
vector<string> DATAFILES_fundings{};

const float FEE = 0.1f; // FEES in %
const float USDT_amount_initial = 1000.0f;
const uint MIN_NUMBER_OF_TRADES = 20;          // minimum number of trades required (to avoid some noise / lucky circunstances)
const float MIN_ALLOWED_MAX_DRAWBACK = -40.0f; // %
std::vector<uint> start_indexes{};

// RANGE OF EMA PERIDOS TO TESTs
std::vector<int> range_ema_fast = generateRange_int(3, 600, 300);
std::vector<int> range_ema_slow = generateRange_int(3, 600, 300);
//////////////////////////

uint last_times[NB_PAIRS];

uint i_print = 0;
uint nb_tested = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_best_res(const RUN_RESULTf &bestt)
{
    std::cout << "\n--------------------------------------------------------------------------" << endl;
    std::cout << "BEST PARAMETER SET FOUND: " << endl;
    std::cout << "--------------------------------------------------------------------------" << endl;
    std::cout << "Time             : " << GREY << GET_CURRENT_TIME_STR() << RESET << std::endl;
    std::cout << "Strategy         : " << BLUE << STRAT_NAME << RESET << endl;
    std::cout << "Parameters       : " << YELLOW << bestt.param_str << RESET << std::endl;
    std::cout << "Max Open Trades  : " << bestt.max_open_trades << endl;
    std::cout << "Gain             : " << bestt.gain_pc << "%" << endl;
    std::cout << "Porfolio         : " << bestt.WALLET_VAL_USDT << "$ (started with 1000$)" << endl;
    std::cout << "Win rate         : " << bestt.win_rate << "%" << endl;
    std::cout << "max DD           : " << bestt.max_DD << "%" << endl;
    std::cout << "Gain/DDC         : " << bestt.gain_over_DDC << endl;
    std::cout << "Score            : " << GREEN << bestt.score << RESET << endl;
    std::cout << "Calmar ratio     : " << bestt.calmar_ratio << endl;
    std::cout << "Number of trades : " << bestt.nb_posi_entered << endl;
    std::cout << "Total fees paid  : " << round(bestt.total_fees_paid * 100.0f) / 100.0f << "$ (started with 1000$)" << endl;

    std::cout << "--------------------------------------------------------------------------" << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(const vector<KLINEf> &df, const std::vector<fundings> &FUNDINGS, const int &ema_f, const int &ema_s, const uint &MAX_OPEN_TRADES)
{
    nb_tested++;

    RUN_RESULTf result{};

    vector<float> USDT_tracking{};
    USDT_tracking.reserve(2000);
    vector<int> USDT_tracking_ts{};
    USDT_tracking_ts.reserve(2000);

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

            // conditions for open / close position
            const bool c_cross = df[ic].high[ii] > df[ic].EMA_1h[ema_f][ii] && df[ic].low[ii] < df[ic].EMA_1h[ema_f][ii];

            OPEN_LONG_CONDI = df[ic].EMA_1h[ema_f][ii] > df[ic].EMA_1h[ema_s][ii] && df[ic].SuperTrend_1h[ii] == 1 && c_cross;
            OPEN_SHORT_CONDI = df[ic].EMA_1h[ema_f][ii] < df[ic].EMA_1h[ema_s][ii] && df[ic].SuperTrend_1h[ii] == -1 && c_cross;

            CLOSE_LONG_CONDI = (df[ic].EMA_1h[ema_f][ii] < df[ic].EMA_1h[ema_s][ii] || df[ic].SuperTrend_1h[ii] == -1) && c_cross;
            CLOSE_SHORT_CONDI = (df[ic].EMA_1h[ema_f][ii] > df[ic].EMA_1h[ema_s][ii] || df[ic].SuperTrend_1h[ii] == 1) && c_cross;

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

                price_position_open[ic] = 0.0f;
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

                price_position_open[ic] = 0.0f;
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
                ACTIVE_POSITIONS++;
                NB_POSI_ENTERED++;

                price_position_open[ic] = df[ic].close[ii];
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
                ACTIVE_POSITIONS++;
                NB_POSI_ENTERED++;

                price_position_open[ic] = df[ic].close[ii];
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
    result.ema1 = ema_f;
    result.ema2 = ema_s;
    result.total_fees_paid = total_fees_paid_USDT;
    result.max_open_trades = MAX_OPEN_TRADES;
    result.param_str = "\n  EMAf: " + std::to_string(ema_f) + " ; EMAs: " + std::to_string(ema_s);

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CALCULATE_INDICATORS(std::vector<KLINEf> &PAIRS, const int ltf_in_minutes, const int htf_in_minutes)
{
    std::cout << "Calculating indicators..." << std::endl;

    const int splitSize = htf_in_minutes / ltf_in_minutes;

    const std::vector<int> ema_values = combineAndRemoveDuplicates(range_ema_fast, range_ema_slow);

    // calculate htf KLines from ltf KLines

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {

        const std::vector<std::vector<float>> s_open = splitVector(PAIRS[ic].open, splitSize);
        const std::vector<std::vector<float>> s_high = splitVector(PAIRS[ic].high, splitSize);
        const std::vector<std::vector<float>> s_low = splitVector(PAIRS[ic].low, splitSize);
        const std::vector<std::vector<float>> s_close = splitVector(PAIRS[ic].close, splitSize);

        std::vector<float> open_htf{};
        std::vector<float> high_htf{};
        std::vector<float> low_htf{};
        std::vector<float> close_htf{};

        const int nb_htf = s_close.size();
        open_htf.reserve(nb_htf);
        high_htf.reserve(nb_htf);
        low_htf.reserve(nb_htf);
        close_htf.reserve(nb_htf);

        for (size_t is = 0; is < s_close.size(); is++)
        {
            const int last_idx = s_close[is].size() - 1;
            open_htf.push_back(s_open[is][0]);
            close_htf.push_back(s_close[is][last_idx]);
            high_htf.push_back(*std::max_element(s_high[is].begin(), s_high[is].end()));
            low_htf.push_back(*std::min_element(s_low[is].begin(), s_low[is].end()));
        }

        /// Supertrend
        std::vector<float> SuperTrend = TALIB_SuperTrend_dir_only(high_htf, low_htf, close_htf, ST_ATRper, ST_mult);

        // resampling to lower time frame by duplicating elements and shifting values to the right (to avoid forward looking bias)
        SuperTrend = duplicateElements(SuperTrend, splitSize);
        SuperTrend.insert(SuperTrend.begin(), splitSize, 0);              // add 12 zeros a the beginning
        SuperTrend.erase(SuperTrend.end() - splitSize, SuperTrend.end()); // remove the 12 last elements

        if (std::abs(int(SuperTrend.size()) - int(PAIRS[ic].open.size())) > splitSize)
        {
            std::cout << "Error when calculating the htf supertrend into the ltf." << std::endl;
            std::cout << ic << " " << SuperTrend.size() << " " << PAIRS[ic].open.size() << std::endl;
            std::abort();
        }

        // remove end elements if necessary

        while (SuperTrend.size() != PAIRS[ic].open.size())
        {
            SuperTrend.pop_back();
        }

        PAIRS[ic].SuperTrend_1h = SuperTrend;

        /// EMAs
        for (const int ema_per : ema_values)
        {
            PAIRS[ic].EMA_1h[ema_per] = TALIB_EMA(close_htf, ema_per);

            // resampling to lower time frame by duplicating elements and shifting values to the right (to avoid forward looking bias)
            PAIRS[ic].EMA_1h[ema_per] = duplicateElements(PAIRS[ic].EMA_1h[ema_per], splitSize);
            PAIRS[ic].EMA_1h[ema_per].insert(PAIRS[ic].EMA_1h[ema_per].begin(), splitSize, 0);                             // add 12 zeros a the beginning
            PAIRS[ic].EMA_1h[ema_per].erase(PAIRS[ic].EMA_1h[ema_per].end() - splitSize, PAIRS[ic].EMA_1h[ema_per].end()); // remove the 12 last elements

            if (std::abs(int(PAIRS[ic].EMA_1h[ema_per].size()) - int(PAIRS[ic].open.size())) > splitSize)
            {
                std::cout << "Error when calculating the htf EMAs into the ltf." << std::endl;
                std::abort();
            }

            while (PAIRS[ic].EMA_1h[ema_per].size() != PAIRS[ic].open.size())
            {
                PAIRS[ic].EMA_1h[ema_per].pop_back();
            }
        }
    }

    std::cout << "Done calculating indicators." << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    const double t_begin = get_wall_time();
    std::cout << "\n--------------------------------------------------------------------------" << endl;
    std::cout << "Strategy to test: " << BLUE << STRAT_NAME << RESET << endl;
    std::cout << "DATA FILES TO PROCESS: " << endl;

    fill_datafile_paths_f(COINS, timeframe, DATAFILES, DATAFILES_fundings);

    RUN_RESULTf best{};

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

    vector<KLINEf> PAIRS;
    PAIRS.reserve(NB_PAIRS);
    for (const string &dataf : DATAFILES)
    {
        PAIRS.push_back(read_input_data_f(dataf, "2023-06-18"));
    }
    std::vector<fundings> FUNDINGS{};
    FUNDINGS.reserve(NB_PAIRS);
    for (const string &dataf : DATAFILES_fundings)
    {
        FUNDINGS.push_back(read_funding_rates_data(dataf));
    }

    start_indexes = INITIALIZE_DATA(PAIRS);
    CALCULATE_INDICATORS(PAIRS, int_ltf, int_htf);

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
    std::cout << "Begin day                   : " << year << "/" << month << "/" << day << endl;
    std::cout << "End day                     : " << last_year << "/" << last_month << "/" << last_day << endl;
    std::cout << "Number of days              : " << YELLOW << days << RESET << std::endl;
    std::cout << "Open/Close FEE              : " << FEE << " %" << endl;
    std::cout << "Minimum number of trades    : " << MIN_NUMBER_OF_TRADES << endl;
    std::cout << "Maximum drawdown allowed    : " << MIN_ALLOWED_MAX_DRAWBACK << " %" << endl;
    std::cout << "EMA short period max tested : " << find_max(range_ema_fast) << endl;
    std::cout << "EMA long period max tested  : " << find_max(range_ema_slow) << endl;
    std::cout << "--------------------------------------------------------------------------" << endl;

    // MAIN LOOP

    std::vector<SR_params> param_list{};
    param_list.reserve(range_ema_slow.size() * range_ema_fast.size() * MAX_OPEN_TRADES_TO_TEST.size());

    random_shuffle_vector(range_ema_slow);

    for (const uint MAX_OPEN_TRADES : MAX_OPEN_TRADES_TO_TEST)
    {
        for (const int ema_s : range_ema_slow)
        {
            for (const int ema_f : range_ema_fast)
            {
                SR_params to_add{ema_f, ema_s, MAX_OPEN_TRADES};
                param_list.push_back(to_add);
            }
        }
    }
    std::cout << "Saved parameter list to test." << std::endl;
    std::cout << "Running all backtests..." << std::endl;

    random_shuffle_vector(param_list);

    uint i_print3 = 0;
    uint nb_done = 0;

    for (const SR_params &par : param_list)
    {
        if (par.ema_fast > par.ema_slow)
        {
            nb_done++;
            continue;
        }
        const RUN_RESULTf res = PROCESS(PAIRS, FUNDINGS, par.ema_fast, par.ema_slow, par.max_open_trades);

        if (res.score > best.score && res.gain_pc < 1000000.0f && res.nb_posi_entered >= MIN_NUMBER_OF_TRADES && res.max_DD > MIN_ALLOWED_MAX_DRAWBACK)
        {
            best = res;
        }

        i_print3++;
        nb_done++;

        if (i_print3 == 100)
        {
            print_best_res(best);
            WRITE_OR_UPDATE_BEST_SCORE_FILE(STRAT_NAME, out_filename, best);
            i_print3 = 0;
            const float pc_done = std::round(float(nb_done) / float(param_list.size()) * 100.0 * 100.0) / 100.0;
            std::cout << "DONE " << nb_done << " / " << param_list.size() << "   = " << pc_done << "%" << std::endl;
        }
    }

    print_best_res(best);
    WRITE_OR_UPDATE_BEST_SCORE_FILE(STRAT_NAME, out_filename, best);

    const double t_end = get_wall_time();

    std::cout << "Number of backtests performed : " << nb_tested << endl;
    std::cout << "Time taken                    : " << t_end - t_begin << " seconds " << endl;
    const double ram_usage = process_mem_usage();
    std::cout << "RAM usage                     : " << std::round(ram_usage * 10.0) / 10.0 << " MB" << endl;
    std::cout << "--------------------------------------------------------------------------" << endl;

    TA_Shutdown();

    return 0;
}
