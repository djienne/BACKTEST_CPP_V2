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

const string STRAT_NAME = "F_BBTREND";
const std::string out_filename = STRAT_NAME + "_best.txt";

static const uint NB_PAIRS = 11;
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
const vector<uint> MAX_OPEN_TRADES_TO_TEST{5, 6, 7, 8, 9, 10, 11, 12};
const string timeframe = "2h";
const bool CAN_SHORT = true;

vector<string> DATAFILES = {};
vector<string> DATAFILES_fundings{};

const float start_year = 2017; // forced year to start (applies if data below is available)
const float FEE = 0.1f;        // FEES in %
const float USDT_amount_initial = 1000.0f;
const uint MIN_NUMBER_OF_TRADES = 100;         // minimum number of trades required (to avoid some noise / lucky circunstances)
const float MIN_ALLOWED_MAX_DRAWBACK = -50.0f; // %

std::vector<uint> start_indexes{};

// RANGE OF EMA PERIDOS TO TESTs
const int period_max_EMA = 600;
const vector<int> range_EMA = generateRange_int(5, period_max_EMA + 2, 120);
const vector<int> range_BBlength = generateRange_int(5, 400, 120);
const vector<float> range_BBstd = float_Nvalues_range(1.0, 4.5, 40);
//////////////////////////

uint i_print = 0;
uint nb_tested = 0;

RUN_RESULTf best{};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_best_res(const RUN_RESULTf &bestt)
{
    std::cout << "\n--------------------------------------------------------------------------" << endl;
    std::cout << "BEST PARAMETER SET FOUND: " << endl;
    std::cout << "--------------------------------------------------------------------------" << endl;
    std::cout << "Time             : " << GREY << GET_CURRENT_TIME_STR() << RESET << std::endl;
    std::cout << "Strategy         : " << BLUE << STRAT_NAME << RESET << endl;
    std::cout << "Parameters       : " << YELLOW << bestt.param_str << RESET << endl;
    std::cout << "Max Open Trades  : " << YELLOW << bestt.max_open_trades << RESET << endl;
    std::cout << "Best Gain        : " << bestt.gain_pc << "%" << endl;
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

RUN_RESULTf PROCESS(vector<KLINEf> &df, const std::vector<fundings> &FUNDINGS, const int &ema_v, const int &BBlength, const float &BBstd, const uint &MAX_OPEN_TRADES)
{
    nb_tested++;

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        TALIB_BBANDS(df[ic].close, BBstd, BBstd, BBlength, df[ic].BollB_U, df[ic].BollB_M, df[ic].BollB_L);
    }

    RUN_RESULTf result{};

    vector<float> USDT_tracking{};
    vector<int> USDT_tracking_ts{};

    const uint nb_max = df[0].nb;

    bool LAST_ITERATION = false;
    bool OPEN_LONG_CONDI = false;
    bool OPEN_SHORT_CONDI = false;
    bool CLOSE_LONG_CONDI = false;
    bool CLOSE_SHORT_CONDI = false;
    uint nb_profit = 0;
    uint nb_loss = 0;
    uint NB_POSI_ENTERED = 0;
    float pc_change_with_max = 0, max_drawdown = 0;
    float USDT_amount = USDT_amount_initial;
    float MAX_WALLET_VAL_USDT = USDT_amount_initial;
    float total_fees_paid_USDT = 0.0f;
    float WALLET_VAL_USDT = USDT_amount_initial;
    array<float, NB_PAIRS> TSL_max_price_increase{};
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

    bool NEW_MONTH = false;

    for (uint ii = ii_begin; ii < nb_max; ii++)
    {
        if (ii == nb_max - 1)
        {
            LAST_ITERATION = true;
        }

        const int month_b = get_month_from_timestamp(df[0].timestamp[ii - 1]);
        const int month = get_month_from_timestamp(df[0].timestamp[ii]);
        // std::cout << month << std::endl;

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

            OPEN_LONG_CONDI = df[ic].close[ii - 1] < df[ic].BollB_U[ii - 1] && df[ic].close[ii] > df[ic].BollB_U[ii] && df[ic].close[ii] > df[ic].EMA[ema_v][ii];
            OPEN_SHORT_CONDI = df[ic].close[ii - 1] > df[ic].BollB_L[ii - 1] && df[ic].close[ii] < df[ic].BollB_L[ii] && df[ic].close[ii] < df[ic].EMA[ema_v][ii];
            CLOSE_LONG_CONDI = df[ic].close[ii] < df[ic].BollB_M[ii];
            CLOSE_SHORT_CONDI = df[ic].close[ii] > df[ic].BollB_M[ii];

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
        if (closed || LAST_ITERATION || NEW_MONTH)
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
    result.ema1 = ema_v;
    result.total_fees_paid = total_fees_paid_USDT;
    result.max_open_trades = MAX_OPEN_TRADES;
    result.param_str = "\n  EMA: " + std::to_string(ema_v) + " ; BBlength: " + std::to_string(BBlength) + " ; BBstd: " + std::to_string(BBstd);

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CALCULATE_INDICATORS(std::vector<KLINEf> &PAIRS)
{

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        std::cout << "Calculating EMAs for " << COINS[ic] << endl;

        for (const uint ema_per : range_EMA)
        {
            PAIRS[ic].EMA[ema_per] = TALIB_EMA(PAIRS[ic].close, ema_per);
        }
    }

    std::cout << "Done. " << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    const double t_begin = get_wall_time();

    fill_datafile_paths_f(COINS, timeframe, DATAFILES, DATAFILES_fundings);

    std::cout << "\n--------------------------------------------------------------------------" << endl;
    std::cout << "Strategy to test: " << BLUE << STRAT_NAME << RESET << endl;
    std::cout << "DATA FILES TO PROCESS: " << endl;
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
    std::cout << "Begin day                : " << year << "/" << month << "/" << day << endl;
    std::cout << "End day                  : " << last_year << "/" << last_month << "/" << last_day << endl;
    std::cout << "Number of days           : " << YELLOW << days << RESET << std::endl;
    std::cout << "OPEN/CLOSE FEE           : " << FEE << " %" << endl;
    std::cout << "Minimum number of trades : " << MIN_NUMBER_OF_TRADES << endl;
    std::cout << "Maximum drawdown allowed : " << MIN_ALLOWED_MAX_DRAWBACK << " %" << endl;
    std::cout << "EMA period max tested    : " << find_max(range_EMA) << endl;
    std::cout << "--------------------------------------------------------------------------" << endl;

    // MAIN LOOP

    std::vector<BBTREND_params> param_list{};
    const int nb_tested = range_EMA.size() * MAX_OPEN_TRADES_TO_TEST.size() * range_BBlength.size() * range_BBstd.size();
    param_list.reserve(nb_tested);

    for (const uint MAX_OPEN_TRADES : MAX_OPEN_TRADES_TO_TEST)
    {
        for (const int ema : range_EMA)
        {
            for (const int BBl : range_BBlength)
            {
                for (const float BBs : range_BBstd)
                {
                    BBTREND_params to_add{ema, BBl, BBs, MAX_OPEN_TRADES};
                    param_list.push_back(to_add);
                }
            }
        }
    }
    std::cout << "Saved parameter list to test." << std::endl;
    std::cout << "Running all backtests..." << std::endl;

    random_shuffle_vector(param_list);

    int i_print3 = 0;
    int nb_done = 0;

    for (const BBTREND_params param : param_list)
    {
        const RUN_RESULTf res = PROCESS(PAIRS, FUNDINGS, param.ema, param.BBlength, param.BBstd, param.max_open_trades);
        i_print3++;
        nb_done++;

        if (res.score > best.score && res.gain_pc < 1000000.0f && res.nb_posi_entered >= MIN_NUMBER_OF_TRADES && res.max_DD > MIN_ALLOWED_MAX_DRAWBACK)
        {
            best = res;
        }

        if (i_print3 == 1000)
        {
            print_best_res(best);
            WRITE_OR_UPDATE_BEST_SCORE_FILE(STRAT_NAME, out_filename, best);
            i_print3 = 0;
            const float pc_done = std::round(float(nb_done) / float(nb_tested) * 100.0 * 100.0) / 100.0;
            std::cout << "DONE " << nb_done << " / " << nb_tested << "   = " << pc_done << "%" << std::endl;
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
