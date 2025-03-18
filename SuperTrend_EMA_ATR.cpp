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

const string STRAT_NAME = "SuperTrend_EMA_ATR";
const std::string out_filename = STRAT_NAME + "_best.txt";

static const uint NB_PAIRS = 11;
const vector<string> COINS = {"BTC", "ETH", "BNB", "XRP", "TRX", "MATIC", "LTC", "XMR", "XLM", "EOS", "ETC"};
const vector<uint> MAX_OPEN_TRADES_TO_TEST{5, 6, 7, 8, 9, 10};
const vector<float> RANGE_UP{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
const vector<float> RANGE_DOWN{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
const string timeframe = "4h";

vector<string> DATAFILES = {};

const float start_year = 2017; // forced year to start (applies if data below is available)
const float FEE = 0.1f;        // FEES in %
const float USDT_amount_initial = 1000.0f;
const uint MIN_NUMBER_OF_TRADES = 100;         // minimum number of trades required (to avoid some noise / lucky circunstances)
const float MIN_ALLOWED_MAX_DRAWBACK = -40.0f; // %

std::vector<uint> start_indexes{};

// RANGE OF EMA PERIDOS TO TESTs
const int period_max_EMA = 600;
vector<int> range_EMA = integer_range(5, period_max_EMA + 2, 1);
//////////////////////////
array<std::unordered_map<string, vector<float>>, NB_PAIRS> EMA_LISTS{};
array<vector<float>, NB_PAIRS> ATR_LISTS{};
array<SuperTrend, NB_PAIRS> BollB{};

uint i_print = 0;
uint nb_tested = 0;

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
    std::cout << "\n-------------------------------------" << endl;
    std::cout << "BEST PARAMETER SET FOUND: " << endl;
    std::cout << "-------------------------------------" << endl;
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
    std::cout << "-------------------------------------" << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(const vector<KLINEf> &PAIRS, const int ema_v, const float UP, const float DOWN, const uint NB_POSITION_MAX)
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
    array<float, NB_PAIRS> TSL_max_price_increase{};
    array<float, NB_PAIRS> price_position_open{};
    array<float, NB_PAIRS> COIN_AMOUNTS{};
    array<float, NB_PAIRS> take_profit{};
    array<float, NB_PAIRS> stop_loss{};
    array<float, NB_PAIRS> stop_loss_at_open{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        COIN_AMOUNTS[ic] = 0.0f;
        TSL_max_price_increase[ic] = 0.0f;
    }
    uint ACTIVE_POSITIONS = 0;

    const uint ii_begin = start_indexes[0];

    for (uint ii = ii_begin; ii < nb_max; ii++)
    {
        if (ii == nb_max - 1)
        {
            LAST_ITERATION = true;
        }

        bool closed = false;
        // For all pairs, check to close / open positions
        for (uint ic = 0; ic < NB_PAIRS; ic++)
        {
            if (ii < start_indexes[ic])
                continue;

            if (COIN_AMOUNTS[ic] > 0.0f)
            {
                const float delta = PAIRS[ic].close[ii] - price_position_open[ic];
                if (delta > TSL_max_price_increase[ic])
                {
                    TSL_max_price_increase[ic] = delta;
                    stop_loss[ic] = stop_loss_at_open[ic] + TSL_max_price_increase[ic];
                }
            }

            // conditions for open / close position

            OPEN_LONG_CONDI = PAIRS[ic].close[ii] > EMA_LISTS[ic]["EMA_" + std::to_string(ema_v)][ii] && BollB[ic].supertrend[ii] == 1;
            CLOSE_LONG_CONDI = PAIRS[ic].low[ii] < EMA_LISTS[ic]["EMA_" + std::to_string(ema_v)][ii] || BollB[ic].supertrend[ii] == -1 || PAIRS[ic].close[ii] > take_profit[ic] || PAIRS[ic].high[ii] < stop_loss[ic];

            // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING POSITION AND ONLY THEN FOR OPENING POSITION

            // CLOSE LONG
            if (COIN_AMOUNTS[ic] > 0.0f && (CLOSE_LONG_CONDI || LAST_ITERATION))
            {
                const float to_add = COIN_AMOUNTS[ic] * PAIRS[ic].close[ii];
                USDT_amount += to_add;
                COIN_AMOUNTS[ic] = 0.0f;
                TSL_max_price_increase[ic] = 0.0f;

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
                stop_loss[ic] = PAIRS[ic].close[ii] - ATR_LISTS[ic][ii] * DOWN;
                stop_loss_at_open[ic] = stop_loss[ic];
                take_profit[ic] = PAIRS[ic].close[ii] + ATR_LISTS[ic][ii] * UP;

                ACTIVE_POSITIONS++;
                NB_POSI_ENTERED++;
            }
        }

        // check wallet status
        if (closed || LAST_ITERATION)
        {
            array<float, NB_PAIRS> closes{};
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

    array<float, NB_PAIRS> last_closes{};
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
    if (i_print == 500)
    {
        i_print = 0;
        std::cout << "DONE: EMA: " << ema_v << endl;
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
    result.ema1 = ema_v;
    result.total_fees_paid = total_fees_paid_USDT;
    result.max_open_trades = NB_POSITION_MAX;
    result.param_str = "\n  EMA: " + std::to_string(ema_v) + "\n  UP: " + std::to_string(UP) + " ; DOWN: " + std::to_string(DOWN);

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CALCULATE_INDICATORS(std::vector<KLINEf> &PAIRS)
{
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        std::cout << "Calculating for " << COINS[ic] << endl;

        ATR_LISTS[ic] = TALIB_ATR(PAIRS[ic].high, PAIRS[ic].low, PAIRS[ic].close, 17);

        BollB[ic] = TALIB_SuperTrend(PAIRS[ic].high, PAIRS[ic].low, PAIRS[ic].close, 10, 3);

        for (const uint ema_per : range_EMA)
        {
            EMA_LISTS[ic]["EMA_" + std::to_string(ema_per)] = TALIB_EMA(PAIRS[ic].close, ema_per);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    const double t_begin = get_wall_time();
    fill_datafile_paths();
    std::cout << "\n-------------------------------------" << endl;
    std::cout << "Strategy to test: " << BLUE << STRAT_NAME << RESET << endl;
    std::cout << "DATA FILES TO PROCESS: " << endl;
    for (const string &dataf : DATAFILES)
    {
        std::cout << YELLOW << "  " << dataf << RESET << endl;
    }

    random_shuffle_vector(range_EMA);

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
    std::cout << "-------------------------------------" << endl;

    // MAIN LOOP

    std::vector<ST_EMA_ATR_params> param_list{};
    const int nb_tested = range_EMA.size() * RANGE_UP.size() * RANGE_DOWN.size() * MAX_OPEN_TRADES_TO_TEST.size();
    param_list.reserve(nb_tested);

    for (const uint MAX_OPEN_TRADES : MAX_OPEN_TRADES_TO_TEST)
    {
        for (const int ema : range_EMA)
        {
            for (const int up : RANGE_UP)
            {
                for (const int down : RANGE_DOWN)
                {
                    ST_EMA_ATR_params to_add{ema, up, down, MAX_OPEN_TRADES};
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

    for (const ST_EMA_ATR_params param : param_list)
    {
        const RUN_RESULTf res = PROCESS(PAIRS, param.ema, param.up, param.down, param.max_open_trades);
        i_print3++;
        nb_done++;

        if (res.score > best.score && res.gain_pc < 1000000.0f && res.nb_posi_entered >= MIN_NUMBER_OF_TRADES && res.max_DD > MIN_ALLOWED_MAX_DRAWBACK)
        {
            best = res;
        }

        if (i_print3 == 500)
        {
            WRITE_OR_UPDATE_BEST_SCORE_FILE(STRAT_NAME, out_filename, best);
            i_print3 = 0;
            const float pc_done = std::round(float(nb_done) / float(nb_tested) * 100.0 * 100.0) / 100.0;
            std::cout << "DONE " << nb_done << " / " << nb_tested << "   = " << pc_done << "%" << std::endl;
        }
    }

    print_best_res(best);

    const double t_end = get_wall_time();

    std::cout << "Number of backtests performed : " << nb_tested << endl;
    std::cout << "Time taken                    : " << t_end - t_begin << " seconds " << endl;
    const double ram_usage = process_mem_usage();
    std::cout << "RAM usage                     : " << std::round(ram_usage * 10.0) / 10.0 << " MB" << endl;
    std::cout << "-------------------------------------" << endl;

    TA_Shutdown();

    return 0;
}
