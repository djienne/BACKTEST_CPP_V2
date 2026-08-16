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
#include "trade_core.hh"
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

    // Bollinger bands depend only on (BBlength, BBstd). This strategy samples the
    // parameter space at random, so only the current band set is held: there are ~4,800
    // reachable (length, std) combinations and caching all of them across 11 pairs would
    // run to gigabytes.
    const std::string bb_key = IndicatorCache::key("BBANDS", BBlength, BBstd);
    static std::string cached_bb_key;
    if (bb_key != cached_bb_key)
    {
        for (uint ic = 0; ic < NB_PAIRS; ic++)
        {
            df[ic].indicators.erase(cached_bb_key + ":U");
            df[ic].indicators.erase(cached_bb_key + ":M");
            df[ic].indicators.erase(cached_bb_key + ":L");

            BandsResult bands = TALIB_BBANDS_R(df[ic].close, BBstd, BBstd, BBlength);
            df[ic].indicators.put(bb_key + ":U", std::move(bands.upper));
            df[ic].indicators.put(bb_key + ":M", std::move(bands.middle));
            df[ic].indicators.put(bb_key + ":L", std::move(bands.lower));
        }
        cached_bb_key = bb_key;
    }

    std::array<const std::vector<float> *, NB_PAIRS> BB_U{};
    std::array<const std::vector<float> *, NB_PAIRS> BB_M{};
    std::array<const std::vector<float> *, NB_PAIRS> BB_L{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        BB_U[ic] = &df[ic].indicators.get(bb_key + ":U");
        BB_M[ic] = &df[ic].indicators.get(bb_key + ":M");
        BB_L[ic] = &df[ic].indicators.get(bb_key + ":L");
    }

    std::array<const std::vector<float> *, NB_PAIRS> EMA{};
    const std::string ema_key = IndicatorCache::key("EMA", ema_v);
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        EMA[ic] = &df[ic].indicators.get(ema_key);
    }

    RUN_RESULTf result{};
    trade_core::WalletTrace wallet_trace{};
    trade_core::TradeStats stats{};
    trade_core::PortfolioState<NB_PAIRS> portfolio(USDT_amount_initial);

    const uint nb_max = df[0].nb;

    bool LAST_ITERATION = false;
    bool OPEN_LONG_CONDI = false;
    bool OPEN_SHORT_CONDI = false;
    bool CLOSE_LONG_CONDI = false;
    bool CLOSE_SHORT_CONDI = false;

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

            trade_core::apply_funding_fee(portfolio, ic, df[ic].close[ii], funding_fee);

            const std::vector<float> &bb_u = *BB_U[ic];
            const std::vector<float> &bb_m = *BB_M[ic];
            const std::vector<float> &bb_l = *BB_L[ic];
            const std::vector<float> &ema = *EMA[ic];

            OPEN_LONG_CONDI = df[ic].close[ii - 1] < bb_u[ii - 1] && df[ic].close[ii] > bb_u[ii] && df[ic].close[ii] > ema[ii];
            OPEN_SHORT_CONDI = df[ic].close[ii - 1] > bb_l[ii - 1] && df[ic].close[ii] < bb_l[ii] && df[ic].close[ii] < ema[ii];
            CLOSE_LONG_CONDI = df[ic].close[ii] < bb_m[ii];
            CLOSE_SHORT_CONDI = df[ic].close[ii] > bb_m[ii];

            // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING POSITION AND ONLY THEN FOR OPENING POSITION

            // CLOSE LONG
            if (portfolio.coin_amounts[ic] > 0.0f && (CLOSE_LONG_CONDI || LAST_ITERATION))
            {
                trade_core::close_futures_long(portfolio, stats, ic, df[ic].close[ii], FEE);
                closed = true;
                portfolio.price_position_open[ic] = 0.0f;
            }

            // CLOSE SHORT
            if (CAN_SHORT && portfolio.coin_amounts[ic] < 0.0f && (CLOSE_SHORT_CONDI || LAST_ITERATION))
            {
                trade_core::close_futures_short(portfolio, stats, ic, df[ic].close[ii], FEE);
                closed = true;
                portfolio.price_position_open[ic] = 0.0f;
            }

            // OPEN LONG
            if (portfolio.coin_amounts[ic] == 0.0f && OPEN_LONG_CONDI && !LAST_ITERATION && portfolio.active_positions < MAX_OPEN_TRADES && portfolio.usdt_amount > 0.0f)
            {
                trade_core::open_futures_long(portfolio, stats, ic, df[ic].close[ii], FEE, MAX_OPEN_TRADES);
            }

            // OPEN SHORT
            if (CAN_SHORT && portfolio.coin_amounts[ic] == 0.0f && OPEN_SHORT_CONDI && !LAST_ITERATION && portfolio.active_positions < MAX_OPEN_TRADES && portfolio.usdt_amount > 0.0f)
            {
                trade_core::open_futures_short(portfolio, stats, ic, df[ic].close[ii], FEE, MAX_OPEN_TRADES);
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

            trade_core::record_futures_snapshot(portfolio, wallet_trace, current_closes, df[0].timestamp[ii]);
        }
    }

    std::array<float, NB_PAIRS> last_closes{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        last_closes[ic] = df[ic].close[nb_max - 1];
    }

    const double wallet_val_usdt = calculate_wallet_val_usdt<NB_PAIRS>(portfolio.usdt_amount, portfolio.coin_amounts, last_closes, portfolio.price_position_open);

    const trade_core::ResultMetrics metrics = trade_core::calculate_result_metrics(wallet_val_usdt, USDT_amount_initial, portfolio.max_drawdown, stats);

    trade_core::populate_common_result(result, metrics, wallet_val_usdt, portfolio.max_drawdown, portfolio.total_fees_paid_usdt, stats, MAX_OPEN_TRADES);
    result.calmar_ratio = calculate_calmar_ratio(wallet_trace.timestamps, wallet_trace.wallet_values, metrics.ddc);
    result.ema1 = ema_v;
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
            PAIRS[ic].indicators.put(IndicatorCache::key("EMA", ema_per), TALIB_EMA(PAIRS[ic].close, ema_per));
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

    std::time_t difference = std::abs(PAIRS[0].timestamp[last_idx] - PAIRS[0].timestamp[0]);
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
