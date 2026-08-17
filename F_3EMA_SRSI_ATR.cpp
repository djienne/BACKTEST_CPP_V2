#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <math.h>
#include <unordered_map>
#include "config.hh"
#include "tools.hh"
#include "trade_core.hh"
#include "indicators.hh"
#include "strategy_runner.hh"
#include <ta-lib/ta_libc.h>
#include <iomanip>
#include <thread>
#include <random>

using namespace std;
using uint = unsigned int;

const string STRAT_NAME = "F_EMA3_SRSI_ATR";
const string out_filename = STRAT_NAME + "_best.txt";

// Filled in main() as 1..NB_PAIRS: the book cannot hold more than one
// position per pair, so larger values just repeat the same backtest.
std::vector<uint> MAX_OPEN_TRADES_TO_TEST{};
// Filled from backtest_config.json in main().
std::vector<std::string> COINS{};
// MAX_PAIRS sizes the fixed per-pair arrays; NB_PAIRS is the count actually
// traded, filled from the config in main().
using backtest_config::MAX_PAIRS;
uint NB_PAIRS = 0;
const bool CAN_SHORT = true;

std::string timeframe{};
const float FEE = 0.1f; // FEES in %
const float USDT_amount_initial = 1000.0f;

const uint MIN_NUMBER_OF_TRADES = 100; // minimum number of trades required (to avoid some noise / lucky circunstances)

const float MIN_ALLOWED_MAX_DRAWBACK = -50.0f; // %

// RANGE OF PARAMETERS TO TEST
const std::vector<int> range_EMA1 = generateRange_int(3, 100, 30);
const std::vector<int> range_EMA2 = generateRange_int(5, 400, 120);
const std::vector<int> range_EMA3 = generateRange_int(5, 580, 250);
const std::vector<float> range_UP = float_Nvalues_range(2.0, 20.0, 10);
const std::vector<float> range_DOWN = float_Nvalues_range(2.0, 20.0, 10);
const std::vector<float> range_STOCH_RSI_LOWER = float_Nvalues_range(0.09, 0.91, 12);
const std::vector<float> range_STOCH_RSI_UPPER = float_Nvalues_range(0.09, 0.91, 12);
//////////////////////////

uint nb_tested = 0;

// How many backtests between progress/score-file updates.
uint refresh_idx = 10;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_best_res(const RUN_RESULTf &bestt)
{
    std::cout << "\n-------------------------------------------------------------------------" << std::endl;
    std::cout << "BEST PARAMETER SET FOUND: " << std::endl;
    std::cout << "-------------------------------------------------------------------------" << std::endl;
    std::cout << "Thread               : " << std::this_thread::get_id() << std::endl;
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

RUN_RESULTf PROCESS(std::vector<KLINEf> &df, const std::vector<fundings> &FUNDINGS, const std::vector<uint> &start_indexes, const int &ema1, const int &ema2, const int &ema3, const float &up, const float &down, const float &STOCH_RSI_LOWER, const float &STOCH_RSI_UPPER, const uint &MAX_OPEN_TRADES)
{

    // The three EMAs depend only on (ema1, ema2, ema3). This strategy samples the
    // parameter space at random, so only the current triple is held: on 5m data a
    // single EMA series is ~2.8 MB per pair, and caching all ~400 sampled periods across
    // 11 pairs would run to double-digit gigabytes.
    //
    // thread_local, not static: mainProgramLogic runs on several worker threads, each
    // with its own copy of df.
    thread_local std::array<std::string, 3> cached_ema_keys;
    const std::array<std::string, 3> ema_keys = {IndicatorCache::key("EMA", ema1),
                                                 IndicatorCache::key("EMA", ema2),
                                                 IndicatorCache::key("EMA", ema3)};
    if (ema_keys != cached_ema_keys)
    {
        const std::array<int, 3> periods = {ema1, ema2, ema3};
        for (uint ic = 0; ic < NB_PAIRS; ic++)
        {
            for (const std::string &old_key : cached_ema_keys)
            {
                // Never evict a key the new triple still needs.
                if (!old_key.empty() && old_key != ema_keys[0] && old_key != ema_keys[1] && old_key != ema_keys[2])
                {
                    df[ic].indicators.erase(old_key);
                }
            }
            for (size_t k = 0; k < periods.size(); ++k)
            {
                if (!df[ic].indicators.has(ema_keys[k]))
                {
                    df[ic].indicators.put(ema_keys[k], TALIB_EMA(df[ic].close, periods[k]));
                }
            }
        }
        cached_ema_keys = ema_keys;
    }

    std::array<const std::vector<float> *, MAX_PAIRS> EMA1{};
    std::array<const std::vector<float> *, MAX_PAIRS> EMA2{};
    std::array<const std::vector<float> *, MAX_PAIRS> EMA3{};
    std::array<const std::vector<float> *, MAX_PAIRS> SRSI_K{};
    std::array<const std::vector<float> *, MAX_PAIRS> SRSI_D{};
    std::array<const std::vector<float> *, MAX_PAIRS> ATR{};
    const std::string k_key = IndicatorCache::key("STOCHRSI_K", 14, 14, 3, 3);
    const std::string d_key = IndicatorCache::key("STOCHRSI_D", 14, 14, 3, 3);
    const std::string atr_key = IndicatorCache::key("ATR", 14);
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        EMA1[ic] = &df[ic].indicators.get(ema_keys[0]);
        EMA2[ic] = &df[ic].indicators.get(ema_keys[1]);
        EMA3[ic] = &df[ic].indicators.get(ema_keys[2]);
        SRSI_K[ic] = &df[ic].indicators.get(k_key);
        SRSI_D[ic] = &df[ic].indicators.get(d_key);
        ATR[ic] = &df[ic].indicators.get(atr_key);
    }

    nb_tested++;

    RUN_RESULTf result{};
    trade_core::WalletTrace wallet_trace{};
    wallet_trace.wallet_values.reserve(2000);
    wallet_trace.timestamps.reserve(2000);
    trade_core::TradeStats stats{};
    trade_core::PortfolioState<MAX_PAIRS> portfolio(USDT_amount_initial, NB_PAIRS);

    const uint nb_max = df[0].close.size();

    bool LAST_ITERATION = false;
    std::array<float, MAX_PAIRS> ATR_AT_OPEN{};
    std::array<int64_t, MAX_PAIRS> OPEN_TS{};

    const uint ii_begin = start_indexes[0];

    bool NEW_MONTH = false;

    for (uint ii = ii_begin + 1; ii < nb_max; ii++)
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
            if (ii < start_indexes[ic] + 1)
            {
                continue;
            }

            // APPLY FUNDING FEES
            const float funding_fee = get_funding_fee_if_any(FUNDINGS[ic], df[ic].timestamp[ii]);

            trade_core::apply_funding_fee(portfolio, ic, df[ic].close[ii], funding_fee);

            // conditions for open / close position

            const std::vector<float> &e1 = *EMA1[ic];
            const std::vector<float> &e2 = *EMA2[ic];
            const std::vector<float> &e3 = *EMA3[ic];
            const std::vector<float> &sk = *SRSI_K[ic];
            const std::vector<float> &sd = *SRSI_D[ic];

            const bool OPEN_LONG_CONDI = e1[ii] > e2[ii] && e2[ii] > e3[ii] && df[ic].close[ii] > e1[ii] && sk[ii] <= STOCH_RSI_LOWER && sd[ii] <= STOCH_RSI_LOWER && sk[ii - 1] >= sd[ii - 1] && sk[ii] <= sd[ii];

            const bool OPEN_SHORT_CONDI = e1[ii] < e2[ii] && e2[ii] < e3[ii] && df[ic].close[ii] < e1[ii] && sk[ii] >= STOCH_RSI_UPPER && sd[ii] >= STOCH_RSI_UPPER && sk[ii - 1] <= sd[ii - 1] && sk[ii] >= sd[ii];

            bool timeout = false;
            if (portfolio.coin_amounts[ic] != 0.0f)
            {
                timeout = (df[ic].timestamp[ii] - OPEN_TS[ic]) >= 2 * 24 * 3600;
            }

            const bool CLOSE_LONG_CONDI = df[ic].high[ii] >= portfolio.price_position_open[ic] + up * ATR_AT_OPEN[ic] || df[ic].low[ii] <= portfolio.price_position_open[ic] - down * ATR_AT_OPEN[ic] || timeout;

            const bool CLOSE_SHORT_CONDI = df[ic].low[ii] <= portfolio.price_position_open[ic] - up * ATR_AT_OPEN[ic] || df[ic].high[ii] >= portfolio.price_position_open[ic] + down * ATR_AT_OPEN[ic] || timeout;

            // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING POSITION AND ONLY THEN FOR OPENING POSITION

            // CLOSE LONG
            if (portfolio.coin_amounts[ic] > 0.0f && (CLOSE_LONG_CONDI || LAST_ITERATION))
            {
                trade_core::close_futures_long(portfolio, stats, ic, df[ic].close[ii], FEE);
                closed = true;
                ATR_AT_OPEN[ic] = 0.0f;
                OPEN_TS[ic] = -10;
                portfolio.price_position_open[ic] = 0.0f;
            }

            // CLOSE SHORT
            if (CAN_SHORT && portfolio.coin_amounts[ic] < 0.0f && (CLOSE_SHORT_CONDI || LAST_ITERATION))
            {
                trade_core::close_futures_short(portfolio, stats, ic, df[ic].close[ii], FEE);
                closed = true;
                ATR_AT_OPEN[ic] = 0.0f;
                OPEN_TS[ic] = -10;
                portfolio.price_position_open[ic] = 0.0f;
            }

            // OPEN LONG
            if (portfolio.coin_amounts[ic] == 0.0f && OPEN_LONG_CONDI && !LAST_ITERATION && portfolio.active_positions < MAX_OPEN_TRADES && portfolio.usdt_amount > 0.0f)
            {
                trade_core::open_futures_long(portfolio, stats, ic, df[ic].close[ii], FEE, MAX_OPEN_TRADES);

                ATR_AT_OPEN[ic] = (*ATR[ic])[ii];
                OPEN_TS[ic] = df[ic].timestamp[ii];
            }

            // OPEN SHORT
            if (CAN_SHORT && portfolio.coin_amounts[ic] == 0.0f && OPEN_SHORT_CONDI && !LAST_ITERATION && portfolio.active_positions < MAX_OPEN_TRADES && portfolio.usdt_amount > 0.0f)
            {
                trade_core::open_futures_short(portfolio, stats, ic, df[ic].close[ii], FEE, MAX_OPEN_TRADES);

                ATR_AT_OPEN[ic] = (*ATR[ic])[ii];
                OPEN_TS[ic] = df[ic].timestamp[ii];
            }
        }

        // check wallet status
        if (closed || LAST_ITERATION || NEW_MONTH)
        {
            std::array<float, MAX_PAIRS> current_closes{};
            for (uint ic = 0; ic < NB_PAIRS; ic++)
            {
                current_closes[ic] = df[ic].close[ii];
            }

            trade_core::record_futures_snapshot(portfolio, wallet_trace, current_closes, df[0].timestamp[ii]);
        }
    }

    std::array<float, MAX_PAIRS> last_closes{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        last_closes[ic] = df[ic].close[nb_max - 1];
    }

    const double wallet_val_usdt = calculate_wallet_val_usdt<MAX_PAIRS>(portfolio.usdt_amount, portfolio.coin_amounts, last_closes, portfolio.price_position_open);

    const trade_core::ResultMetrics metrics = trade_core::calculate_result_metrics(wallet_val_usdt, USDT_amount_initial, portfolio.max_drawdown, stats);

    trade_core::populate_common_result(result, metrics, wallet_val_usdt, portfolio.max_drawdown, portfolio.total_fees_paid_usdt, stats, MAX_OPEN_TRADES);
    result.calmar_ratio_monthly = calculate_calmar_ratio_monthly(wallet_trace.timestamps, wallet_trace.wallet_values, metrics.ddc);
    result.calmar_ratio = calculate_calmar_ratio(wallet_trace.timestamps, wallet_trace.wallet_values, metrics.ddc);
    result.ema1 = ema1;
    result.ema2 = ema2;
    result.param_str = "\n  EMA1: " + std::to_string(ema1) + " ; EMA2: " + std::to_string(ema2) + " ; EMA3: " + std::to_string(ema3) +
                       "\n  up: " + std::to_string(up) + " ; down: " + std::to_string(down) + " ; STOCH_RSI_LOWER: " + std::to_string(STOCH_RSI_LOWER) + " ; STOCH_RSI_UPPER: " + std::to_string(STOCH_RSI_UPPER);

    // The EMAs stay in the cache rather than being freed here: they are keyed by period
    // and this sweep revisits the same periods across many parameter combinations, so
    // holding them is what turns a recompute-per-run into a compute-once.

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CALCULATE_SOME_INDICATORS(std::vector<KLINEf> &PAIRS)
{
    std::vector<int> ema_values = combineAndRemoveDuplicates(range_EMA1, range_EMA2);
    ema_values = combineAndRemoveDuplicates(ema_values, range_EMA3);

    for (uint ic = 0; ic < COINS.size(); ic++)
    {
        std::cout << "Calculating some indicators for " << COINS[ic] << std::endl;

        PAIRS[ic].indicators.put(IndicatorCache::key("ATR", 14),
                                 TALIB_ATR(PAIRS[ic].high, PAIRS[ic].low, PAIRS[ic].close, 14));
        PAIRS[ic].indicators.put(IndicatorCache::key("STOCHRSI_K", 14, 14, 3, 3),
                                 TALIB_STOCHRSI_K(PAIRS[ic].close, 14, 14, 3, 3));
        PAIRS[ic].indicators.put(IndicatorCache::key("STOCHRSI_D", 14, 14, 3, 3),
                                 TALIB_STOCHRSI_D(PAIRS[ic].close, 14, 14, 3, 3));

        std::cout << "Done." << std::endl;
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void mainProgramLogic()
{
    RandomNumberGenerator rng1;
    
    const double t_begin = get_wall_time();
    std::cout << "\n-------------------------------------------------------------------------" << endl;
    std::cout << "Strategy to test: " << BLUE << STRAT_NAME << RESET << endl;
    std::cout << "DATA FILES TO PROCESS: " << endl;
    
    std::vector<std::string> DATAFILES{};
    std::vector<std::string> DATAFILES_fundings{};

    // Coins, market and timeframe come from backtest_config.json; load() aborts if this
    // strategy has no entry, rather than guessing a universe nobody downloaded.
    const backtest_config::StrategyConfig CFG = backtest_config::load(STRAT_NAME);
    NB_PAIRS = CFG.nb_pairs();
    COINS = CFG.coins;
    timeframe = CFG.timeframe;
    backtest_config::print_summary(CFG);

    for (uint n = 1; n <= NB_PAIRS; ++n)
    {
        MAX_OPEN_TRADES_TO_TEST.push_back(n);
    }

    DATAFILES = backtest_config::futures_paths(CFG);
    DATAFILES_fundings = backtest_config::funding_paths(CFG);

    if (timeframe == "4h")
    {
        refresh_idx = 1000;
    }
    else if (timeframe == "1h")
    {
        refresh_idx = 100;
    }

    for (const string &dataf : DATAFILES)
    {
        std::cout << YELLOW << "  " << dataf << RESET << endl;
    }

    // TA-Lib is initialized once in main(), not here: this function runs on every
    // worker thread, and one thread's TA_Shutdown() used to land while another was
    // still computing.

    std::vector<KLINEf> PAIRS{};
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

    const std::vector<uint> start_indexes = INITIALIZE_DATA(PAIRS); // this function modifies PAIRS
    CALCULATE_SOME_INDICATORS(PAIRS);
    
    RUN_RESULTf best{};

    best.gain_over_DDC = -100.0f;
    best.calmar_ratio = -100.0f;
    best.calmar_ratio_monthly = -100.0f;
    best.score = -100.0f;

    const uint last_idx = PAIRS[0].close.size() - 1;

    const int year = get_year_from_timestamp(PAIRS[0].timestamp[0]);
    const int month = get_month_from_timestamp(PAIRS[0].timestamp[0]);
    const int day = get_day_from_timestamp(PAIRS[0].timestamp[0]);

    const int last_year = get_year_from_timestamp(PAIRS[0].timestamp[last_idx]);
    const int last_month = get_month_from_timestamp(PAIRS[0].timestamp[last_idx]);
    const int last_day = get_day_from_timestamp(PAIRS[0].timestamp[last_idx]);

    std::time_t difference = std::abs(PAIRS[0].timestamp[last_idx] - PAIRS[0].timestamp[0]);
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
    const long int nb_total = range_EMA1.size() * range_EMA2.size() * range_EMA3.size() * range_UP.size() * range_DOWN.size() * range_STOCH_RSI_LOWER.size() * range_STOCH_RSI_UPPER.size() * MAX_OPEN_TRADES_TO_TEST.size();

    std::cout << "Running all backtests..." << std::endl;

    uint nb_done = 0;
    uint i_print3 = 0;

    for (long int ii = 0; ii < nb_total; ii++)
    {
        const int i1 = rng1.getRandomNumber(range_EMA1.size() - 1);
        const int i2 = rng1.getRandomNumber(range_EMA2.size() - 1);
        const int i3 = rng1.getRandomNumber(range_EMA3.size() - 1);
        const int i4 = rng1.getRandomNumber(range_UP.size() - 1);
        const int i5 = rng1.getRandomNumber(range_DOWN.size() - 1);
        const int i6 = rng1.getRandomNumber(range_STOCH_RSI_LOWER.size() - 1);
        const int i7 = rng1.getRandomNumber(range_STOCH_RSI_UPPER.size() - 1);
        const int i8 = rng1.getRandomNumber(MAX_OPEN_TRADES_TO_TEST.size() - 1);

        const EMA3_params para{range_EMA1[i1],
                               range_EMA2[i2],
                               range_EMA3[i3],
                               range_UP[i4],
                               range_DOWN[i5],
                               range_STOCH_RSI_LOWER[i6],
                               range_STOCH_RSI_UPPER[i7],
                               MAX_OPEN_TRADES_TO_TEST[i8]};

        if (para.ema3 < para.ema2 || para.ema3 < para.ema1 || para.ema2 < para.ema1 || std::abs(para.ema2 - para.ema1) < 5 || std::abs(para.ema3 - para.ema1) < 5 || std::abs(para.ema2 - para.ema3) < 5)
        {
            nb_done++;
            continue;
        }

        const RUN_RESULTf res = PROCESS(PAIRS, FUNDINGS, start_indexes, para.ema1, para.ema2, para.ema3, para.up, para.down, para.SRSIL, para.SRSIU, para.max_open_trades);

        if (res.score > best.score && res.gain_pc < 1000000.0f && res.nb_posi_entered >= int(MIN_NUMBER_OF_TRADES) && res.max_DD > MIN_ALLOWED_MAX_DRAWBACK)
        {
            best = res;
        }

        nb_done++;
        i_print3++;

        if (i_print3 == refresh_idx)
        {
            print_best_res(best);
            WRITE_OR_UPDATE_BEST_SCORE_FILE(STRAT_NAME, out_filename, best);
            i_print3 = 0;
            const float pc_done = std::round(float(nb_done) / float(nb_total) * 100.0 * 100.0) / 100.0;
            std::cout << "DONE " << nb_done << " / " << nb_total << "   = " << pc_done << "%" << std::endl;
            std::cout << "Last done: " << res.param_str << "\n  " << i1 << " " << i2 << " " << i3 << " " << i4 << " " << i5 << " " << i6 << " " << i7 << " " << i8 << std::endl;
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

}

int main()
{
    const int N = 2; // Number of threads

    // TA-Lib is a process-wide library: initialize and shut it down exactly once,
    // around the workers rather than inside each of them.
    strategy_runner::init_talib();

    std::vector<std::thread> threads;
    threads.reserve(N);

    // Each worker samples the parameter space independently, so they only need distinct
    // RNG seeds -- which RandomNumberGenerator now provides by mixing in the thread id.
    for (int i = 0; i < N; ++i)
    {
        threads.emplace_back(mainProgramLogic);
    }

    for (auto &thread : threads)
    {
        thread.join();
    }

    TA_Shutdown();
    return 0;
}
