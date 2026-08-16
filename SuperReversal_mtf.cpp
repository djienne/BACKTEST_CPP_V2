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
#include "strategy_runner.hh"
#include <ta-lib/ta_libc.h>
using namespace std;
using uint = unsigned int;

const std::string STRAT_NAME = "SuperReversal_mtf";
const std::string out_filename = STRAT_NAME + "_best.txt";

const std::vector<uint> MAX_OPEN_TRADES_TO_TEST{1, 2, 3, 4};
const std::vector<std::string> COINS = {"BTC", "ETH"};
static const uint NB_PAIRS = 2;

const std::string timeframe_1 = "1h";
const std::string timeframe = "5m";
const int int_htf = 60; // minutes
const int int_ltf = 5;  // minutes


const float start_year = 2017; // forced year to start (applies if data below is available)
const float FEE = 0.1f;        // FEES in %
const float USDT_amount_initial = 1000.0f;
const uint MIN_NUMBER_OF_TRADES = 200;         // minimum number of trades required (to avoid some noise / lucky circunstances)
const float MIN_ALLOWED_MAX_DRAWBACK = -40.0f; // %
std::vector<uint> start_indexes{};

// RANGE OF EMA PERIODS TO TEST
std::vector<int> range_ema_fast = integer_range(3, 200 + 4, 2);
std::vector<int> range_ema_slow = integer_range(70, 590, 5);
std::array<std::vector<float>, NB_PAIRS> SuperTrend_ltf{};
std::array<std::unordered_map<string, vector<float>>, NB_PAIRS> EMA_LISTS{};
//////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(vector<KLINEf> &PAIRS, const int &ema_f, const int &ema_s, const uint &MAX_OPEN_TRADES)
{
    RUN_RESULTf result{};

    trade_core::WalletTrace wallet_trace{};
    wallet_trace.wallet_values.reserve(1000);
    wallet_trace.timestamps.reserve(1000);
    trade_core::TradeStats stats{};
    trade_core::PortfolioState<NB_PAIRS> portfolio(USDT_amount_initial);

    const uint nb_max = PAIRS[0].nb;

    bool LAST_ITERATION = false;
    bool OPEN_LONG_CONDI = false;
    bool CLOSE_LONG_CONDI = false;

    const uint ii_begin = start_indexes[0];

    // Hoist EMA lookup keys out of the per-bar loop (was re-built 2x/bar previously).
    const std::string ema_f_str = "EMA_" + std::to_string(ema_f) + "_1h";
    const std::string ema_s_str = "EMA_" + std::to_string(ema_s) + "_1h";
    std::array<const std::vector<float> *, NB_PAIRS> ema_fast_ptr{};
    std::array<const std::vector<float> *, NB_PAIRS> ema_slow_ptr{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        ema_fast_ptr[ic] = &EMA_LISTS[ic].at(ema_f_str);
        ema_slow_ptr[ic] = &EMA_LISTS[ic].at(ema_s_str);
    }

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

            const std::vector<float> &ema_fast_v = *ema_fast_ptr[ic];
            const std::vector<float> &ema_slow_v = *ema_slow_ptr[ic];

            OPEN_LONG_CONDI = ema_fast_v[ii] > ema_slow_v[ii]
                                && SuperTrend_ltf[ic][ii] == 1
                                && PAIRS[ic].high[ii] > ema_fast_v[ii]
                                && PAIRS[ic].low[ii] < ema_fast_v[ii];

            CLOSE_LONG_CONDI = (ema_fast_v[ii] < ema_slow_v[ii] || SuperTrend_ltf[ic][ii] == -1)
                                && PAIRS[ic].high[ii] > ema_fast_v[ii]
                                && PAIRS[ic].low[ii] < ema_fast_v[ii];

            // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING POSITION AND ONLY THEN FOR OPENING POSITION

            // CLOSE LONG
            if (portfolio.coin_amounts[ic] > 0.0f && (CLOSE_LONG_CONDI || LAST_ITERATION))
            {
                trade_core::close_spot_long(portfolio, stats, ic, PAIRS[ic].close[ii], FEE);
                closed = true;
            }

            // OPEN LONG
            if (portfolio.coin_amounts[ic] == 0.0f && OPEN_LONG_CONDI && LAST_ITERATION == false && portfolio.active_positions < MAX_OPEN_TRADES)
            {
                trade_core::open_spot_long(portfolio, stats, ic, PAIRS[ic].close[ii], FEE, MAX_OPEN_TRADES);
            }
        }

        // check wallet status
        if (closed || LAST_ITERATION)
        {
            std::array<float, NB_PAIRS> closes{};
            for (uint ic = 0; ic < NB_PAIRS; ic++)
            {
                closes[ic] = PAIRS[ic].close[ii];
            }
            trade_core::record_spot_snapshot(portfolio, wallet_trace, closes, PAIRS[0].timestamp[ii]);
        }
    }

    std::array<float, NB_PAIRS> last_closes{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        last_closes[ic] = PAIRS[ic].close[nb_max - 1];
    }

    const double wallet_val_usdt = trade_core::calculate_spot_wallet_val_usdt(portfolio, last_closes);

    const trade_core::ResultMetrics metrics = trade_core::calculate_result_metrics(wallet_val_usdt, USDT_amount_initial, portfolio.max_drawdown, stats);

    trade_core::populate_common_result(result, metrics, wallet_val_usdt, portfolio.max_drawdown, portfolio.total_fees_paid_usdt, stats, MAX_OPEN_TRADES);
    result.calmar_ratio = calculate_calmar_ratio(wallet_trace.timestamps, wallet_trace.wallet_values, metrics.ddc);
    result.ema1 = ema_f;
    result.ema2 = ema_s;
    result.param_str = "\n  EMAf: " + std::to_string(ema_f) + " ; EMAs: " + std::to_string(ema_s);

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CALCULATE_INDICATORS(std::vector<KLINEf> &PAIRS, const int ltf_in_minutes, const int htf_in_minutes)
{
    std::cout << "Calculating indicators..." << std::endl;

    const int splitSize = htf_in_minutes / ltf_in_minutes;

    std::vector<int> ema_values = combineAndRemoveDuplicates(range_ema_fast, range_ema_slow);

    // Aggregate to the higher timeframe, compute there, then project back down.
    // RESAMPLE_TIMEFRAME / PROJECT_HTF_TO_LTF own the boundary alignment and the
    // anti-lookahead shift that this function used to hand-roll per indicator.
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        const Resampled htf = RESAMPLE_TIMEFRAME(PAIRS[ic], splitSize, ltf_in_minutes, htf_in_minutes);
        const size_t ltf_size = PAIRS[ic].close.size();

        /// Supertrend
        const std::vector<float> st_htf = TALIB_SuperTrend_dir_only(htf.kline.high, htf.kline.low, htf.kline.close, 15, 5.0f);
        SuperTrend_ltf[ic] = PROJECT_HTF_TO_LTF(st_htf, splitSize, ltf_size, htf.ltf_offset, -777.0f);

        /// EMAs
        for (const int ema_per : ema_values)
        {
            const std::string str = "EMA_" + std::to_string(ema_per) + "_1h";
            const std::vector<float> ema_htf = TALIB_EMA(htf.kline.close, ema_per);
            EMA_LISTS[ic][str] = PROJECT_HTF_TO_LTF(ema_htf, splitSize, ltf_size, htf.ltf_offset, -777.0f);
        }
    }

    std::cout << "Done calculating indicators." << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    const double t_begin = get_wall_time();
    strategy_runner::init_talib();

    const std::vector<std::string> DATAFILES = strategy_runner::build_spot_datafile_paths(COINS, timeframe);

    vector<KLINEf> PAIRS;
    PAIRS.reserve(NB_PAIRS);
    for (const string &dataf : DATAFILES)
    {
        PAIRS.push_back(read_input_data(dataf));
    }

    start_indexes = INITIALIZE_DATA(PAIRS);
    CALCULATE_INDICATORS(PAIRS, int_ltf, int_htf);

    strategy_runner::print_banner(STRAT_NAME, DATAFILES, PAIRS[0], FEE, MIN_NUMBER_OF_TRADES, MIN_ALLOWED_MAX_DRAWBACK);
    std::cout << "EMA short period max tested : " << find_max(range_ema_fast) << endl;
    std::cout << "EMA long period max tested  : " << find_max(range_ema_slow) << endl;

    std::vector<SR_params> param_list{};
    param_list.reserve(range_ema_slow.size() * range_ema_fast.size() * MAX_OPEN_TRADES_TO_TEST.size());
    for (const uint MAX_OPEN_TRADES : MAX_OPEN_TRADES_TO_TEST)
    {
        for (const int ema_s : range_ema_slow)
        {
            for (const int ema_f : range_ema_fast)
            {
                param_list.push_back({ema_f, ema_s, MAX_OPEN_TRADES});
            }
        }
    }
    std::cout << "Running all backtests..." << std::endl;

    strategy_runner::SweepConfig cfg;
    cfg.strategy_name = STRAT_NAME;
    cfg.out_filename = out_filename;
    cfg.min_trades = MIN_NUMBER_OF_TRADES;
    cfg.min_dd = MIN_ALLOWED_MAX_DRAWBACK;
    cfg.print_every = 100;

    uint nb_tested = 0;
    strategy_runner::sweep(cfg, param_list, [&](const SR_params &p) {
        nb_tested++;
        return PROCESS(PAIRS, p.ema_fast, p.ema_slow, p.max_open_trades);
    });

    strategy_runner::print_timing_and_ram(t_begin, nb_tested);
    TA_Shutdown();
    return 0;
}
