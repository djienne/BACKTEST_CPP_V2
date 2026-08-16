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

// RANGE OF PERIODS TO TEST
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
uint nb_tested = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fill_datafile_paths()
{
    for (uint i = 0; i < COINS.size(); i++)
    {
        DATAFILES.push_back("./data/data/binance/" + timeframe + "/" + COINS[i] + "-USDT.csv");
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(const vector<KLINEf> &PAIRS, const int fast, const int slow, int ema_fast, int ema_slow, const uint MAX_OPEN_TRADES)
{
    nb_tested++;

    const std::string emaf_str = "EMA_" + std::to_string(ema_fast);
    const std::string emas_str = "EMA_" + std::to_string(ema_slow);

    RUN_RESULTf result{};

    trade_core::WalletTrace wallet_trace{};
    trade_core::TradeStats stats{};
    trade_core::PortfolioState<NB_PAIRS> portfolio(USDT_amount_initial);

    const uint nb_max = PAIRS[0].nb;

    bool LAST_ITERATION = false;
    bool OPEN_LONG_CONDI = false;
    bool CLOSE_LONG_CONDI = false;

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        AO[ic] = TALIB_AO(PAIRS[ic].high, PAIRS[ic].low, fast, slow);
    }

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
            if (portfolio.coin_amounts[ic] > 0)
            {
                const float pc_gain = (PAIRS[ic].high[ii] - portfolio.price_position_open[ic]) / portfolio.price_position_open[ic] * 100.0f;
                TP_condition = pc_gain > HARD_TP_PC;
            }

            OPEN_LONG_CONDI = EMA_LISTS[ic][emaf_str][ii] >= EMA_LISTS[ic][emas_str][ii] && WILLR[ic][ii] < WillOverSold && AO[ic][ii] > 0.0f && AO[ic][ii - 1] > AO[ic][ii];
            CLOSE_LONG_CONDI = (AO[ic][ii] < 0.0f && StochRSI[ic][ii] > stochOverSold) || WILLR[ic][ii] > WillOverBought || TP_condition;

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
            array<float, NB_PAIRS> closes{};
            for (uint ic = 0; ic < NB_PAIRS; ic++)
                closes[ic] = PAIRS[ic].close[ii];

            trade_core::record_spot_snapshot(portfolio, wallet_trace, closes, PAIRS[0].timestamp[ii]);
        }
    }

    array<float, NB_PAIRS> last_closes{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
        last_closes[ic] = PAIRS[ic].close[nb_max - 1];

    const double WALLET_VAL_USDT = trade_core::calculate_spot_wallet_val_usdt(portfolio, last_closes);

    const trade_core::ResultMetrics metrics = trade_core::calculate_result_metrics(WALLET_VAL_USDT, USDT_amount_initial, portfolio.max_drawdown, stats);

    trade_core::populate_common_result(result, metrics, WALLET_VAL_USDT, portfolio.max_drawdown, portfolio.total_fees_paid_usdt, stats, MAX_OPEN_TRADES);
    result.calmar_ratio = calculate_calmar_ratio(wallet_trace.timestamps, wallet_trace.wallet_values, metrics.ddc);
    result.ema1 = fast;
    result.ema2 = slow;
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
    strategy_runner::init_talib();

    fill_datafile_paths();

    vector<KLINEf> PAIRS;
    PAIRS.reserve(NB_PAIRS);
    for (const string &dataf : DATAFILES)
    {
        PAIRS.push_back(read_input_data(dataf));
    }

    start_indexes = INITIALIZE_DATA(PAIRS);
    CALCULATE_INDICATORS(PAIRS);

    strategy_runner::print_banner(STRAT_NAME, DATAFILES, PAIRS[0], FEE, MIN_NUMBER_OF_TRADES, MIN_ALLOWED_MAX_DRAWBACK);
    std::cout << "fast period max tested : " << find_max(range_AO_fast) << std::endl;
    std::cout << "slow period max tested : " << find_max(range_AO_slow) << std::endl;
    std::cout << "StochRSI Upper Band   : " << stockOverBought << std::endl;
    std::cout << "StochRSI Lower Band   : " << stochOverSold << std::endl;
    std::cout << "WillR Over Bought     : " << WillOverBought << std::endl;
    std::cout << "WillR Over Sold       : " << WillOverSold << std::endl;

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
                        if (std::abs(fast - slow) < 7) continue;
                        param_list.push_back({fast, slow, ema_f, ema_s, max_op_tr});
                    }
                }
            }
        }
    }

    strategy_runner::SweepConfig cfg;
    cfg.strategy_name = STRAT_NAME;
    cfg.out_filename = out_filename;
    cfg.min_trades = MIN_NUMBER_OF_TRADES;
    cfg.min_dd = MIN_ALLOWED_MAX_DRAWBACK;
    cfg.min_reasonable_gain = 50.0f;
    cfg.print_every = 1000;

    strategy_runner::sweep(cfg, param_list, [&](const BigWill_params &p) {
        return PROCESS(PAIRS, p.AO_fast, p.AO_slow, p.ema_f, p.ema_s, p.max_open_trades);
    });

    strategy_runner::print_timing_and_ram(t_begin, nb_tested);
    TA_Shutdown();
    return 0;
}
