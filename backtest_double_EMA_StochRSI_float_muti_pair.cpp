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

const string STRAT_NAME = "2EMA_crossover_StochRSI";
const string out_filename = STRAT_NAME + "_best.txt";

static const uint NB_PAIRS = 2;
const uint NB_POSITION_MAX = 2;
const vector<string> COINS = {"BTC", "ETH"};
const string timeframe = "1h";

const vector<string> DATAFILES = {"./data/data/binance/" + timeframe + "/" + COINS[0] + "-USDT.csv",
                                  "./data/data/binance/" + timeframe + "/" + COINS[1] + "-USDT.csv"};

const float FEE = 0.1f;        // FEES in %
const float USDT_amount_initial = 1000.0f;
const uint MIN_NUMBER_OF_TRADES = 100;         // minimum number of trades required (to avoid some noise / lucky circunstances)
const float MIN_ALLOWED_MAX_DRAWBACK = -50.0f; // %
const float STOCH_RSI_UPPER = 0.800;
const float STOCH_RSI_LOWER = 0.200;

std::vector<uint> start_indexes{};

// RANGE OF EMA PERIODS TO TEST
const int period_max_EMA = 600;
vector<int> range_EMA_short = integer_range(2, 300 + 4, 1);
vector<int> range_EMA_long = integer_range(70, period_max_EMA + 4, 1);
//////////////////////////
array<std::unordered_map<string, vector<float>>, NB_PAIRS> EMA_LISTS{};
array<vector<float>, NB_PAIRS> StochRSI{};

uint nb_tested = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(const vector<KLINEf> &PAIRS, const int ema_s, const int ema_l)
{
    nb_tested++;

    RUN_RESULTf result{};

    trade_core::WalletTrace wallet_trace{};
    trade_core::TradeStats stats{};
    trade_core::PortfolioState<NB_PAIRS> portfolio(USDT_amount_initial);

    const uint nb_max = PAIRS[0].nb;

    bool LAST_ITERATION = false;
    bool OPEN_LONG_CONDI = false;
    bool CLOSE_LONG_CONDI = false;

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
                continue;

            // conditions for open / close position

            OPEN_LONG_CONDI = EMA_LISTS[ic]["EMA_" + std::to_string(ema_s)][ii] >= EMA_LISTS[ic]["EMA_" + std::to_string(ema_l)][ii] && StochRSI[ic][ii] < STOCH_RSI_UPPER;
            CLOSE_LONG_CONDI = EMA_LISTS[ic]["EMA_" + std::to_string(ema_s)][ii] <= EMA_LISTS[ic]["EMA_" + std::to_string(ema_l)][ii] && StochRSI[ic][ii] > STOCH_RSI_LOWER;

            // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING POSITION AND ONLY THEN FOR OPENING POSITION

            // CLOSE LONG
            if (portfolio.coin_amounts[ic] > 0.0f && (CLOSE_LONG_CONDI || LAST_ITERATION))
            {
                trade_core::close_spot_long(portfolio, stats, ic, PAIRS[ic].close[ii], FEE);
                closed = true;
            }

            // OPEN LONG
            if (portfolio.coin_amounts[ic] == 0.0f && OPEN_LONG_CONDI && LAST_ITERATION == false && portfolio.active_positions < NB_POSITION_MAX)
            {
                trade_core::open_spot_long(portfolio, stats, ic, PAIRS[ic].close[ii], FEE, NB_POSITION_MAX);
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

    result.param_str = "\n  EMA_short: " + std::to_string(ema_s) + " ; EMA_long: " + std::to_string(ema_l) +
                       "\n  STOCH_RSI_LOWER: " + std::to_string(STOCH_RSI_LOWER) + " ; STOCH_RSI_UPPER: " + std::to_string(STOCH_RSI_UPPER);
    trade_core::populate_common_result(result, metrics, WALLET_VAL_USDT, portfolio.max_drawdown, portfolio.total_fees_paid_usdt, stats, NB_POSITION_MAX);
    result.calmar_ratio = calculate_calmar_ratio(wallet_trace.timestamps, wallet_trace.wallet_values, metrics.ddc);
    result.ema1 = ema_s;
    result.ema2 = ema_l;

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CALCULATE_INDICATORS(std::vector<KLINEf> &PAIRS)
{
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        std::cout << "Calculating for " << COINS[ic] << std::endl;

        // std::cout << "Calculated STOCHRSI." << std::endl;

        std::vector<int> ema_values = combineAndRemoveDuplicates(range_EMA_long, range_EMA_short);

        for (const int ema_per : ema_values)
        {
            EMA_LISTS[ic]["EMA_" + std::to_string(ema_per)] = TALIB_EMA(PAIRS[ic].close, ema_per);
        }
        // std::cout << "Calculated EMAs." << std::endl;
    }
    cout << "Calculated EMAs." << std::endl;

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        StochRSI[ic] = TALIB_STOCHRSI_not_averaged(PAIRS[ic].close, 14, 14);
    }
    cout << "Calculated STOCHRSI." << std::endl;

    std::cout << "Initialized calculations." << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    const double t_begin = get_wall_time();
    strategy_runner::init_talib();

    vector<KLINEf> PAIRS;
    PAIRS.reserve(NB_PAIRS);
    for (const string &dataf : DATAFILES)
    {
        PAIRS.push_back(read_input_data(dataf));
    }

    random_shuffle_vector(range_EMA_long);

    start_indexes = INITIALIZE_DATA(PAIRS);
    CALCULATE_INDICATORS(PAIRS);

    strategy_runner::print_banner(STRAT_NAME, DATAFILES, PAIRS[0], FEE, MIN_NUMBER_OF_TRADES, MIN_ALLOWED_MAX_DRAWBACK);
    std::cout << "Maximum number open trades        : " << NB_POSITION_MAX << std::endl;
    std::cout << "StochRSI Upper Band               : " << STOCH_RSI_UPPER << std::endl;
    std::cout << "StochRSI Lower Band               : " << STOCH_RSI_LOWER << std::endl;
    std::cout << "EMA short period max tested       : " << find_max(range_EMA_short) << std::endl;
    std::cout << "EMA long period max tested        : " << find_max(range_EMA_long) << std::endl;

    std::vector<doubleEMA_params> param_list{};
    param_list.reserve(range_EMA_short.size() * range_EMA_long.size());
    for (const int ema1 : range_EMA_long)
    {
        for (const int ema2 : range_EMA_short)
        {
            if (std::abs(ema1 - ema2) < 5) continue;
            param_list.push_back({ema1, ema2});
        }
    }

    strategy_runner::SweepConfig cfg;
    cfg.strategy_name = STRAT_NAME;
    cfg.out_filename = out_filename;
    cfg.min_trades = MIN_NUMBER_OF_TRADES;
    cfg.min_dd = MIN_ALLOWED_MAX_DRAWBACK;
    cfg.print_every = 1000;

    strategy_runner::sweep(cfg, std::move(param_list), [&](const doubleEMA_params &p) {
        return PROCESS(PAIRS, p.ema1, p.ema2);
    });

    strategy_runner::print_timing_and_ram(t_begin, nb_tested);
    TA_Shutdown();
    return 0;
}
