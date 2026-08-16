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

const string STRAT_NAME = "BBTREND";
const std::string out_filename = STRAT_NAME + "_best.txt";

static const uint NB_PAIRS = 11;
const vector<string> COINS = {"BTC", "ETH", "BNB", "XRP", "TRX", "MATIC", "LTC", "XMR", "XLM", "EOS", "ETC"};
const vector<uint> MAX_OPEN_TRADES_TO_TEST{5, 6, 7, 8, 9, 10};
const string timeframe = "2h";

vector<string> DATAFILES = {};

const float FEE = 0.1f;        // FEES in %
const float USDT_amount_initial = 1000.0f;
const uint MIN_NUMBER_OF_TRADES = 100;         // minimum number of trades required (to avoid some noise / lucky circunstances)
const float MIN_ALLOWED_MAX_DRAWBACK = -40.0f; // %

std::vector<uint> start_indexes{};

// RANGE OF EMA PERIDOS TO TESTs
const int period_max_EMA = 600;
const vector<int> range_EMA = integer_range(5, period_max_EMA + 2, 5);
const vector<int> range_BBlength = integer_range(5, 250, 5);
const vector<float> range_BBstd = float_Nvalues_range(0.25, 4.5, 10);
//////////////////////////
array<std::unordered_map<string, vector<float>>, NB_PAIRS> EMA_LISTS{};
array<vector<float>, NB_PAIRS> BollB_U{};
array<vector<float>, NB_PAIRS> BollB_M{};
array<vector<float>, NB_PAIRS> BollB_L{};

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

RUN_RESULTf PROCESS(const vector<KLINEf> &PAIRS, const int &ema_v, const int &BBlength, const float &BBstd, const uint &NB_POSITION_MAX)
{
    nb_tested++;

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        TALIB_BBANDS(PAIRS[ic].close, BBstd, BBstd, BBlength, BollB_U[ic], BollB_M[ic], BollB_L[ic]);
    }

    RUN_RESULTf result{};
    trade_core::WalletTrace wallet_trace{};
    trade_core::TradeStats stats{};
    trade_core::PortfolioState<NB_PAIRS> portfolio(USDT_amount_initial);
    array<float, NB_PAIRS> TSL_max_price_increase{};
    TSL_max_price_increase.fill(0.0f);

    const uint nb_max = PAIRS[0].nb;

    bool LAST_ITERATION = false;
    bool OPEN_LONG_CONDI = false;
    bool CLOSE_LONG_CONDI = false;

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
            if (ii < start_indexes[ic] + 2)
            {
                continue;
            }

            // conditions for open / close position

            // (dataframe['n1_close'] < dataframe['n1_higher_band'])
            //     &
            //     (dataframe['close'] > dataframe['higher_band'])
            //     &
            //     ((dataframe['n1_higher_band'] - dataframe['n1_lower_band']) / dataframe['n1_lower_band'] > self.min_bol_spread.value)
            //     &
            //     (dataframe['close'] > dataframe['long_ma'])

            // (
            //     (dataframe['close'] < dataframe['ma_band'])
            // ),
            // 'exit_long'] = 1

            OPEN_LONG_CONDI = PAIRS[ic].close[ii - 1] < BollB_U[ic][ii - 1] && PAIRS[ic].close[ii] > BollB_U[ic][ii] && PAIRS[ic].close[ii] > EMA_LISTS[ic]["EMA_" + std::to_string(ema_v)][ii];

            CLOSE_LONG_CONDI = PAIRS[ic].close[ii] < BollB_M[ic][ii];

            // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING POSITION AND ONLY THEN FOR OPENING POSITION

            // CLOSE LONG
            if (portfolio.coin_amounts[ic] > 0.0f && (CLOSE_LONG_CONDI || LAST_ITERATION))
            {
                trade_core::close_spot_long(portfolio, stats, ic, PAIRS[ic].close[ii], FEE);
                TSL_max_price_increase[ic] = 0.0f;
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
            {
                closes[ic] = PAIRS[ic].close[ii];
            }

            trade_core::record_spot_snapshot(portfolio, wallet_trace, closes, PAIRS[0].timestamp[ii]);
        }
    }

    array<float, NB_PAIRS> last_closes{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        last_closes[ic] = PAIRS[ic].close[nb_max - 1];
    }

    const double wallet_val_usdt = trade_core::calculate_spot_wallet_val_usdt(portfolio, last_closes);

    const trade_core::ResultMetrics metrics = trade_core::calculate_result_metrics(wallet_val_usdt, USDT_amount_initial, portfolio.max_drawdown, stats);

    trade_core::populate_common_result(result, metrics, wallet_val_usdt, portfolio.max_drawdown, portfolio.total_fees_paid_usdt, stats, NB_POSITION_MAX);
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
            EMA_LISTS[ic]["EMA_" + std::to_string(ema_per)] = TALIB_EMA(PAIRS[ic].close, ema_per);
        }
    }

    std::cout << "Done. " << endl;
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
    std::cout << "EMA period max tested    : " << find_max(range_EMA) << endl;

    std::vector<BBTREND_params> param_list{};
    param_list.reserve(range_EMA.size() * MAX_OPEN_TRADES_TO_TEST.size() * range_BBlength.size() * range_BBstd.size());
    for (const uint MAX_OPEN_TRADES : MAX_OPEN_TRADES_TO_TEST)
    {
        for (const int ema : range_EMA)
        {
            for (const int BBl : range_BBlength)
            {
                for (const float BBs : range_BBstd)
                {
                    param_list.push_back({ema, BBl, BBs, MAX_OPEN_TRADES});
                }
            }
        }
    }

    strategy_runner::SweepConfig cfg;
    cfg.strategy_name = STRAT_NAME;
    cfg.out_filename = out_filename;
    cfg.min_trades = MIN_NUMBER_OF_TRADES;
    cfg.min_dd = MIN_ALLOWED_MAX_DRAWBACK;
    cfg.print_every = 1000;

    strategy_runner::sweep(cfg, std::move(param_list), [&](const BBTREND_params &p) {
        return PROCESS(PAIRS, p.ema, p.BBlength, p.BBstd, p.max_open_trades);
    });

    strategy_runner::print_timing_and_ram(t_begin, nb_tested);
    TA_Shutdown();
    return 0;
}
