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

const string STRAT_NAME = "TRIX";
const string out_filename = STRAT_NAME + "_best.txt";

static const uint NB_PAIRS = 11;
const vector<uint> MAX_OPEN_TRADES_TO_TEST{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
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

const string timeframe = "1h";

vector<string> DATAFILES = {};

const float start_year = 2017; // forced year to start (applies if data below is available)
const float FEE = 0.1f;        // FEES in %
const float USDT_amount_initial = 1000.0f;
const uint MIN_NUMBER_OF_TRADES = 200;
const float MIN_ALLOWED_MAX_DRAWBACK = -40.0f; // %
const float STOCH_RSI_UPPER = 0.800f;
const float STOCH_RSI_LOWER = 0.200f;

std::vector<uint> start_indexes{};

// RANGE OF EMA PERIODS TO TEST
const int period_max_EMA = 600;
const vector<int> range_EMA = integer_range(40, period_max_EMA + 2, 3); // best calmar from 40 to 122: 2.34
const vector<int> range_trixLength = integer_range(2, 100, 2);
const vector<int> range_trixSignal = integer_range(10, 100, 2);
//////////////////////////
array<std::unordered_map<string, vector<float>>, NB_PAIRS> EMA_LISTS{};
array<vector<float>, NB_PAIRS> StochRSI_LISTS{};

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

RUN_RESULTf PROCESS(const vector<KLINEf> &PAIRS, const int ema_v, const int trixLength_v, const int trixSignal_v, const uint MAX_OPEN_TRADESS)
{
    nb_tested++;

    RUN_RESULTf result{};

    trade_core::WalletTrace wallet_trace{};
    trade_core::TradeStats stats{};
    trade_core::PortfolioState<NB_PAIRS> portfolio(USDT_amount_initial);

    array<vector<float>, NB_PAIRS> TRIX_HISTO{};

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        TRIX_HISTO[ic] = TALIB_TRIX(PAIRS[ic].close, trixLength_v, trixSignal_v);
    }

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
        // For all paits, check to close and open positions
        for (uint ic = 0; ic < NB_PAIRS; ic++)
        {
            if (ii < start_indexes[ic])
                continue;

            // conditions for open / close position

            OPEN_LONG_CONDI = PAIRS[ic].close[ii] > EMA_LISTS[ic]["EMA_" + std::to_string(ema_v)][ii] && TRIX_HISTO[ic][ii] > 0.0f && StochRSI_LISTS[ic][ii] < STOCH_RSI_UPPER;
            CLOSE_LONG_CONDI = TRIX_HISTO[ic][ii] < 0.0f && StochRSI_LISTS[ic][ii] > STOCH_RSI_LOWER;

            // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING POSITION AND ONLY THEN FOR OPENING POSITION

            // CLOSE LONG
            if (portfolio.coin_amounts[ic] > 0.0f && (CLOSE_LONG_CONDI || LAST_ITERATION))
            {
                trade_core::close_spot_long(portfolio, stats, ic, PAIRS[ic].close[ii], FEE);
                closed = true;
            }

            // OPEN LONG
            if (portfolio.coin_amounts[ic] == 0.0f && OPEN_LONG_CONDI && LAST_ITERATION == false && portfolio.active_positions < MAX_OPEN_TRADESS)
            {
                trade_core::open_spot_long(portfolio, stats, ic, PAIRS[ic].close[ii], FEE, MAX_OPEN_TRADESS);
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
    {
        last_closes[ic] = PAIRS[ic].close[nb_max - 1];
    }

    const float wallet_val_usdt = trade_core::calculate_spot_wallet_val_usdt(portfolio, last_closes);

    const trade_core::ResultMetrics metrics = trade_core::calculate_result_metrics(wallet_val_usdt, USDT_amount_initial, portfolio.max_drawdown, stats);

    result.param_str = "\n  EMA: " + std::to_string(ema_v) + " ; trixLength: " + std::to_string(trixLength_v) + " ; trixSignal: " + std::to_string(trixSignal_v) +
                       "\n  STOCH_RSI_LOWER: " + std::to_string(STOCH_RSI_LOWER) + " ; STOCH_RSI_UPPER: " + std::to_string(STOCH_RSI_UPPER);

    trade_core::populate_common_result(result, metrics, wallet_val_usdt, portfolio.max_drawdown, portfolio.total_fees_paid_usdt, stats, MAX_OPEN_TRADESS);
    result.calmar_ratio = calculate_calmar_ratio(wallet_trace.timestamps, wallet_trace.wallet_values, metrics.ddc);
    result.ema1 = ema_v;

    return result;
}

void CALCULATE_INDICATORS(const std::vector<KLINEf> &PAIRS)
{
    for (uint ic = 0; ic < COINS.size(); ic++)
    {
        std::cout << "Calculating Indicators for " << COINS[ic] << std::endl;

        // StochRSI = TALIB_STOCHRSI_K(kline.d_close,14,3,3);
        StochRSI_LISTS[ic] = TALIB_STOCHRSI_not_averaged(PAIRS[ic].close, 14, 14);
        // std::cout << "Calculated STOCHRSI." << std::endl;

        for (const int ema_per : range_EMA)
        {
            EMA_LISTS[ic]["EMA_" + std::to_string(ema_per)] = TALIB_EMA(PAIRS[ic].close, ema_per);
        }
        // std::cout << "Calculated EMAs." << std::endl;
    }
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

    std::cout << "Total number of BTC KLines: " << PAIRS[0].timestamp.size() << std::endl;

    start_indexes = INITIALIZE_DATA(PAIRS);
    CALCULATE_INDICATORS(PAIRS);

    strategy_runner::print_banner(STRAT_NAME, DATAFILES, PAIRS[0], FEE, MIN_NUMBER_OF_TRADES, MIN_ALLOWED_MAX_DRAWBACK);
    std::cout << "StochRSI Upper Band      : " << STOCH_RSI_UPPER << std::endl;
    std::cout << "StochRSI Lower Band      : " << STOCH_RSI_LOWER << std::endl;
    std::cout << "EMA period max tested    : " << find_max(range_EMA) << std::endl;
    std::cout << "trixLength max tested    : " << find_max(range_trixLength) << std::endl;
    std::cout << "trixSignal max tested    : " << find_max(range_trixSignal) << std::endl;

    std::vector<trix_params> param_list{};
    param_list.reserve(range_EMA.size() * range_trixLength.size() * range_trixSignal.size() * MAX_OPEN_TRADES_TO_TEST.size());
    for (const uint MAX_OPEN_TRADES : MAX_OPEN_TRADES_TO_TEST)
    {
        for (const int ema : range_EMA)
        {
            for (const int trixL : range_trixLength)
            {
                for (const int trixS : range_trixSignal)
                {
                    param_list.push_back({ema, trixL, trixS, MAX_OPEN_TRADES});
                }
            }
        }
    }

    strategy_runner::SweepConfig cfg;
    cfg.strategy_name = STRAT_NAME;
    cfg.out_filename = out_filename;
    cfg.min_trades = MIN_NUMBER_OF_TRADES;
    cfg.min_dd = MIN_ALLOWED_MAX_DRAWBACK;
    cfg.print_every = 100;

    strategy_runner::sweep(cfg, param_list, [&](const trix_params &p) {
        return PROCESS(PAIRS, p.ema1, p.trixLength, p.trixSignal, p.max_open_trades);
    });

    strategy_runner::print_timing_and_ram(t_begin, nb_tested);
    TA_Shutdown();
    return 0;
}
