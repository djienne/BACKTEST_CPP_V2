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
#include "custom_talib_wrapper.hh"
#include "strategy_runner.hh"
#include <ta-lib/ta_libc.h>
using namespace std;
using uint = unsigned int;

const string STRAT_NAME = "SuperTrend_EMA_ATR";
const std::string out_filename = STRAT_NAME + "_best.txt";

// MAX_PAIRS sizes the fixed per-pair arrays; NB_PAIRS is the count actually
// traded, filled from the config in main().
using backtest_config::MAX_PAIRS;
uint NB_PAIRS = 0;
// Filled from backtest_config.json in main().
std::vector<std::string> COINS{};
// Filled in main() as 1..NB_PAIRS: the book cannot hold more than one
// position per pair, so larger values just repeat the same backtest.
std::vector<uint> MAX_OPEN_TRADES_TO_TEST{};
const vector<float> RANGE_UP{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
const vector<float> RANGE_DOWN{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
std::string timeframe{};
vector<string> DATAFILES = {};

const float FEE = 0.1f;        // FEES in %
const float USDT_amount_initial = 1000.0f;
const uint MIN_NUMBER_OF_TRADES = 100;         // minimum number of trades required (to avoid some noise / lucky circunstances)
const float MIN_ALLOWED_MAX_DRAWBACK = -40.0f; // %

std::vector<uint> start_indexes{};

// RANGE OF EMA PERIDOS TO TESTs
const int period_max_EMA = 600;
vector<int> range_EMA = integer_range(5, period_max_EMA + 2, 1);
//////////////////////////
array<std::unordered_map<string, vector<float>>, MAX_PAIRS> EMA_LISTS{};
array<vector<float>, MAX_PAIRS> ATR_LISTS{};
array<SuperTrend, MAX_PAIRS> BollB{};

uint nb_tested = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(const vector<KLINEf> &PAIRS, const int ema_v, const float UP, const float DOWN, const uint NB_POSITION_MAX)
{
    nb_tested++;

    RUN_RESULTf result{};

    trade_core::WalletTrace wallet_trace{};
    trade_core::TradeStats stats{};
    trade_core::PortfolioState<MAX_PAIRS> portfolio(USDT_amount_initial, NB_PAIRS);

    const uint nb_max = PAIRS[0].nb;

    bool LAST_ITERATION = false;
    bool OPEN_LONG_CONDI = false;
    bool CLOSE_LONG_CONDI = false;
    array<float, MAX_PAIRS> TSL_max_price_increase{};
    array<float, MAX_PAIRS> take_profit{};
    array<float, MAX_PAIRS> stop_loss{};
    array<float, MAX_PAIRS> stop_loss_at_open{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        TSL_max_price_increase[ic] = 0.0f;
    }

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

            if (portfolio.coin_amounts[ic] > 0.0f)
            {
                const float delta = PAIRS[ic].close[ii] - portfolio.price_position_open[ic];
                if (delta > TSL_max_price_increase[ic])
                {
                    TSL_max_price_increase[ic] = delta;
                    stop_loss[ic] = stop_loss_at_open[ic] + TSL_max_price_increase[ic];
                }
            }

            // conditions for open / close position

            const std::vector<float> &ema_series = EMA_LISTS[ic]["EMA_" + std::to_string(ema_v)];

            OPEN_LONG_CONDI = PAIRS[ic].close[ii] > ema_series[ii] && BollB[ic].supertrend[ii] == 1;

            // Exit on trend reversal, or on the ATR take-profit / stop-loss being reached
            // during the bar.
            //
            // The TP and SL comparisons were inverted: take-profit tested `close > TP`
            // and the stop tested `high < SL`. For a long, `high < SL` requires the whole
            // bar -- including its high -- to sit below the stop, and since high >= low it
            // is strictly harder to satisfy than the price actually touching the stop. The
            // stop therefore only fired once the market had already traded clean through
            // it, one or more bars late.
            //
            // Corrected to the intrabar convention this codebase already uses elsewhere
            // (see the CLOSE_LONG_CONDI in F_3EMA_SRSI_ATR.cpp): a long take-profit is
            // reached when the bar's high gets there, and a long stop when the bar's low
            // does.
            // TP/SL are only meaningful while a position is open: they are set on entry
            // and keep stale values afterwards. Gating on the position also avoids the
            // trap that a zeroed take_profit would make `high >= take_profit` trivially
            // true for every bar.
            const bool has_position = portfolio.coin_amounts[ic] > 0.0f;
            const bool trend_exit = PAIRS[ic].low[ii] < ema_series[ii] || BollB[ic].supertrend[ii] == -1;
            const bool take_profit_hit = has_position && PAIRS[ic].high[ii] >= take_profit[ic];
            const bool stop_loss_hit = has_position && PAIRS[ic].low[ii] <= stop_loss[ic];

            CLOSE_LONG_CONDI = trend_exit || take_profit_hit || stop_loss_hit;

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
                stop_loss[ic] = PAIRS[ic].close[ii] - ATR_LISTS[ic][ii] * DOWN;
                stop_loss_at_open[ic] = stop_loss[ic];
                take_profit[ic] = PAIRS[ic].close[ii] + ATR_LISTS[ic][ii] * UP;
            }
        }

        // check wallet status
        if (closed || LAST_ITERATION)
        {
            array<float, MAX_PAIRS> closes{};
            for (uint ic = 0; ic < NB_PAIRS; ic++)
            {
                closes[ic] = PAIRS[ic].close[ii];
            }

            trade_core::record_spot_snapshot(portfolio, wallet_trace, closes, PAIRS[0].timestamp[ii]);
        }
    }

    array<float, MAX_PAIRS> last_closes{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        last_closes[ic] = PAIRS[ic].close[nb_max - 1];
    }

    const double wallet_val_usdt = trade_core::calculate_spot_wallet_val_usdt(portfolio, last_closes);

    const trade_core::ResultMetrics metrics = trade_core::calculate_result_metrics(wallet_val_usdt, USDT_amount_initial, portfolio.max_drawdown, stats);

    trade_core::populate_common_result(result, metrics, wallet_val_usdt, portfolio.max_drawdown, portfolio.total_fees_paid_usdt, stats, NB_POSITION_MAX);
    result.calmar_ratio = calculate_calmar_ratio(wallet_trace.timestamps, wallet_trace.wallet_values, metrics.ddc);
    result.ema1 = ema_v;
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
    strategy_runner::init_talib();

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

    DATAFILES = backtest_config::spot_paths(CFG);
    random_shuffle_vector(range_EMA);

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

    std::vector<ST_EMA_ATR_params> param_list{};
    param_list.reserve(range_EMA.size() * RANGE_UP.size() * RANGE_DOWN.size() * MAX_OPEN_TRADES_TO_TEST.size());
    for (const uint MAX_OPEN_TRADES : MAX_OPEN_TRADES_TO_TEST)
    {
        for (const int ema : range_EMA)
        {
            for (const float up : RANGE_UP)
            {
                for (const float down : RANGE_DOWN)
                {
                    param_list.push_back({ema, up, down, MAX_OPEN_TRADES});
                }
            }
        }
    }

    strategy_runner::SweepConfig cfg;
    cfg.strategy_name = STRAT_NAME;
    cfg.out_filename = out_filename;
    cfg.min_trades = MIN_NUMBER_OF_TRADES;
    cfg.min_dd = MIN_ALLOWED_MAX_DRAWBACK;
    cfg.print_every = 500;

    strategy_runner::sweep(cfg, std::move(param_list), [&](const ST_EMA_ATR_params &p) {
        return PROCESS(PAIRS, p.ema, p.up, p.down, p.max_open_trades);
    });

    strategy_runner::print_timing_and_ram(t_begin, nb_tested);
    TA_Shutdown();
    return 0;
}
