// ---------------------------------------------------------------------------
//  Strategy template -- a complete, working, minimal multi-pair spot strategy.
//
//  Copy this file, rename it, and edit the four marked sections. `make` discovers
//  *.cpp automatically, so there is no build list to update:
//
//      make MyStrategy && ./MyStrategy.exe
//
//  As written it is an RSI mean-reversion strategy: buy when RSI drops below a
//  threshold while price is above a long EMA (trend filter), sell when RSI
//  recovers past an upper threshold. It runs against the bundled 1h spot data.
//
//  See CONTRIBUTING.md for the conventions this file demonstrates.
// ---------------------------------------------------------------------------

#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "config.hh"
#include "indicators.hh"
#include "strategy_runner.hh"
#include "tools.hh"
#include "trade_core.hh"
#include <ta-lib/ta_libc.h>

using namespace std;

// --- 1. IDENTITY AND UNIVERSE ----------------------------------------------

// STRAT_NAME is the key looked up in backtest_config.json -- the coins and timeframe
// come from there, so a delisted symbol is a one-line config edit rather than an edit
// to every strategy that traded it.
const string STRAT_NAME = "TEMPLATE_RSI_EMA";
const string out_filename = STRAT_NAME + "_best.txt";

// Compile-time ceiling for the fixed-size per-pair arrays; the traded count below is
// filled from the config in main(). See the note on backtest_config::MAX_PAIRS.
using backtest_config::MAX_PAIRS;
uint NB_PAIRS = 0;

const float FEE = 0.1f;                        // per side, in %
const double USDT_amount_initial = 1000.0;
const uint MIN_NUMBER_OF_TRADES = 100;         // reject lucky, low-sample results
const float MIN_ALLOWED_MAX_DRAWBACK = -50.0f; // %

// --- 2. PARAMETER SPACE ----------------------------------------------------
// Everything swept goes here. Fixed settings stay as constants above.

const vector<int> range_rsi_period = integer_range(7, 29, 7);   // 7, 14, 21, 28
const vector<int> range_ema_period = integer_range(50, 350, 50);
const vector<float> range_rsi_buy = float_Nvalues_range(20.0f, 40.0f, 5);
const vector<float> range_rsi_sell = float_Nvalues_range(60.0f, 80.0f, 5);
// Filled in main() as 1..NB_PAIRS: sweeping more concurrent positions than there are
// pairs just repeats the same backtest.
vector<uint> MAX_OPEN_TRADES_TO_TEST{};

struct TemplateParams
{
    int rsi_period;
    int ema_period;
    float rsi_buy;
    float rsi_sell;
    uint max_open_trades;
};

vector<uint> start_indexes{};
uint nb_tested = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// --- 3. INDICATORS ---------------------------------------------------------
// Computed once, before the sweep, into each pair's IndicatorCache. Anything that
// depends on a swept parameter but not on all of them can instead be cached lazily
// inside PROCESS -- see the AO handling in BigWill.cpp.

// Worst warmup across every indicator computed below, so PROCESS can start its loop
// past the zero padding instead of guessing a constant.
size_t g_worst_warmup = 0;

void CALCULATE_INDICATORS(vector<KLINEf> &PAIRS)
{
    std::cout << "Calculating indicators..." << std::endl;

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        for (const int period : range_rsi_period)
        {
            size_t warmup = 0;
            PAIRS[ic].indicators.put(IndicatorCache::key("RSI", period),
                                     TALIB_RSI(PAIRS[ic].close, period, &warmup));
            g_worst_warmup = std::max(g_worst_warmup, warmup);
        }
        for (const int period : range_ema_period)
        {
            size_t warmup = 0;
            PAIRS[ic].indicators.put(IndicatorCache::key("EMA", period),
                                     TALIB_EMA(PAIRS[ic].close, period, &warmup));
            g_worst_warmup = std::max(g_worst_warmup, warmup);
        }
    }

    std::cout << "Done. Worst indicator warmup: " << g_worst_warmup << " bars." << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// --- 4. SIGNAL LOGIC -------------------------------------------------------
// One backtest for one parameter set.

RUN_RESULTf PROCESS(vector<KLINEf> &PAIRS, const TemplateParams &p)
{
    nb_tested++;

    trade_core::WalletTrace wallet_trace{};
    trade_core::TradeStats stats{};
    trade_core::PortfolioState<MAX_PAIRS> portfolio(USDT_amount_initial, NB_PAIRS);

    // Hoist the series out of the bar loop -- get() hashes a string, so calling it
    // per bar would dominate the runtime.
    array<const vector<float> *, MAX_PAIRS> RSI{};
    array<const vector<float> *, MAX_PAIRS> EMA{};
    const string rsi_key = IndicatorCache::key("RSI", p.rsi_period);
    const string ema_key = IndicatorCache::key("EMA", p.ema_period);
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        RSI[ic] = &PAIRS[ic].indicators.get(rsi_key);
        EMA[ic] = &PAIRS[ic].indicators.get(ema_key);
    }

    const uint nb_max = PAIRS[0].nb;
    // Start past both the indicator warmup and the pair's own history, rather than a
    // hard-coded constant that a larger indicator period could silently outgrow.
    const uint ii_begin = strategy_runner::first_tradable_index({g_worst_warmup}, start_indexes[0]);

    for (uint ii = ii_begin; ii < nb_max; ii++)
    {
        const bool LAST_ITERATION = (ii == nb_max - 1);
        bool closed = false;

        for (uint ic = 0; ic < NB_PAIRS; ic++)
        {
            // Each pair's history starts at a different index.
            if (ii < start_indexes[ic])
            {
                continue;
            }

            const vector<float> &rsi = *RSI[ic];
            const vector<float> &ema = *EMA[ic];
            const float close = PAIRS[ic].close[ii];

            const bool OPEN_LONG_CONDI = rsi[ii] < p.rsi_buy && close > ema[ii];
            const bool CLOSE_LONG_CONDI = rsi[ii] > p.rsi_sell;

            // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING AND ONLY THEN FOR OPENING:
            // a slot freed on this bar should be reusable on this bar.
            if (portfolio.coin_amounts[ic] > 0.0 && (CLOSE_LONG_CONDI || LAST_ITERATION))
            {
                trade_core::close_spot_long(portfolio, stats, ic, close, FEE);
                closed = true;
            }

            if (portfolio.coin_amounts[ic] == 0.0 && OPEN_LONG_CONDI && !LAST_ITERATION &&
                portfolio.active_positions < p.max_open_trades)
            {
                trade_core::open_spot_long(portfolio, stats, ic, close, FEE, p.max_open_trades);
            }
        }

        // Snapshot the equity curve whenever the book changes, which is what the
        // drawdown and Calmar figures are computed from.
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
    const trade_core::ResultMetrics metrics =
        trade_core::calculate_result_metrics(wallet_val_usdt, USDT_amount_initial, portfolio.max_drawdown, stats);

    RUN_RESULTf result{};
    trade_core::populate_common_result(result, metrics, wallet_val_usdt, portfolio.max_drawdown,
                                       portfolio.total_fees_paid_usdt, stats, p.max_open_trades);
    result.calmar_ratio = calculate_calmar_ratio(wallet_trace.timestamps, wallet_trace.wallet_values, metrics.ddc);
    // Parameters go in param_str rather than as new fields on RUN_RESULTf.
    result.param_str = "\n  RSI: " + std::to_string(p.rsi_period) +
                       " ; EMA: " + std::to_string(p.ema_period) +
                       "\n  buy<: " + std::to_string(p.rsi_buy) +
                       " ; sell>: " + std::to_string(p.rsi_sell);
    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    const double t_begin = get_wall_time();
    strategy_runner::init_talib();

    // Coins, market and timeframe come from backtest_config.json. load() aborts if this
    // strategy has no entry there, rather than guessing a universe nobody downloaded.
    const backtest_config::StrategyConfig CFG = backtest_config::load(STRAT_NAME);
    NB_PAIRS = CFG.nb_pairs();
    backtest_config::print_summary(CFG);

    const vector<string> DATAFILES = backtest_config::data_paths(CFG);

    vector<KLINEf> PAIRS;
    PAIRS.reserve(NB_PAIRS);
    for (const string &dataf : DATAFILES)
    {
        PAIRS.push_back(read_input_data(dataf));
    }

    // 1..NB_PAIRS concurrent positions; more than that is the same backtest repeated.
    for (uint n = 1; n <= NB_PAIRS; ++n)
    {
        MAX_OPEN_TRADES_TO_TEST.push_back(n);
    }

    // Aligns the pairs onto a common timestamp axis and returns each pair's start index.
    start_indexes = INITIALIZE_DATA(PAIRS);
    CALCULATE_INDICATORS(PAIRS);

    strategy_runner::print_banner(STRAT_NAME, DATAFILES, PAIRS[0], FEE, MIN_NUMBER_OF_TRADES, MIN_ALLOWED_MAX_DRAWBACK);

    vector<TemplateParams> param_list{};
    param_list.reserve(range_rsi_period.size() * range_ema_period.size() * range_rsi_buy.size() *
                       range_rsi_sell.size() * MAX_OPEN_TRADES_TO_TEST.size());
    for (const uint max_open : MAX_OPEN_TRADES_TO_TEST)
    {
        for (const int rsi_p : range_rsi_period)
        {
            for (const int ema_p : range_ema_period)
            {
                for (const float buy : range_rsi_buy)
                {
                    for (const float sell : range_rsi_sell)
                    {
                        param_list.push_back({rsi_p, ema_p, buy, sell, max_open});
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
    cfg.print_every = 200;

    // std::move: sweep() shuffles the list in place, so moving avoids copying it.
    strategy_runner::sweep(cfg, std::move(param_list), [&](const TemplateParams &p) {
        return PROCESS(PAIRS, p);
    });

    strategy_runner::print_timing_and_ram(t_begin, nb_tested);
    TA_Shutdown();
    return 0;
}
