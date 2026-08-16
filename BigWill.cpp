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
// Indicator series live in each pair's IndicatorCache (see Klinef.hh), not in
// strategy-level globals.
uint nb_tested = 0;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(vector<KLINEf> &PAIRS, const int fast, const int slow, int ema_fast, int ema_slow, const uint MAX_OPEN_TRADES)
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

    // AO depends only on (fast, slow) -- roughly 1,050 distinct pairs -- but the sweep
    // also varies ema_fast, ema_slow and max_open_trades, so recomputing it here ran the
    // same full-series work millions of times over.
    //
    // main() orders the parameter list so all entries sharing a (fast, slow) pair run
    // consecutively, which means holding just the current pair's series is enough: AO is
    // computed ~1,050 times instead of ~5,000,000. Caching every combination instead
    // would work too, but at ~2.2 MB per pair it would grow past 2 GB.
    static std::string cached_ao_key;
    const std::string ao_key = IndicatorCache::key("AO", fast, slow);
    if (ao_key != cached_ao_key)
    {
        for (uint ic = 0; ic < NB_PAIRS; ic++)
        {
            PAIRS[ic].indicators.erase(cached_ao_key);
            PAIRS[ic].indicators.put(ao_key, TALIB_AO(PAIRS[ic].high, PAIRS[ic].low, fast, slow));
        }
        cached_ao_key = ao_key;
    }

    std::array<const std::vector<float> *, NB_PAIRS> AO{};
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        AO[ic] = &PAIRS[ic].indicators.get(ao_key);
    }

    // Hoist the per-pair series out of the bar loop. These used to be looked up through
    // a string-keyed map on every bar for every pair.
    std::array<const std::vector<float> *, NB_PAIRS> EMA_F{};
    std::array<const std::vector<float> *, NB_PAIRS> EMA_S{};
    std::array<const std::vector<float> *, NB_PAIRS> SRSI{};
    std::array<const std::vector<float> *, NB_PAIRS> WILL{};
    const std::string emaf_key = IndicatorCache::key("EMA", ema_fast);
    const std::string emas_key = IndicatorCache::key("EMA", ema_slow);
    const std::string srsi_key = IndicatorCache::key("STOCHRSI", 14, 14);
    const std::string willr_key = IndicatorCache::key("WILLR", 14);
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        EMA_F[ic] = &PAIRS[ic].indicators.get(emaf_key);
        EMA_S[ic] = &PAIRS[ic].indicators.get(emas_key);
        SRSI[ic] = &PAIRS[ic].indicators.get(srsi_key);
        WILL[ic] = &PAIRS[ic].indicators.get(willr_key);
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
                const double pc_gain = (PAIRS[ic].high[ii] - portfolio.price_position_open[ic]) / portfolio.price_position_open[ic] * 100.0;
                TP_condition = pc_gain > HARD_TP_PC;
            }

            const std::vector<float> &ao = *AO[ic];
            const std::vector<float> &willr = *WILL[ic];
            const std::vector<float> &srsi = *SRSI[ic];

            OPEN_LONG_CONDI = (*EMA_F[ic])[ii] >= (*EMA_S[ic])[ii] && willr[ii] < WillOverSold && ao[ii] > 0.0f && ao[ii - 1] > ao[ii];
            CLOSE_LONG_CONDI = (ao[ii] < 0.0f && srsi[ii] > stochOverSold) || willr[ii] > WillOverBought || TP_condition;

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
        PAIRS[ic].indicators.put(IndicatorCache::key("STOCHRSI", 14, 14),
                                 TALIB_STOCHRSI_not_averaged(PAIRS[ic].close, 14, 14));
    }
    cout << "Calculated STOCHRSI." << std::endl;

    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        PAIRS[ic].indicators.put(IndicatorCache::key("WILLR", 14),
                                 TALIB_WILLR(PAIRS[ic].high, PAIRS[ic].low, PAIRS[ic].close, 14));
    }
    cout << "Calculated WILLR." << std::endl;

    std::vector<int> ema_values = combineAndRemoveDuplicates(range_EMA_fast, range_EMA_slow);
    for (uint ic = 0; ic < NB_PAIRS; ic++)
    {
        for (const int ema_per : ema_values)
        {
            PAIRS[ic].indicators.put(IndicatorCache::key("EMA", ema_per), TALIB_EMA(PAIRS[ic].close, ema_per));
        }
    }
    cout << "Calculated EMAs." << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    const double t_begin = get_wall_time();
    strategy_runner::init_talib();

    // Shared path helper rather than a per-strategy copy of the same loop.
    DATAFILES = strategy_runner::build_spot_datafile_paths(COINS, timeframe);

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

    // The parameter list is built grouped by (AO_fast, AO_slow), because AO is the one
    // expensive indicator here and depends only on that pair. Keeping a group contiguous
    // lets PROCESS hold a single AO series at a time -- ~1,050 computations instead of
    // ~5,000,000, without the ~2.2 GB a cache-everything approach would need.
    //
    // The group order is shuffled (and each group internally too), so an interrupted run
    // has still sampled the whole space -- the property the sweep's own shuffle provides.
    // cfg.shuffle is therefore turned off, since reshuffling would break the grouping.
    std::vector<std::pair<int, int>> ao_groups{};
    ao_groups.reserve(range_AO_fast.size() * range_AO_slow.size());
    for (const int fast : range_AO_fast)
    {
        for (const int slow : range_AO_slow)
        {
            if (std::abs(fast - slow) < 7) continue;
            ao_groups.emplace_back(fast, slow);
        }
    }
    random_shuffle_vector(ao_groups);

    std::vector<BigWill_params> inner{};
    inner.reserve(range_EMA_fast.size() * range_EMA_slow.size() * MAX_OPEN_TRADES_TO_TEST.size());

    std::vector<BigWill_params> param_list{};
    // Reserve the whole product. This used to omit both EMA ranges, under-reserving by
    // a factor of ~550 and forcing repeated reallocation of a list that grows past five
    // million entries.
    param_list.reserve(ao_groups.size() * range_EMA_fast.size() * range_EMA_slow.size() *
                       MAX_OPEN_TRADES_TO_TEST.size());
    for (const auto &g : ao_groups)
    {
        inner.clear();
        for (const uint max_op_tr : MAX_OPEN_TRADES_TO_TEST)
        {
            for (const int ema_f : range_EMA_fast)
            {
                for (const int ema_s : range_EMA_slow)
                {
                    inner.push_back({g.first, g.second, ema_f, ema_s, max_op_tr});
                }
            }
        }
        random_shuffle_vector(inner, false); // quiet: this runs once per group
        param_list.insert(param_list.end(), inner.begin(), inner.end());
    }

    strategy_runner::SweepConfig cfg;
    cfg.strategy_name = STRAT_NAME;
    cfg.out_filename = out_filename;
    cfg.min_trades = MIN_NUMBER_OF_TRADES;
    cfg.min_dd = MIN_ALLOWED_MAX_DRAWBACK;
    cfg.min_reasonable_gain = 50.0f;
    cfg.print_every = 1000;
    cfg.shuffle = false; // ordering above is deliberate; see the comment on the grouping

    strategy_runner::sweep(cfg, std::move(param_list), [&](const BigWill_params &p) {
        return PROCESS(PAIRS, p.AO_fast, p.AO_slow, p.ema_f, p.ema_s, p.max_open_trades);
    });

    strategy_runner::print_timing_and_ram(t_begin, nb_tested);
    TA_Shutdown();
    return 0;
}
