#pragma once

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "custom_talib_wrapper.hh"
#include "tools.hh"
#include <ta-lib/ta_libc.h>

// Shared plumbing for parameter-sweep strategies. Strategy files now only own:
//   - their signal logic (the PROCESS lambda passed to sweep())
//   - their indicator precomputation (CALCULATE_INDICATORS)
//   - their parameter ranges
// Everything else — TA-Lib init, banner, result printing, the best-tracking loop,
// the periodic output cadence, the final report — lives here.
//
// This header stays header-only: sweep() is a template and everything else is a
// short free function. No extra Makefile target needed.
namespace strategy_runner
{
struct SweepConfig
{
    std::string strategy_name;
    std::string out_filename;
    uint min_trades = 100;
    float min_dd = -40.0f;            // keep if res.max_DD > min_dd (closer to zero wins)
    float min_reasonable_gain = 0.0f; // reject tiny or negative gains
    float max_reasonable_gain = 1.0e6f;
    uint print_every = 1000;

    // Shuffle the parameter list before sweeping, so an interrupted run has still
    // sampled the whole space rather than a prefix of it.
    //
    // Set false when the caller has already ordered the list deliberately -- typically
    // to group parameters that share an expensive indicator, so a strategy can cache
    // that indicator for one group at a time instead of accumulating every
    // combination. BigWill does exactly this for its Awesome Oscillator.
    bool shuffle = true;
};

// First bar index at which every listed indicator holds a real value rather than warmup
// padding, and at which the pair's own history has begun.
//
// Indicator series are left-padded with zeros through their warmup (see indicators.hh),
// so a loop that starts too early compares prices against zeros -- which silently
// produces trades rather than an error. Strategies traditionally hard-coded a constant
// (600) that happened to exceed the warmups in use; this derives it instead, so raising
// an indicator period cannot quietly outgrow it.
//
//   const uint ii_begin = strategy_runner::first_tradable_index({warm_ema, warm_atr},
//                                                               start_indexes[0],
//                                                               1); // reads [ii - 1]
//
// `extra_lookback` covers strategies that read backwards from the current bar.
inline uint first_tradable_index(const std::initializer_list<size_t> indicator_warmups,
                                 const uint data_start_index = 0,
                                 const uint extra_lookback = 0)
{
    size_t worst = 0;
    for (const size_t w : indicator_warmups)
    {
        worst = std::max(worst, w);
    }
    return std::max(static_cast<uint>(worst), data_start_index) + extra_lookback;
}

inline void init_talib()
{
    const TA_RetCode ret = TA_Initialize();
    if (ret != TA_SUCCESS)
    {
        std::cout << "Cannot initialize TA-Lib !\n"
                  << ret << "\n";
    }
    else
    {
        std::cout << "Initialized TA-Lib !\n";
    }
}

inline std::vector<std::string> build_spot_datafile_paths(const std::vector<std::string> &coins, const std::string &timeframe)
{
    std::vector<std::string> out;
    out.reserve(coins.size());
    for (const std::string &coin : coins)
    {
        out.push_back("./data/data/binance/" + timeframe + "/" + coin + "-USDT.csv");
    }
    return out;
}

inline void print_date_range(const KLINEf &reference)
{
    if (reference.nb == 0)
    {
        return;
    }
    const uint last_idx = reference.nb - 1;
    const int y0 = get_year_from_timestamp(reference.timestamp[0]);
    const int m0 = get_month_from_timestamp(reference.timestamp[0]);
    const int d0 = get_day_from_timestamp(reference.timestamp[0]);
    const int y1 = get_year_from_timestamp(reference.timestamp[last_idx]);
    const int m1 = get_month_from_timestamp(reference.timestamp[last_idx]);
    const int d1 = get_day_from_timestamp(reference.timestamp[last_idx]);
    const std::time_t difference = std::abs(reference.timestamp[last_idx] - reference.timestamp[0]);
    const int days = difference / (24 * 60 * 60);
    std::cout << "Begin day      : " << y0 << "/" << m0 << "/" << d0 << std::endl;
    std::cout << "End day        : " << y1 << "/" << m1 << "/" << d1 << std::endl;
    std::cout << "Number of days : " << YELLOW << days << RESET << std::endl;
}

inline void print_banner(const std::string &name, const std::vector<std::string> &datafiles,
                         const KLINEf &reference, float fee, uint min_trades, float min_dd)
{
    std::cout << "\n--------------------------------------------------------------------------" << std::endl;
    std::cout << "Strategy to test: " << BLUE << name << RESET << std::endl;
    std::cout << "DATA FILES TO PROCESS: " << std::endl;
    for (const std::string &f : datafiles)
    {
        std::cout << "  " << YELLOW << f << RESET << std::endl;
    }
    print_date_range(reference);
    std::cout << "OPEN/CLOSE FEE : " << fee << " %" << std::endl;
    std::cout << "Minimum number of trades required    : " << min_trades << std::endl;
    std::cout << "Maximum drawback (=drawdown) allowed : " << min_dd << " %" << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
}

inline void print_best_res(const std::string &name, const RUN_RESULTf &r)
{
    std::cout << "\n--------------------------------------------------------------------------" << std::endl;
    std::cout << "BEST PARAMETER SET FOUND: " << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
    std::cout << "Time             : " << GREY << GET_CURRENT_TIME_STR() << RESET << std::endl;
    std::cout << "Strategy         : " << BLUE << name << RESET << std::endl;
    std::cout << "Parameters       : " << YELLOW << r.param_str << RESET << std::endl;
    std::cout << "Max Open Trades  : " << YELLOW << r.max_open_trades << RESET << std::endl;
    std::cout << "Gain             : " << r.gain_pc << "%" << std::endl;
    std::cout << "Porfolio         : " << r.WALLET_VAL_USDT << "$ (started with 1000$)" << std::endl;
    std::cout << "Win rate         : " << r.win_rate << "%" << std::endl;
    std::cout << "max DD           : " << r.max_DD << "%" << std::endl;
    std::cout << "Gain/DDC         : " << r.gain_over_DDC << std::endl;
    std::cout << "Score            : " << GREEN << r.score << RESET << std::endl;
    std::cout << "Calmar ratio     : " << r.calmar_ratio << std::endl;
    if (r.calmar_ratio_monthly != 0.0f)
    {
        std::cout << "Calmar monthly   : " << r.calmar_ratio_monthly << std::endl;
    }
    std::cout << "Number of trades : " << r.nb_posi_entered << std::endl;
    std::cout << "Total fees paid  : " << std::round(r.total_fees_paid * 100.0f) / 100.0f << "$ (started with 1000$)" << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
}

inline void print_timing_and_ram(double t_begin, uint nb_tested)
{
    const double t_end = get_wall_time();
    std::cout << "Number of backtests performed : " << nb_tested << std::endl;
    std::cout << "Time taken                    : " << t_end - t_begin << " seconds " << std::endl;
    const double ram_usage = process_mem_usage();
    std::cout << "RAM usage                     : " << std::round(ram_usage * 10.0) / 10.0 << " MB" << std::endl;
    std::cout << "--------------------------------------------------------------------------" << std::endl;
}

// Runs a parameter sweep. `process_fn` is called once per parameter; `print_fn`
// customizes result-printing (default matches the canonical multi-pair format).
//
// `params` is taken by value because it is shuffled in place, so callers should
// std::move() their list in -- BigWill's is around 5.7 million entries and copying it
// cost over 100 MB of pointless allocation and memcpy per run.
//
// Single-threaded by design, for now. process_fn typically populates each pair's
// IndicatorCache lazily on first sight of a parameter combination (see BigWill), so
// parallelising this loop needs that cache made thread-safe first. The lazy caching is
// worth far more than the core count: it removes a redundancy factor in the thousands,
// where threading would give at most the number of cores.
template <typename ParamsT, typename ProcessFn,
          typename PrintFn = std::function<void(const RUN_RESULTf &)>>
RUN_RESULTf sweep(const SweepConfig &cfg, std::vector<ParamsT> params, ProcessFn process_fn,
                  PrintFn print_fn = nullptr)
{
    auto do_print = [&](const RUN_RESULTf &r) {
        if (print_fn)
        {
            print_fn(r);
        }
        else
        {
            print_best_res(cfg.strategy_name, r);
        }
    };

    if (cfg.shuffle)
    {
        random_shuffle_vector(params);
    }

    RUN_RESULTf best{};
    best.score = -std::numeric_limits<float>::infinity();
    best.gain_over_DDC = -100.0f;
    best.calmar_ratio = -100.0f;

    uint nb_done = 0;
    uint tick = 0;
    for (const ParamsT &p : params)
    {
        const RUN_RESULTf res = process_fn(p);
        ++nb_done;
        ++tick;

        if (res.score > best.score && res.gain_pc > cfg.min_reasonable_gain &&
            res.gain_pc < cfg.max_reasonable_gain &&
            res.nb_posi_entered >= int(cfg.min_trades) && res.max_DD > cfg.min_dd)
        {
            best = res;
        }

        if (tick >= cfg.print_every)
        {
            tick = 0;
            do_print(best);
            WRITE_OR_UPDATE_BEST_SCORE_FILE(cfg.strategy_name, cfg.out_filename, best);
            const float pc_done = std::round(float(nb_done) / float(params.size()) * 100.0f * 100.0f) / 100.0f;
            std::cout << "DONE " << nb_done << " / " << params.size() << "   = " << pc_done << "%" << std::endl;
        }
    }

    do_print(best);
    WRITE_OR_UPDATE_BEST_SCORE_FILE(cfg.strategy_name, cfg.out_filename, best);
    return best;
}
} // namespace strategy_runner
