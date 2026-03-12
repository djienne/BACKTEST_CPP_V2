#include <iomanip>
#include <iostream>
#include <unordered_map>

#include "custom_talib_wrapper.hh"
#include "trade_core.hh"
#include "tools.hh"

using namespace std;
using uint = unsigned int;

namespace
{
static const uint NB_PAIRS = 2;
const vector<string> DATAFILES = {
    "./data/data/binance/1h/BTC-USDT.csv",
    "./data/data/binance/1h/ETH-USDT.csv",
};
constexpr float FEE = 0.1f;
constexpr float USDT_amount_initial = 1000.0f;
constexpr float STOCH_RSI_UPPER = 0.800f;
constexpr float STOCH_RSI_LOWER = 0.200f;
constexpr int EMA_SHORT = 89;
constexpr int EMA_LONG = 279;

trade_core::ResultMetrics run_multi_pair_case()
{
    vector<KLINEf> pairs{};
    pairs.reserve(NB_PAIRS);
    for (const string &data_file : DATAFILES)
    {
        pairs.push_back(read_input_data(data_file));
    }

    const vector<uint> start_indexes = INITIALIZE_DATA(pairs);

    array<unordered_map<string, vector<float>>, NB_PAIRS> ema_lists{};
    array<vector<float>, NB_PAIRS> stoch_rsi{};
    for (uint ic = 0; ic < NB_PAIRS; ++ic)
    {
        ema_lists[ic]["EMA_" + std::to_string(EMA_SHORT)] = TALIB_EMA(pairs[ic].close, EMA_SHORT);
        ema_lists[ic]["EMA_" + std::to_string(EMA_LONG)] = TALIB_EMA(pairs[ic].close, EMA_LONG);
        stoch_rsi[ic] = TALIB_STOCHRSI_not_averaged(pairs[ic].close, 14, 14);
    }

    trade_core::WalletTrace wallet_trace{};
    trade_core::TradeStats stats{};
    trade_core::PortfolioState<NB_PAIRS> portfolio(USDT_amount_initial);

    const uint nb_max = pairs[0].nb;
    const uint ii_begin = start_indexes[0];

    for (uint ii = ii_begin; ii < nb_max; ++ii)
    {
        const bool last_iteration = (ii == nb_max - 1);
        bool closed = false;

        for (uint ic = 0; ic < NB_PAIRS; ++ic)
        {
            if (ii < start_indexes[ic])
            {
                continue;
            }

            const bool open_long = ema_lists[ic]["EMA_" + std::to_string(EMA_SHORT)][ii] >= ema_lists[ic]["EMA_" + std::to_string(EMA_LONG)][ii] &&
                                   stoch_rsi[ic][ii] < STOCH_RSI_UPPER;
            const bool close_long = ema_lists[ic]["EMA_" + std::to_string(EMA_SHORT)][ii] <= ema_lists[ic]["EMA_" + std::to_string(EMA_LONG)][ii] &&
                                    stoch_rsi[ic][ii] > STOCH_RSI_LOWER;

            if (portfolio.coin_amounts[ic] > 0.0f && (close_long || last_iteration))
            {
                trade_core::close_spot_long(portfolio, stats, ic, pairs[ic].close[ii], FEE);
                closed = true;
            }

            if (portfolio.coin_amounts[ic] == 0.0f && open_long && !last_iteration && portfolio.active_positions < NB_PAIRS)
            {
                trade_core::open_spot_long(portfolio, stats, ic, pairs[ic].close[ii], FEE, NB_PAIRS);
            }
        }

        if (closed || last_iteration)
        {
            array<float, NB_PAIRS> closes{};
            for (uint ic = 0; ic < NB_PAIRS; ++ic)
            {
                closes[ic] = pairs[ic].close[ii];
            }

            trade_core::record_spot_snapshot(portfolio, wallet_trace, closes, pairs[0].timestamp[ii]);
        }
    }

    array<float, NB_PAIRS> last_closes{};
    for (uint ic = 0; ic < NB_PAIRS; ++ic)
    {
        last_closes[ic] = pairs[ic].close[nb_max - 1];
    }

    const float wallet_val_usdt = trade_core::calculate_spot_wallet_val_usdt(portfolio, last_closes);
    const trade_core::ResultMetrics metrics = trade_core::calculate_result_metrics(wallet_val_usdt, USDT_amount_initial, portfolio.max_drawdown, stats);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "MULTI_SPOT_EMA_SHORT " << EMA_SHORT << "\n";
    std::cout << "MULTI_SPOT_EMA_LONG " << EMA_LONG << "\n";
    std::cout << "MULTI_SPOT_GAIN " << metrics.gain << "\n";
    std::cout << "MULTI_SPOT_WIN_RATE " << metrics.win_rate << "\n";
    std::cout << "MULTI_SPOT_MAX_DD " << portfolio.max_drawdown << "\n";
    std::cout << "MULTI_SPOT_GAIN_OVER_DDC " << (metrics.gain / metrics.ddc) << "\n";
    std::cout << "MULTI_SPOT_SCORE " << metrics.score << "\n";
    std::cout << "MULTI_SPOT_TRADES " << stats.nb_positions_entered << "\n";

    return metrics;
}
} // namespace

int main()
{
    TA_RetCode retCode = TA_Initialize();
    if (retCode != TA_SUCCESS)
    {
        std::cout << "Cannot initialize TA-Lib !\n"
                  << retCode << "\n";
        return 1;
    }

    run_multi_pair_case();
    return 0;
}
