#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include "custom_talib_wrapper.hh"
#include "tools.hh"

namespace
{
float tail_sum(const std::vector<float> &values, const size_t count)
{
    if (values.empty())
    {
        return 0.0f;
    }

    const size_t begin = values.size() > count ? values.size() - count : 0;
    return std::accumulate(values.begin() + begin, values.end(), 0.0f);
}
} // namespace

int main()
{
    std::cout << std::fixed << std::setprecision(6);

    const KLINEf btc_spot = read_input_data("./data/data/binance/1h/BTC-USDT.csv");
    const std::vector<float> ema_fast = TALIB_EMA(btc_spot.close, 11);
    const std::vector<float> ema_slow = TALIB_EMA(btc_spot.close, 448);
    std::vector<float> bb_upper{};
    std::vector<float> bb_middle{};
    std::vector<float> bb_lower{};
    TALIB_BBANDS(btc_spot.close, 2.0f, 2.0f, 20, bb_upper, bb_middle, bb_lower);

    std::cout << "SPOT_ROWS " << btc_spot.nb << "\n";
    std::cout << "EMA11_TAIL_SUM " << tail_sum(ema_fast, 5) << "\n";
    std::cout << "EMA448_TAIL_SUM " << tail_sum(ema_slow, 5) << "\n";
    std::cout << "BBANDS_TAIL_SUM " << tail_sum(bb_upper, 5) << " " << tail_sum(bb_middle, 5) << " " << tail_sum(bb_lower, 5) << "\n";

    std::vector<KLINEf> aligned_pairs{};
    aligned_pairs.push_back(read_input_data("./data/data/binance/1h/BTC-USDT.csv"));
    aligned_pairs.push_back(read_input_data("./data/data/binance/1h/ETH-USDT.csv"));
    const std::vector<uint> start_indexes = INITIALIZE_DATA(aligned_pairs);

    std::cout << "ALIGN_STARTS " << start_indexes[0] << " " << start_indexes[1] << "\n";
    std::cout << "ALIGN_LAST_TS " << aligned_pairs[0].timestamp.back() << " " << aligned_pairs[1].timestamp.back() << "\n";
    std::cout << "ALIGN_LAST_CLOSE_SUM " << aligned_pairs[0].close.back() + aligned_pairs[1].close.back() << "\n";

    const KLINEf btc_futures = read_input_data_f("./data/data/futures/BTC_USDT-1h-futures.json", "2023-06-18");
    const fundings btc_funding = read_funding_rates_data("./data/data/futures/BTC_USDT-8h-funding_rate.json");

    std::cout << "FUTURES_ROWS " << btc_futures.nb << "\n";
    std::cout << "FUNDING_ROWS " << btc_funding.nb << "\n";
    std::cout << "FUNDING_FIRST " << get_funding_fee_if_any(btc_funding, btc_funding.timestamp.front()) << "\n";
    std::cout << "FUNDING_SECOND " << get_funding_fee_if_any(btc_funding, btc_funding.timestamp[1]) << "\n";
    std::cout << "FUNDING_MISS " << get_funding_fee_if_any(btc_funding, btc_funding.timestamp.front() + 3600) << "\n";

    return 0;
}
