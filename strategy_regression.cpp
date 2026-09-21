// Calls the same EMA evaluator as the production binary, without duplicating signals.
#define main ema_program_main
#include "backtest_double_EMA_float.cpp"
#undef main
int main(int argc, char **argv)
{
    try
    {
        if (argc != 4)
            throw std::runtime_error("Usage: strategy_regression candles.csv slow fast");
        MarketData d;
        d.signal = {read_input_data(argv[1])};
        d.execution = d.signal;
        d.funding.resize(1);
        d.signal_seconds = d.execution_seconds = static_cast<int>(d.signal[0].timestamp[1] - d.signal[0].timestamp[0]);
        validate_kline(d.signal[0], d.signal_seconds);
        Params p{double(std::stoi(argv[2])), double(std::stoi(argv[3]))};
        const size_t begin = static_cast<size_t>(std::max(p[0], p[1])) + 1;
        std::vector<IndicatorCache> cache(1);
        init_talib();
        auto result = evaluate_double_ema(d, cache, p, {begin, d.execution[0].nb}, true);
        std::cout << result_json(result).dump() << "\n";
        TA_Shutdown();
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
