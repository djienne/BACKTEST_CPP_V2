#include "data_io.hh"
#include <iostream>
#include <numeric>
int main(int argc, char **argv)
{
    try
    {
        if (argc != 2)
            throw std::runtime_error("Usage: verification_regression.exe candles.csv");
        const auto k = read_input_data(argv[1]);
        if (k.nb < 2)
            throw std::runtime_error("Need at least two candles");
        validate_kline(k, static_cast<int>(k.timestamp[1] - k.timestamp[0]));
        std::cout << nlohmann::json({{"rows", k.nb},
                                     {"close_sum", std::accumulate(k.close.begin(), k.close.end(), 0.0)},
                                     {"first", k.timestamp.front()},
                                     {"last", k.timestamp.back()}})
                         .dump()
                  << "\n";
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << "\n";
        return 1;
    }
}
