#pragma once
#include "config.hh"
#include "Klinef.hh"
#include <nlohmann/json.hpp>
struct FundingEvent
{
    int64_t timestamp = 0;
    double rate = 0, mark_price = 0;
    unsigned millisecond = 0; // retain exchange settlement jitter after a bar opens
};
struct MarketData
{
    std::vector<KLINEf> signal, execution;
    std::vector<KLINEf> higher;
    std::vector<size_t> higher_offset;
    int higher_ratio = 1;
    std::vector<std::vector<FundingEvent>> funding;
    int signal_seconds = 0, execution_seconds = 0;
    int64_t start = 0, end = 0;
    nlohmann::json provenance;
};
KLINEf read_input_data(const std::string &path);
KLINEf read_input_data_f(const std::string &path, const std::string &max_time = "2099-06-20");
fundings read_funding_rates_data(const std::string &path);
void validate_kline(const KLINEf &k, int step = 0);
void align_data(std::vector<KLINEf> &pairs, int step, int64_t begin, int64_t end);
std::vector<uint> INITIALIZE_DATA(std::vector<KLINEf> &pairs);
bool check_timestamp_consistencies(const std::vector<KLINEf> &pairs);
void realign_timestamps(const KLINEf &reference, KLINEf &pair);
MarketData load_market(const backtest_config::StrategyConfig &cfg, size_t warmup_bars);
