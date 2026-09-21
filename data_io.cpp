#include "data_io.hh"
#include "tools.hh"
#include "indicators.hh"
#include <algorithm>
#include <cerrno>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <spawn.h>
#include <sstream>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>
extern char **environ;
namespace
{
void append(KLINEf &k, int64_t ms, float o, float h, float l, float c, float v)
{
    if (ms % 1000)
        throw std::runtime_error("Candle timestamp is not in whole seconds");
    if (!k.timestamp.empty() && ms / 1000 == k.timestamp.back())
    {
        if (o == k.open.back() && h == k.high.back() && l == k.low.back() && c == k.close.back() &&
            v == k.volume.back())
            return;
        throw std::runtime_error("Conflicting duplicate candle in " + k.name);
    }
    k.timestamp.push_back(ms / 1000);
    k.open.push_back(o);
    k.high.push_back(h);
    k.low.push_back(l);
    k.close.push_back(c);
    k.volume.push_back(v);
}
template <class T> void slice(std::vector<T> &v, size_t first, size_t last)
{
    v = std::vector<T>(v.begin() + first, v.begin() + last);
}
void crop(KLINEf &k, int64_t begin, int64_t end)
{
    size_t a = std::lower_bound(k.timestamp.begin(), k.timestamp.end(), begin) - k.timestamp.begin();
    size_t b = std::lower_bound(k.timestamp.begin(), k.timestamp.end(), end) - k.timestamp.begin();
    slice(k.timestamp, a, b);
    slice(k.open, a, b);
    slice(k.high, a, b);
    slice(k.low, a, b);
    slice(k.close, a, b);
    slice(k.volume, a, b);
    k.nb = static_cast<uint>(k.close.size());
    k.start_idx = 0;
    k.indicators.clear();
}
} // namespace
void validate_kline(const KLINEf &k, int step)
{
    const size_t n = k.close.size();
    if (!n || k.timestamp.size() != n || k.open.size() != n || k.high.size() != n || k.low.size() != n ||
        k.volume.size() != n)
        throw std::runtime_error("Empty/inconsistent OHLCV arrays: " + k.name);
    for (size_t i = 0; i < n; ++i)
    {
        if (!(std::isfinite(k.open[i]) && std::isfinite(k.high[i]) && std::isfinite(k.low[i]) &&
              std::isfinite(k.close[i]) && std::isfinite(k.volume[i])) ||
            std::min({k.open[i], k.high[i], k.low[i], k.close[i]}) <= 0 || k.volume[i] < 0 ||
            k.low[i] > std::min(k.open[i], k.close[i]) || k.high[i] < std::max(k.open[i], k.close[i]) ||
            (i && k.timestamp[i] <= k.timestamp[i - 1]))
            throw std::runtime_error("Invalid OHLCV row in " + k.name + " at " + std::to_string(i));
        if (step && (k.timestamp[i] % step || (i && k.timestamp[i] - k.timestamp[i - 1] != step)))
            throw std::runtime_error("Unrepaired candle gap in " + k.name + " at " + std::to_string(k.timestamp[i]));
    }
}
KLINEf read_input_data(const std::string &path)
{
    KLINEf k{};
    k.name = std::filesystem::path(path).stem().string();
    std::ifstream f(path);
    if (!f)
        throw std::runtime_error("Cannot open " + path);
    std::string line;
    std::getline(f, line);
    if (!line.empty() && line.back() == '\r')
        line.pop_back();
    if (line != "date,open,high,low,close,volume")
        throw std::runtime_error("Invalid CSV header: " + path);
    while (std::getline(f, line))
    {
        if (line.empty())
            continue;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream row(line);
        int64_t t;
        float o, h, l, c, v;
        std::string extra;
        if (!(row >> t >> o >> h >> l >> c >> v) || row >> extra)
            throw std::runtime_error("Invalid CSV row: " + path);
        append(k, t, o, h, l, c, v);
    }
    k.nb = static_cast<uint>(k.close.size());
    validate_kline(k);
    return k;
}
KLINEf read_input_data_f(const std::string &path, const std::string &max_time)
{
    KLINEf k{};
    k.name = std::filesystem::path(path).stem().string();
    std::ifstream f(path);
    if (!f)
        throw std::runtime_error("Cannot open " + path);
    const auto rows = nlohmann::json::parse(f);
    const auto limit = convertToUnixTimestamp(max_time);
    for (const auto &r : rows)
    {
        if (!r.is_array() || r.size() != 6)
            throw std::runtime_error("Invalid futures row");
        const auto t = r[0].get<int64_t>();
        if (t / 1000 <= limit)
            append(k, t, r[1], r[2], r[3], r[4], r[5]);
    }
    k.nb = static_cast<uint>(k.close.size());
    validate_kline(k);
    return k;
}
fundings read_funding_rates_data(const std::string &path)
{
    // Legacy rate-only fixtures are readable for loader tests, not research settlement.
    std::ifstream f(path);
    if (!f)
        throw std::runtime_error("Cannot open " + path);
    auto rows = nlohmann::json::parse(f);
    fundings out{};
    out.name = path;
    for (const auto &r : rows)
    {
        const auto t = r.at(0).get<int64_t>() / 1000;
        const auto rate = r.at(1).get<float>();
        if (!std::isfinite(rate) || (!out.timestamp.empty() && t <= out.timestamp.back()))
            throw std::runtime_error("Invalid funding history");
        out.timestamp.push_back(t);
        out.funding.push_back(rate);
        out.funding_by_timestamp.emplace(t, rate);
    }
    out.nb = static_cast<uint>(out.timestamp.size());
    return out;
}
void align_data(std::vector<KLINEf> &pairs, int step, int64_t begin, int64_t end)
{
    if (pairs.empty() || begin >= end)
        throw std::runtime_error("Empty alignment interval");
    for (auto &k : pairs)
    {
        crop(k, begin, end);
        validate_kline(k, step);
        if (k.timestamp.front() != begin || k.timestamp.back() + step != end)
            throw std::runtime_error("Incomplete alignment interval: " + k.name);
    }
}
bool check_timestamp_consistencies(const std::vector<KLINEf> &pairs)
{
    if (pairs.empty())
        return false;
    for (const auto &k : pairs)
        if (k.timestamp != pairs.front().timestamp)
            return false;
    return true;
}
void realign_timestamps(const KLINEf &reference, KLINEf &pair)
{
    validate_kline(reference);
    validate_kline(pair);
    if (reference.nb < 2)
        throw std::runtime_error("Need two reference timestamps");
    const int step = static_cast<int>(reference.timestamp[1] - reference.timestamp[0]);
    crop(pair, reference.timestamp.front(), reference.timestamp.back() + step);
    validate_kline(pair, step);
    if (pair.timestamp != reference.timestamp)
        throw std::runtime_error("Cannot align missing timestamps");
}
std::vector<uint> INITIALIZE_DATA(std::vector<KLINEf> &pairs)
{
    if (pairs.empty())
        throw std::runtime_error("No pairs");
    for (const auto &p : pairs)
        validate_kline(p);
    if (pairs[0].nb < 2)
        throw std::runtime_error("Need two timestamps");
    const int step = static_cast<int>(pairs[0].timestamp[1] - pairs[0].timestamp[0]);
    int64_t begin = pairs[0].timestamp.front(), end = pairs[0].timestamp.back() + step;
    for (const auto &p : pairs)
    {
        begin = std::max(begin, p.timestamp.front());
        end = std::min(end, p.timestamp.back() + step);
    }
    align_data(pairs, step, begin, end);
    return std::vector<uint>(pairs.size(), 0);
}
MarketData load_market(const backtest_config::StrategyConfig &cfg, size_t warmup)
{
    const auto temp = std::filesystem::temp_directory_path() / ("backtest-data-" + std::to_string(getpid()) + ".json");
    std::vector<std::string> args = {
        "python3",       "tools/download_data.py", "--config",   cfg.path,     "--ensure", cfg.name,
        "--warmup-bars", std::to_string(warmup),   "--manifest", temp.string()};
    if (cfg.offline)
        args.push_back("--offline");
    std::vector<char *> argv;
    for (auto &a : args)
        argv.push_back(a.data());
    argv.push_back(nullptr);
    pid_t pid;
    const int error = posix_spawnp(&pid, "python3", nullptr, nullptr, argv.data(), environ);
    if (error)
        throw std::runtime_error("Cannot launch data repair; use the documented Docker environment");
    int status = 0;
    pid_t waited;
    do
    {
        waited = waitpid(pid, &status, 0);
    } while (waited < 0 && errno == EINTR);
    if (waited < 0 || !WIFEXITED(status) || WEXITSTATUS(status))
    {
        std::filesystem::remove(temp);
        throw std::runtime_error("Data preparation failed; see missing intervals above");
    }
    nlohmann::json j;
    try
    {
        std::ifstream f(temp);
        j = nlohmann::json::parse(f);
    }
    catch (...)
    {
        std::filesystem::remove(temp);
        throw;
    }
    std::filesystem::remove(temp);
    MarketData d;
    d.provenance = j;
    d.start = j.at("start_ms").get<int64_t>() / 1000;
    d.end = j.at("end_ms").get<int64_t>() / 1000;
    d.signal_seconds = j.at("signal_seconds");
    d.execution_seconds = j.at("execution_seconds");
    for (const auto &path : j.at("signal_files"))
        d.signal.push_back(cfg.is_futures() ? read_input_data_f(path) : read_input_data(path));
    align_data(d.signal, d.signal_seconds, j.at("history_start_ms").get<int64_t>() / 1000, d.end);
    if (!cfg.htf.empty())
    {
        d.higher_ratio = backtest_config::timeframe_seconds(cfg.htf) / d.signal_seconds;
        for (const auto &k : d.signal)
        {
            auto r = RESAMPLE_TIMEFRAME(k, d.higher_ratio, d.signal_seconds / 60,
                                        backtest_config::timeframe_seconds(cfg.htf) / 60);
            d.higher.push_back(std::move(r.kline));
            d.higher_offset.push_back(r.ltf_offset);
        }
    }
    for (const auto &path : j.at("execution_files"))
        d.execution.push_back(cfg.is_futures() ? read_input_data_f(path) : read_input_data(path));
    align_data(d.execution, d.execution_seconds, d.start, d.end);
    d.funding.resize(d.signal.size());
    for (size_t p = 0; p < j.at("funding_files").size(); ++p)
    {
        std::ifstream f(j["funding_files"][p].get<std::string>());
        const auto ledger = nlohmann::json::parse(f);
        if (ledger.at("start_ms").get<int64_t>() / 1000 > d.start || ledger.at("end_ms").get<int64_t>() / 1000 < d.end)
            throw std::runtime_error("Funding coverage changed during load");
        for (const auto &e : ledger.at("events"))
        {
            const auto ms = e.at("timestamp_ms").get<int64_t>();
            FundingEvent event{ms / 1000, e.at("rate"), e.at("mark_price"), static_cast<unsigned>(ms % 1000)};
            if (event.timestamp >= d.start && event.timestamp < d.end)
                d.funding[p].push_back(event);
        }
    }
    return d;
}
