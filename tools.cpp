#include "tools.hh"

#include <cassert>
#include <numeric>

using json = nlohmann::json;

float find_min(const std::vector<float> &vec)
{
    if (vec.empty())
    {
        return std::numeric_limits<float>::max();
    }
    return *std::min_element(vec.begin(), vec.end());
}

// Note: returns lowest() for an empty vector and so correctly handles all-negative inputs (prior
// versions initialized the running max with 0 and silently masked negative-only data).
float find_max(const std::vector<float> &vec)
{
    if (vec.empty())
    {
        return std::numeric_limits<float>::lowest();
    }
    return *std::max_element(vec.begin(), vec.end());
}

// Reentrant UTC calendar conversion.
namespace
{
std::tm utc_tm_from_timestamp(const int64_t timestamp)
{
    const std::time_t raw = static_cast<std::time_t>(timestamp);
    std::tm out{};
    gmtime_r(&raw, &out);
    return out;
}
} // namespace

int get_hour_from_timestamp(const int64_t timestamp)
{
    return utc_tm_from_timestamp(timestamp).tm_hour;
}

int get_year_from_timestamp(const int64_t timestamp)
{
    return utc_tm_from_timestamp(timestamp).tm_year + 1900; // tm_year counts from 1900
}

int get_month_from_timestamp(const int64_t timestamp)
{
    return utc_tm_from_timestamp(timestamp).tm_mon + 1; // tm_mon ranges 0..11
}

int get_day_from_timestamp(const int64_t timestamp)
{
    return utc_tm_from_timestamp(timestamp).tm_mday;
}

double get_wall_time()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// Current resident set size in MiB on Linux.
double process_mem_usage()
{
#if defined(__linux__)
    unsigned long vsize = 0;
    long rss_pages = 0;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >>
            ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >>
            ignore >> ignore >> vsize >> rss_pages;
    }
    // /proc reports resident memory in pages.
    const long page_size_bytes = sysconf(_SC_PAGE_SIZE);
    return double(rss_pages) * double(page_size_bytes) / 1024.0 / 1024.0;
#else
    return 0.0;
#endif
}

// Half-open: [min, max) with custom step. Intentionally different endpoint than the
// two-arg overload below; changing this would shift every existing parameter sweep.
std::vector<int> integer_range(const int min, const int max, const int step)
{
    if (min > max)
        throw std::invalid_argument("Invalid range bounds");
    if (step <= 0)
        throw std::invalid_argument("Range step must be positive");
    std::vector<int> the_range;
    the_range.reserve(static_cast<size_t>((max - min + step - 1) / step));
    for (int i = min; i < max; i += step)
    {
        the_range.push_back(i);
    }
    return the_range;
}

// Closed: [min, max] with step 1. Kept distinct from the three-arg overload above for
// backwards compatibility with existing sweeps.
std::vector<int> integer_range(const int min, const int max)
{
    if (min > max)
        throw std::invalid_argument("Invalid range bounds");
    std::vector<int> the_range;
    the_range.reserve(static_cast<size_t>(max - min + 1));
    for (int i = min; i <= max; i++)
    {
        the_range.push_back(i);
    }
    return the_range;
}

std::vector<float> float_Nvalues_range(const float &vmin, const float &vmax, const int &N)
{
    if (N < 1)
        throw std::invalid_argument("Range needs at least one value");
    std::vector<float> result;
    result.reserve(static_cast<size_t>(N));
    if (N == 1)
    {
        result.push_back(vmin);
        return result;
    }

    const float step = (vmax - vmin) / (float(N) - 1.0f);
    for (int i = 0; i < N; i++)
    {
        result.push_back(vmin + step * i);
    }
    return result;
}

float get_funding_fee_if_any(const fundings &f, const int64_t timestamp)
{
    const auto it = f.funding_by_timestamp.find(timestamp);
    return it == f.funding_by_timestamp.end() ? 0.0f : it->second;
}

// Rounded evenly spaced values, with duplicates removed.
std::vector<int> generateRange_int(const int &vmin, const int &vmax, const int &N)
{
    if (N < 1)
        throw std::invalid_argument("Range needs at least one value");
    if (vmin > vmax)
        throw std::invalid_argument("Invalid range bounds");

    std::vector<int> result;
    result.reserve(static_cast<size_t>(N));

    if (N == 1)
    {
        result.push_back(vmin);
        return result;
    }

    for (int i = 0; i < N; i++)
    {
        const double t = static_cast<double>(i) / static_cast<double>(N - 1);
        const int value =
            static_cast<int>(std::llround(vmin + t * (static_cast<double>(vmax) - static_cast<double>(vmin))));
        if (result.empty() || result.back() != value)
        {
            result.push_back(value);
        }
    }

    return result;
}
// Parse YYYY-MM-DD at midnight UTC; -1 means invalid input.
time_t convertToUnixTimestamp(const std::string &dateString)
{
    std::tm timeStruct = {};
    std::istringstream iss(dateString);
    iss >> std::get_time(&timeStruct, "%Y-%m-%d");

    if (iss.fail())
    {
        // Failed to parse the date string
        return -1;
    }

    timeStruct.tm_hour = 0;
    timeStruct.tm_min = 0;
    timeStruct.tm_sec = 0;

    return timegm(&timeStruct);
}

#include <iostream>
#include <chrono>
#include <ctime>
#include <sstream>

// "YYYY-MM-DD" for two days ago, in UTC. Pairs with convertToUnixTimestamp, so both ends
// of a max_time cutoff agree; localtime() here would have reintroduced the same
// timezone dependence, and is not reentrant.
std::string getCurrentDateMinusTwoDays()
{
    const std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    const std::chrono::hours twoDays(48);
    const std::time_t time = std::chrono::system_clock::to_time_t(now - twoDays);

    std::tm timeStruct{};
    gmtime_r(&time, &timeStruct);

    std::stringstream ss;
    ss << std::put_time(&timeStruct, "%Y-%m-%d");
    return ss.str();
}
