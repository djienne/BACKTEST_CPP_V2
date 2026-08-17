#include "tools.hh"

#include <cassert>
#include <numeric>

using json = nlohmann::json;

namespace
{
bool parse_kline_csv_row(const std::string &line, long int &ts, float &op, float &hi, float &lo, float &cl, float &vol)
{
    return std::sscanf(line.c_str(), "%ld,%f,%f,%f,%f,%f", &ts, &op, &hi, &lo, &cl, &vol) == 6;
}
} // namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float find_average(const std::vector<float> &vec)
{
    if (vec.empty())
    {
        return 0.0f;
    }
    return std::accumulate(vec.begin(), vec.end(), 0.0f) / static_cast<float>(vec.size());
}

double find_average(const std::vector<double> &vec)
{
    if (vec.empty())
    {
        return 0.0;
    }
    return std::accumulate(vec.begin(), vec.end(), 0.0) / static_cast<double>(vec.size());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float find_min(const std::vector<float> &vec)
{
    if (vec.empty())
    {
        return std::numeric_limits<float>::max();
    }
    return *std::min_element(vec.begin(), vec.end());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int find_max(const std::vector<int> &vec)
{
    if (vec.empty())
    {
        return std::numeric_limits<int>::lowest();
    }
    return *std::max_element(vec.begin(), vec.end());
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// All calendar fields are derived in UTC via the reentrant gmtime_r.
//
// Two bugs are fixed here at once. localtime()/gmtime() return a pointer to a shared
// static tm, so the multithreaded F_* strategies raced on it through
// calculate_calmar_ratio. And localtime() made every result depend on the host's
// timezone: the same data and parameters produced different year/month/day boundaries
// -- and so different Calmar ratios and month-rollover snapshots -- on two machines.
// Exchange kline timestamps are UTC, so UTC is the correct reading of them.
//
// Reading the tm fields directly also replaces a strftime()-then-std::stoi() round-trip
// per call, which these helpers do once or twice per candle in the hot loops.
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int get_year_from_timestamp(const int64_t timestamp)
{
    return utc_tm_from_timestamp(timestamp).tm_year + 1900; // tm_year counts from 1900
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int get_month_from_timestamp(const int64_t timestamp)
{
    return utc_tm_from_timestamp(timestamp).tm_mon + 1; // tm_mon ranges 0..11
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int get_day_from_timestamp(const int64_t timestamp)
{
    return utc_tm_from_timestamp(timestamp).tm_mday;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double get_wall_time()
{
    std::chrono::high_resolution_clock m_clock;
    double time = std::chrono::duration_cast<std::chrono::seconds>(m_clock.now().time_since_epoch()).count();
    return time;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Resident set size in MB -- the memory actually held in RAM.
//
// This used to read the same /proc/self/stat line but return `vsize`, the *virtual*
// size, while every caller printed it as "RAM usage". Virtual size counts reservations
// the process never touches and typically overstates the real footprint several times
// over. The rss field and the page size were both read and then ignored.
double process_mem_usage()
{
#if defined(__linux__)
    unsigned long vsize = 0;
    long rss_pages = 0;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >>
            ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> vsize >> rss_pages;
    }
    // rss is reported in pages; the page size is read rather than assumed, since x86-64
    // can be configured with 2 MB pages.
    const long page_size_bytes = sysconf(_SC_PAGE_SIZE);
    return double(rss_pages) * double(page_size_bytes) / 1024.0 / 1024.0;
#else
    return 0.0;
#endif
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Half-open: [min, max) with custom step. Intentionally different endpoint than the
// two-arg overload below; changing this would shift every existing parameter sweep.
std::vector<int> integer_range(const int min, const int max, const int step)
{
    assert(min <= max);
    assert(step > 0);
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
    assert(min <= max);
    std::vector<int> the_range;
    the_range.reserve(static_cast<size_t>(max - min + 1));
    for (int i = min; i <= max; i++)
    {
        the_range.push_back(i);
    }
    return the_range;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> float_Nvalues_range(const float &vmin, const float &vmax, const int &N)
{
    assert(N >= 1);
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string ReplaceAll(std::string str, const std::string &from, const std::string &to)
{
    size_t start_pos = 0;
    while ((start_pos = str.find(from, start_pos)) != std::string::npos)
    {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double calculate_calmar_ratio(const std::vector<int64_t> &times, const std::vector<double> &wallet_vals, const double drawdown_normalizer)
{

    if (times.size() <= 4)
        return -100.0;

    const uint last_idx = times.size() - 1;
    const int first_month = get_month_from_timestamp(times[0]);
    const int first_day = get_day_from_timestamp(times[0]);
    const int last_month = get_month_from_timestamp(times[last_idx]);
    const int last_day = get_day_from_timestamp(times[last_idx]);

    // The first and last calendar years are usually partial, so their returns are
    // pro-rated by the fraction of the year actually covered.
    const double factor_first_year = (365.0 - double(first_month) * 30.0 - double(first_day)) / 365.0;
    const double factor_last_year = (double(last_month) * 30.0 + double(last_day)) / 365.0;

    std::vector<double> vals_begin_years{};
    vals_begin_years.reserve(10);

    vals_begin_years.push_back(1000.0);

    for (uint ii = 1; ii < times.size(); ii++)
    {
        const int year_b = get_year_from_timestamp(times[ii - 1]);
        const int year = get_year_from_timestamp(times[ii]);

        if (year_b != year || ii == times.size() - 1)
        {
            vals_begin_years.push_back(wallet_vals[ii]);
        }
    }

    std::vector<double> yearly_pc_changes{};
    yearly_pc_changes.reserve(10);
    for (uint iy = 1; iy < vals_begin_years.size(); iy++)
    {
        yearly_pc_changes.push_back((vals_begin_years[iy] - vals_begin_years[iy - 1]) / vals_begin_years[iy - 1] * 100.0);
    }

    if (yearly_pc_changes.empty() || !(drawdown_normalizer > 0.0))
    {
        return -100.0;
    }

    yearly_pc_changes[0] = yearly_pc_changes[0] * factor_first_year;
    yearly_pc_changes[yearly_pc_changes.size() - 1] = yearly_pc_changes[yearly_pc_changes.size() - 1] * factor_last_year;

    return find_average(yearly_pc_changes) / drawdown_normalizer;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double calculate_calmar_ratio_monthly(const std::vector<int64_t> &times, const std::vector<double> &wallet_vals, const double drawdown_normalizer)
{

    if (times.size() <= 4)
        return -100.0;

    const uint last_idx = times.size() - 1;
    const int first_day = get_day_from_timestamp(times[0]);
    const int last_day = get_day_from_timestamp(times[last_idx]);

    const double factor_first_month = (12.0 - double(first_day) / 30.0) / 12.0;
    const double factor_last_month = (double(last_day) / 30.0) / 12.0;

    std::vector<double> vals_begin_months{};
    vals_begin_months.reserve(50);

    vals_begin_months.push_back(1000.0);

    for (uint ii = 1; ii < times.size(); ii++)
    {
        const int month_b = get_month_from_timestamp(times[ii - 1]);
        const int month = get_month_from_timestamp(times[ii]);

        if (month_b != month || ii == times.size() - 1)
        {
            vals_begin_months.push_back(wallet_vals[ii]);
        }
    }

    std::vector<double> monthly_pc_changes{};
    monthly_pc_changes.reserve(10);
    for (uint iy = 1; iy < vals_begin_months.size(); iy++)
    {
        monthly_pc_changes.push_back((vals_begin_months[iy] - vals_begin_months[iy - 1]) / vals_begin_months[iy - 1] * 100.0);
    }

    if (monthly_pc_changes.empty() || !(drawdown_normalizer > 0.0))
    {
        return -100.0;
    }

    monthly_pc_changes[0] = monthly_pc_changes[0] * factor_first_month;
    monthly_pc_changes[monthly_pc_changes.size() - 1] = monthly_pc_changes[monthly_pc_changes.size() - 1] * factor_last_month;

    return find_average(monthly_pc_changes) / drawdown_normalizer * 12.0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string GET_CURRENT_TIME_STR()
{
    std::time_t currentTime = std::time(nullptr);
    std::string currentTimeString;
    char buffer[80];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", std::localtime(&currentTime));
    currentTimeString = buffer;
    return currentTimeString;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool f_exists(const std::string &name)
{
    std::ifstream f(name.c_str());
    return f.good();
}

double read_score(const std::string &out_filename)
{
    std::ifstream file(out_filename, std::fstream::in);

    if (!file.is_open())
    {
        BACKTEST_FATAL("Failed to open score file: " + out_filename);
    }

    std::regex scoreRegex(R"(Score\s*:\s*([+-]?\d*\.?\d+))");
    std::smatch match;
    std::string line;
    double score = 0.0;

    while (std::getline(file, line))
    {
        if (std::regex_search(line, match, scoreRegex))
        {
            score = std::stod(match[1]);
            break;
        }
    }

    file.close();

    return score;
}

void WRITE_OR_UPDATE_BEST_SCORE_FILE(const std::string &STRAT_NAME, const std::string &out_filename, const RUN_RESULTf &result)
{
    // Read-modify-write on a shared file. The F_* strategies run this from several
    // worker threads at once and the multithreaded sweep runner does too, so without
    // this lock two writers could interleave a truncating open with the other's read
    // and lose the best result -- or leave the file half-written.
    static std::mutex score_file_mutex;
    const std::lock_guard<std::mutex> lock(score_file_mutex);

    std::cout << "Reading or updating SCORE FILE..." << std::endl;
    float file_score = 0.0;
    bool file_exists = false;

    // read previous best score in file

    if (f_exists(out_filename))
    {
        file_exists = true;
        std::cout << "File " << out_filename << " exists." << std::endl;
        file_score = read_score(out_filename);
        if (result.score <= file_score)
        {
            std::cout << "Read score in file: " << GREEN << file_score << RESET << " (the old, inside the file, is better or the same)" << std::endl;
        }
        else
        {
            std::cout << "Read score in file: " << GREEN << file_score << RESET << " (this is a new high score)" << std::endl;
        }
        std::cout << "--------------------------------------------------------------------------" << std::endl;
    }
    else
    {
        file_exists = false;
        std::cout << "File " << out_filename << " does not exist." << std::endl;
    }

    // if this score is better than the one in file, overwrite the file

    if (result.score > file_score || !file_exists)
    {
        std::ofstream myfile(out_filename, std::ios::trunc);
        if (myfile.is_open())
        {
            myfile << "\n-------------------------------------" << std::endl;
            myfile << "BEST PARAMETER SET FOUND: " << std::endl;
            myfile << "-------------------------------------" << std::endl;
            myfile << "Time             : " << GET_CURRENT_TIME_STR() << std::endl;
            myfile << "Strategy         : " << STRAT_NAME << std::endl;
            myfile << "Parameters       : " << result.param_str << std::endl;
            myfile << "Max Open Trades  : " << result.max_open_trades << std::endl;
            myfile << "Gain             : " << result.gain_pc << "%" << std::endl;
            myfile << "Porfolio         : " << result.WALLET_VAL_USDT << "$ (started with 1000$)" << std::endl;
            myfile << "Win rate         : " << result.win_rate << "%" << std::endl;
            myfile << "max DD           : " << result.max_DD << "%" << std::endl;
            myfile << "Gain/DDC         : " << result.gain_over_DDC << std::endl;
            myfile << "Score            : " << result.score << std::endl;
            myfile << "Calmar ratio     : " << result.calmar_ratio << std::endl;
            myfile << "Number of trades : " << result.nb_posi_entered << std::endl;
            myfile << "Total fees paid  : " << round(result.total_fees_paid * 100.0f) / 100.0f << "$ (started with 1000$)" << std::endl;
            myfile << "-------------------------------------" << std::endl;
            myfile.close();
        }
        else
        {
            std::cout << "Unable to open file" << std::endl;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool check_if_necessary(const KLINEf &klines_btc, const KLINEf &klines2)
{
    bool result = false;

    for (uint i = klines2.start_idx + 1; i < klines_btc.nb; i++)
    {
        if (klines_btc.timestamp[i] != klines2.timestamp[i])
        {
            result = true;
            std::cout << i << std::endl;
            break;
        }
    }

    return result;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool check_timestamp_consistencies(const std::vector<KLINEf> &PAIRS)
{

    bool is_ok = true;
    // index 0 must be BTC (with largest timestamp list with non-zero values of open, high,...)
    for (size_t i = 1; i < PAIRS.size(); i++)
    {
        std::cout << "Checking timestamp consistencies of " << PAIRS[i].name << std::endl;

        const bool val = check_if_necessary(PAIRS[0], PAIRS[i]);

        if (val)
        {
            is_ok = false;

            std::cout << "Timestamps are not consistent" << std::endl;
            const size_t start = PAIRS[0].timestamp.size() > 100 ? PAIRS[0].timestamp.size() - 100 : 0;
            for (size_t j = start; j < PAIRS[0].timestamp.size(); ++j)
            {
                std::cout << PAIRS[0].timestamp[j] << " " << PAIRS[i].timestamp[j] << " " << PAIRS[0].timestamp[j] - PAIRS[i].timestamp[j] << std::endl;
            }

            BACKTEST_FATAL("Inconsistent timestamps between " + PAIRS[0].name + " and " + PAIRS[i].name);
        }
    }

    if (is_ok)
    {
        std::cout << "OK." << std::endl;
    }

    return is_ok;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void realign_timestamps(const KLINEf &klines_btc, KLINEf &klines2) // klines2 will be modified
{

    std::cout << "Re-aligning timestamps of " << klines2.name << " if necessary..." << std::endl;

    if (!check_if_necessary(klines_btc, klines2))
    {
        std::cout << "It is not necessary." << std::endl;
        return;
    }

    std::cout << "Re-aligningment is necessary. Doing it..." << std::endl;

    KLINEf output_KLINE = klines_btc;
    output_KLINE.name = klines2.name;
    output_KLINE.start_idx = klines2.start_idx;

    std::unordered_map<long int, uint> timestamp_to_index{};
    timestamp_to_index.reserve(klines2.nb);
    for (uint j = 0; j < klines2.nb; ++j)
    {
        timestamp_to_index.emplace(klines2.timestamp[j], j);
    }

    for (uint i = klines2.start_idx; i < klines_btc.nb; i++)
    {
        const auto match = timestamp_to_index.find(klines_btc.timestamp[i]);
        if (match == timestamp_to_index.end())
        {
            BACKTEST_FATAL("realign_timestamps: no match for timestamp in " + klines2.name);
        }

        const uint j = match->second;
        output_KLINE.open[i] = klines2.open[j];
        output_KLINE.high[i] = klines2.high[j];
        output_KLINE.low[i] = klines2.low[j];
        output_KLINE.close[i] = klines2.close[j];
        output_KLINE.volume[i] = klines2.volume[j];
    }

    klines2 = output_KLINE;
    std::cout << " Done." << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<uint> INITIALIZE_DATA(std::vector<KLINEf> &PAIRS)
// will modify PAIRS since it is passed as reference
{
    std::cout << "Running INITIALIZE_DATA..." << std::endl;

    const uint NB_PAIRS = static_cast<uint>(PAIRS.size());

    std::vector<uint> start_indexes{};
    start_indexes.reserve(NB_PAIRS);

    start_indexes.push_back(600);

    const long int timeframe = (PAIRS[0].timestamp[11] - PAIRS[0].timestamp[10]) / 60000;

    if (timeframe <= 5) // skip first ~10 minutes for the 5 min TF to avoid bug
    {
        for (uint ic = 1; ic < NB_PAIRS; ic++)
        {
            PAIRS[ic].timestamp.erase(PAIRS[ic].timestamp.begin(), PAIRS[ic].timestamp.begin() + 3000);
            PAIRS[ic].open.erase(PAIRS[ic].open.begin(), PAIRS[ic].open.begin() + 3000);
            PAIRS[ic].high.erase(PAIRS[ic].high.begin(), PAIRS[ic].high.begin() + 3000);
            PAIRS[ic].low.erase(PAIRS[ic].low.begin(), PAIRS[ic].low.begin() + 3000);
            PAIRS[ic].close.erase(PAIRS[ic].close.begin(), PAIRS[ic].close.begin() + 3000);
            PAIRS[ic].volume.erase(PAIRS[ic].volume.begin(), PAIRS[ic].volume.begin() + 3000);
            PAIRS[ic].nb = PAIRS[ic].close.size();
        }
    }

    // Find each pair's first index within the reference (index 0) series. Pairs list
    // later than BTC, so each one starts at a different offset.
    //
    // The scan stops at the first hit. It previously ran to the end of the reference
    // series and pushed an entry per match, so a repeated timestamp would append a
    // second start index and shift every later pair's entry -- silently desynchronising
    // start_indexes[ic] from PAIRS[ic].
    for (uint ic = 1; ic < NB_PAIRS; ic++)
    {
        bool found = false;
        for (uint i = 0; i < PAIRS[0].timestamp.size(); i++)
        {
            if (PAIRS[ic].timestamp[0] == PAIRS[0].timestamp[i])
            {
                start_indexes.push_back(start_indexes[0] + i);
                std::cout << "Start for " + PAIRS[ic].name << " : " << start_indexes[0] + i << std::endl;
                PAIRS[ic].start_idx = start_indexes[0] + i;
                found = true;
                break;
            }
        }

        if (!found)
        {
            BACKTEST_FATAL("INITIALIZE_DATA: could not find start index for pair " + PAIRS[ic].name);
        }
    }

    for (uint ic = 1; ic < NB_PAIRS; ic++)
    {
        while (PAIRS[ic].close.size() != PAIRS[0].close.size())
        {
            const int nb_to_add = PAIRS[0].close.size() - PAIRS[ic].close.size();
            PAIRS[ic].timestamp = add_zeros(PAIRS[ic].timestamp, nb_to_add);
            PAIRS[ic].open = add_zeros(PAIRS[ic].open, nb_to_add);
            PAIRS[ic].high = add_zeros(PAIRS[ic].high, nb_to_add);
            PAIRS[ic].low = add_zeros(PAIRS[ic].low, nb_to_add);
            PAIRS[ic].close = add_zeros(PAIRS[ic].close, nb_to_add);
            PAIRS[ic].volume = add_zeros(PAIRS[ic].volume, nb_to_add);
        }

        PAIRS[ic].nb = PAIRS[0].close.size();

        for (uint ii = 0; ii < start_indexes[ic] + 5; ii++)
        {
            PAIRS[ic].timestamp[ii] = PAIRS[0].timestamp[ii];
        }
    }

    for (uint ic = 1; ic < NB_PAIRS; ic++)
    {
        realign_timestamps(PAIRS[0], PAIRS[ic]);
    }

    for (uint ic = 1; ic < NB_PAIRS; ic++)
    {
        if (PAIRS[ic].nb != PAIRS[0].open.size() || PAIRS[ic].open.size() != PAIRS[0].open.size())
        {
            BACKTEST_FATAL("INITIALIZE_DATA: inconsistent size for pair " + PAIRS[ic].name);
        }
    }

    if (!check_timestamp_consistencies(PAIRS))
    {
        BACKTEST_FATAL("Last timestamps are not consistent after realignment");
    }

    std::cout << "Initialized calculations." << std::endl;

    return start_indexes;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

KLINEf read_input_data(const std::string &input_file_path)
{
    KLINEf kline;

    std::ifstream myfile(input_file_path);
    if (!myfile.is_open())
    {
        BACKTEST_FATAL("Failed to open CSV data file: " + input_file_path);
    }

    std::string line;
    long int ts;
    float op, hi, lo, cl, vol;
    long int previous_ts = 0;
    int nb_read = 0;
    bool skipped_first = false;

    while (std::getline(myfile, line))
    {
        if (!skipped_first)
        {
            skipped_first = true;
            continue;
        }

        if (line.empty())
        {
            continue;
        }

        if (!parse_kline_csv_row(line, ts, op, hi, lo, cl, vol))
        {
            BACKTEST_FATAL("Malformed CSV row in " + input_file_path + " at data row " + std::to_string(nb_read + 1));
        }

        if (previous_ts == ts)
        {
            std::cout << "FOUND DUPLICATE TS at index " << nb_read << "; IGNORING. " << std::endl;
            previous_ts = ts;
            continue;
        }
        previous_ts = ts;
        kline.timestamp.push_back(ts / long(1000));
        kline.open.push_back(op);
        kline.high.push_back(hi);
        kline.low.push_back(lo);
        kline.close.push_back(cl);
        kline.volume.push_back(vol);
        nb_read++;
    }

    myfile.close();

    kline.nb = int(kline.close.size());

    std::filesystem::path filePath(input_file_path);
    std::string filenameWithExtension = filePath.filename().string();
    std::string filenameWithoutExtension = filenameWithExtension.substr(0, filenameWithExtension.find_last_of('.'));
    kline.name = filenameWithoutExtension;

    std::cout << "Loaded data file: " << YELLOW << input_file_path << RESET << std::endl;

    return kline;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

KLINEf read_input_data_f(const std::string &input_file_path, const std::string &max_time)
{
    const time_t unixTimestamp = convertToUnixTimestamp(max_time);

    KLINEf kline{};

    long int ts;
    float op, hi, lo, cl, vol;
    long int previous_ts = 0;
    int nb_read = 0;

    std::ifstream jsonFile(input_file_path);
    if (!jsonFile.is_open())
    {
        BACKTEST_FATAL("Failed to open JSON data file: " + input_file_path);
    }
    const json jsonData = json::parse(jsonFile);
    kline.timestamp.reserve(jsonData.size());
    kline.open.reserve(jsonData.size());
    kline.high.reserve(jsonData.size());
    kline.low.reserve(jsonData.size());
    kline.close.reserve(jsonData.size());
    kline.volume.reserve(jsonData.size());

    // Iterate over the elements
    for (const auto &item : jsonData)
    {
        // Iterate over the values within each item
        ts = item[0];
        op = item[1];
        hi = item[2];
        lo = item[3];
        cl = item[4];
        vol = item[5];

        if (previous_ts == ts)
        {
            std::cout << "FOUND DUPLICATE TS at index " << nb_read << "; IGNORING. " << std::endl;
            previous_ts = ts;
            continue;
        }
        previous_ts = ts;
        if (ts / long(1000) > unixTimestamp)
        {
            continue;
        }
        kline.timestamp.push_back(ts / long(1000));
        kline.open.push_back(op);
        kline.high.push_back(hi);
        kline.low.push_back(lo);
        kline.close.push_back(cl);
        kline.volume.push_back(vol);
        nb_read++;
    }

    kline.nb = int(kline.close.size());

    std::filesystem::path filePath(input_file_path);
    std::string filenameWithExtension = filePath.filename().string();
    std::string filenameWithoutExtension = filenameWithExtension.substr(0, filenameWithExtension.find_last_of('.'));
    kline.name = filenameWithoutExtension;

    std::cout << "Loaded data file: " << YELLOW << input_file_path << RESET << std::endl;

    // std::cout << kline.timestamp.back() << " " << kline.open.back() << " " << kline.high.back() << " " << kline.low.back() << " " << kline.close.back() << std::endl;
    // std::cout << nb_read << " " << kline.nb << std::endl;

    return kline;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

fundings read_funding_rates_data(const std::string &input_file_path)
{
    fundings FR{};

    long int ts;
    float fu;
    long int previous_ts = 0;
    int nb_read = 0;

    std::ifstream jsonFile(input_file_path);
    if (!jsonFile.is_open())
    {
        BACKTEST_FATAL("Failed to open JSON data file: " + input_file_path);
    }
    const json jsonData = json::parse(jsonFile);
    FR.timestamp.reserve(jsonData.size());
    FR.funding.reserve(jsonData.size());
    FR.funding_by_timestamp.reserve(jsonData.size());

    // Iterate over the elements
    for (const auto &item : jsonData)
    {
        // Iterate over the values within each item
        ts = item[0];
        fu = item[1];

        if (previous_ts == ts)
        {
            std::cout << "FOUND DUPLICATE TS at index " << nb_read << "; IGNORING. " << std::endl;
            previous_ts = ts;
            continue;
        }
        previous_ts = ts;
        const long int timestamp_seconds = ts / long(1000);
        FR.timestamp.push_back(timestamp_seconds);
        FR.funding.push_back(fu);
        FR.funding_by_timestamp.emplace(timestamp_seconds, fu);
        nb_read++;
    }

    FR.nb = int(FR.funding.size());

    std::filesystem::path filePath(input_file_path);
    std::string filenameWithExtension = filePath.filename().string();
    std::string filenameWithoutExtension = filenameWithExtension.substr(0, filenameWithExtension.find_last_of('.'));
    FR.name = filenameWithoutExtension;

    std::cout << "Loaded data file: " << YELLOW << input_file_path << RESET << std::endl;

    return FR;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// splitVector and duplicateElements lived here. They existed only to hand-roll
// timeframe resampling inside the two mtf strategies; RESAMPLE_TIMEFRAME and
// PROJECT_HTF_TO_LTF in indicators.* replaced both, with boundary alignment and
// anti-lookahead handling those helpers did not have.

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Perpetual funding settles at 00:00, 08:00 and 16:00 UTC. gmtime_r rather than gmtime:
// the latter returns a pointer into a shared static tm, which the multithreaded F_*
// strategies raced on -- this is called once per candle per pair.
bool isTimestampAtFundingTimes(const int64_t timestamp)
{
    const std::tm target_time = utc_tm_from_timestamp(timestamp);

    if (target_time.tm_min != 0 || target_time.tm_sec != 0)
    {
        return false;
    }

    return target_time.tm_hour == 0 || target_time.tm_hour == 8 || target_time.tm_hour == 16;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float get_funding_fee_if_any(const fundings &FUND, const int64_t current_timestamp)
{
    if (!isTimestampAtFundingTimes(current_timestamp))
    {
        return 0.0f;
    }

    const auto match = FUND.funding_by_timestamp.find(current_timestamp);
    if (match != FUND.funding_by_timestamp.end())
    {
        return match->second;
    }

    return 0.0f;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fill_datafile_paths_f(const std::vector<std::string> &COINS, const std::string &timeframe, std::vector<std::string> &DATAFILES, std::vector<std::string> &DATAFILES_fundings)
{
    for (uint i = 0; i < COINS.size(); i++)
    {
        DATAFILES.push_back("./data/data/futures/" + COINS[i] + "_USDT-" + timeframe + "-futures.json");
        DATAFILES_fundings.push_back("./data/data/futures/" + COINS[i] + "_USDT-8h-funding_rate.json");
        // BTC_USDT-5m-futures.json
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// N values spread evenly across the closed interval [vmin, vmax], endpoints included.
//
// The step used to be computed as the integer division (vmax - vmin) / (N - 1), which
// truncated and so never reached vmax. generateRange_int(3, 600, 300) produced step 1
// and therefore 3..302 -- half the intended search space, silently. N == 1 also divided
// by zero. Interpolating in floating point and rounding each point fixes both.
//
// Duplicates are removed, so the result may be shorter than N when the interval is
// narrower than N (e.g. 5 points requested across [1, 3]).
std::vector<int> generateRange_int(const int &vmin, const int &vmax, const int &N)
{
    assert(N >= 1);
    assert(vmin <= vmax);

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
        const int value = static_cast<int>(std::llround(vmin + t * (static_cast<double>(vmax) - static_cast<double>(vmin))));
        if (result.empty() || result.back() != value)
        {
            result.push_back(value);
        }
    }

    return result;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parses "YYYY-MM-DD" as midnight **UTC**.
//
// This used std::mktime, which interprets the tm as local time. read_input_data_f uses
// the result as a cutoff, so the number of rows loaded from a futures file depended on
// the machine's timezone -- the same data and the same max_time gave a different row
// count on a UTC+2 laptop than on a UTC CI runner, which is how this surfaced. timegm is
// the UTC counterpart and matches the UTC convention the rest of the codebase uses.
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
