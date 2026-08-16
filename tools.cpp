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

int get_hour_from_timestamp(const int timestamp)
{
    time_t rawtime = timestamp;
    struct tm ts;
    char buf[80];

    // Format time, "ddd yyyy-mm-dd hh:mm:ss zzz" ->  "%a %Y-%m-%d %H:%M:%S %Z"
    ts = *localtime(&rawtime);
    strftime(buf, sizeof(buf), "%H", &ts);

    return std::stoi(buf);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int get_year_from_timestamp(const int timestamp)
{
    time_t rawtime = timestamp;
    struct tm ts;
    char buf[80];

    // Format time, "ddd yyyy-mm-dd hh:mm:ss zzz"
    ts = *localtime(&rawtime);
    strftime(buf, sizeof(buf), "%Y", &ts);

    return std::stoi(buf);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int get_month_from_timestamp(const int &timestamp)
{
    const std::time_t timestamp_t = static_cast<std::time_t>(timestamp);
    // Convert the timestamp to a struct tm
    const std::tm *timeinfo = std::localtime(&timestamp_t);
    // Extract the month number
    const int month = timeinfo->tm_mon + 1; // Adding 1 since tm_mon ranges from 0 to 11
    return month;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int get_day_from_timestamp(const int timestamp)
{
    time_t rawtime = timestamp;
    struct tm ts;
    char buf[80];

    // Format time, "ddd yyyy-mm-dd hh:mm:ss zzz"
    ts = *localtime(&rawtime);
    strftime(buf, sizeof(buf), "%d", &ts);

    return std::stoi(buf);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double get_wall_time()
{
    std::chrono::high_resolution_clock m_clock;
    double time = std::chrono::duration_cast<std::chrono::seconds>(m_clock.now().time_since_epoch()).count();
    return time;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double process_mem_usage()
{
    double vm_usage = 0.0;
#if defined(__linux__)

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >>
            ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> vsize >> rss;
    }
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0 / 1024.0;
#endif
    return vm_usage;
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

float calculate_calmar_ratio(const std::vector<int> &times, const std::vector<float> &wallet_vals, const float &max_DD)
{

    if (times.size() <= 4)
        return -100.0;

    const uint last_idx = times.size() - 1;
    const int first_year = get_year_from_timestamp(times[0]);
    const int first_month = get_month_from_timestamp(times[0]);
    const int first_day = get_day_from_timestamp(times[0]);
    const int last_year = get_year_from_timestamp(times[last_idx]);
    const int last_month = get_month_from_timestamp(times[last_idx]);
    const int last_day = get_day_from_timestamp(times[last_idx]);

    const float factor_first_year = (365.0f - float(first_month) * 30.0f - float(first_day)) / 365.0f;
    const float factor_last_year = (float(last_month) * 30.0f + float(last_day)) / 365.0f;

    std::vector<float> vals_begin_years{};
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

    std::vector<float> yearly_pc_changes{};
    yearly_pc_changes.reserve(10);
    for (uint iy = 1; iy < vals_begin_years.size(); iy++)
    {
        yearly_pc_changes.push_back((vals_begin_years[iy] - vals_begin_years[iy - 1]) / vals_begin_years[iy - 1] * 100.0f);
    }

    yearly_pc_changes[0] = yearly_pc_changes[0] * factor_first_year;
    yearly_pc_changes[yearly_pc_changes.size() - 1] = yearly_pc_changes[yearly_pc_changes.size() - 1] * factor_last_year;

    return find_average(yearly_pc_changes) / max_DD;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float calculate_calmar_ratio_monthly(const std::vector<int> &times, const std::vector<float> &wallet_vals, const float &max_DD)
{

    if (times.size() <= 4)
        return -100.0;

    const uint last_idx = times.size() - 1;
    const int first_year = get_year_from_timestamp(times[0]);
    const int first_month = get_month_from_timestamp(times[0]);
    const int first_day = get_day_from_timestamp(times[0]);
    const int last_year = get_year_from_timestamp(times[last_idx]);
    const int last_month = get_month_from_timestamp(times[last_idx]);
    const int last_day = get_day_from_timestamp(times[last_idx]);

    const float factor_first_month = (12.0f - float(first_day) / 30.0f) / 12.0f;
    const float factor_last_month = (float(last_day) / 30.0f) / 12.0f;

    std::vector<float> vals_begin_months{};
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

    std::vector<float> monthly_pc_changes{};
    monthly_pc_changes.reserve(10);
    for (uint iy = 1; iy < vals_begin_months.size(); iy++)
    {
        monthly_pc_changes.push_back((vals_begin_months[iy] - vals_begin_months[iy - 1]) / vals_begin_months[iy - 1] * 100.0f);
    }

    monthly_pc_changes[0] = monthly_pc_changes[0] * factor_first_month;
    monthly_pc_changes[monthly_pc_changes.size() - 1] = monthly_pc_changes[monthly_pc_changes.size() - 1] * factor_last_month;

    return find_average(monthly_pc_changes) / max_DD * 12.0;
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
    }

    klines2 = output_KLINE;
    std::cout << " Done." << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<uint> INITIALIZE_DATA(std::vector<KLINEf> &PAIRS)
// will modify PAIRS since it is passed as reference
{
    std::cout << "Running INITIALIZE_DATA..." << std::endl;

    const int NB_PAIRS = PAIRS.size();

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

std::vector<std::vector<float>> splitVector(const std::vector<float> &inputVector, const int splitSize)
{
    std::vector<std::vector<float>> result;
    const int numSubVectors = inputVector.size() / splitSize;
    const int remainingElements = inputVector.size() % splitSize;
    int currentIndex = 0;

    // std::cout << numSubVectors << std::endl;
    // std::cout << remainingElements << std::endl;

    for (int i = 0; i < numSubVectors; ++i)
    {
        std::vector<float> subVector(inputVector.begin() + currentIndex, inputVector.begin() + currentIndex + splitSize);
        result.push_back(subVector);
        currentIndex += splitSize;
    }

    if (remainingElements > 0)
    {
        std::vector<float> subVector(inputVector.begin() + currentIndex, inputVector.begin() + currentIndex + remainingElements);
        result.push_back(subVector);
    }

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> duplicateElements(const std::vector<float> &inputVector, const int n)
{
    std::vector<float> outputVector;
    outputVector.reserve(inputVector.size() * n);
    for (const auto &element : inputVector)
    {
        for (int i = 0; i < n; ++i)
        {
            outputVector.push_back(element);
        }
    }
    return outputVector;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool isTimestampAtFundingTimes(const int64_t &timestamp)
{
    // Check if the timestamp falls within the desired times
    tm *target_time = gmtime(&timestamp);

    return (target_time->tm_hour == 0 && target_time->tm_min == 0 && target_time->tm_sec == 0) ||
           (target_time->tm_hour == 8 && target_time->tm_min == 0 && target_time->tm_sec == 0) ||
           (target_time->tm_hour == 16 && target_time->tm_min == 0 && target_time->tm_sec == 0);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float get_funding_fee_if_any(const fundings &FUND, const int &current_timestamp)
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

std::vector<int> generateRange_int(const int &vmin, const int &vmax, const int &N)
{
    std::vector<int> result;
    const int step = (vmax - vmin) / (N - 1);

    for (int i = 0; i < N; i++)
    {
        const int value = vmin + step * i;
        result.push_back(value);
    }

    return result;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

    time_t unixTimestamp = std::mktime(&timeStruct);
    return unixTimestamp;
}

#include <iostream>
#include <chrono>
#include <ctime>
#include <sstream>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string getCurrentDateMinusTwoDays()
{
    // Get the current time
    const std::chrono::system_clock::time_point now = std::chrono::system_clock::now();

    // Subtract two days from the current time
    const std::chrono::hours twoDays(48);
    const std::chrono::system_clock::time_point twoDaysAgo = now - twoDays;

    // Convert the time point to a time_t
    const std::time_t time = std::chrono::system_clock::to_time_t(twoDaysAgo);

    // Convert the time_t to a struct tm
    const std::tm *timeStruct = std::localtime(&time);

    // Format the date as a string
    std::stringstream ss;
    ss << std::put_time(timeStruct, "%Y-%m-%d");
    const std::string date = ss.str();

    return date;
}
