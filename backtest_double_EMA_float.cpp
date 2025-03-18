#include <iostream>
#include <ctime>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <unordered_map>
#include "tools.hh"
#include "custom_talib_wrapper.hh"
#include <ta-lib/ta_libc.h>

using namespace std;

const string STRAT_NAME = "2EMA_crossover";
const string out_filename = STRAT_NAME + "_best.txt";

const float FRACTION_PER_POSI = 1.0; // FRACTION OF CAPITAL PER POSITION
const float LEV = 1.0;               // LEVERAGE
const float FEE = 0.1;               // FEES
const float FUNDING_FEE = 0.00;      // FUNDING FEE APPLIED EVERY 8 hours
const bool CAN_LONG = true;          // LONG ON OR OFF
const bool CAN_SHORT = false;        // SHORT ON OR OFF
const double MAX_ALLOWED_DD = -40.0;
const std::string PAIR = "BTC-USDT";
const std::string DATAFILE = "./data/data/binance/1h/" + PAIR + ".csv";

// RANGE OF EMA PERIDOS TO TEST
const int MIN_EMA = 3;
std::vector<int> range1 = integer_range(MIN_EMA, 600);
std::vector<int> range2 = integer_range(MIN_EMA, 600);
//////////////////////////

uint i_print = 0;

std::unordered_map<std::string, std::vector<float>> EMA_LISTS{};
std::vector<float> year{};
std::vector<float> hour{};
std::vector<float> month{};
std::vector<float> day{};

KLINEf kline{};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_best_res(const RUN_RESULTf best)
{
    const int year = get_year_from_timestamp(kline.timestamp[0]);
    const int month = get_month_from_timestamp(kline.timestamp[0]);
    const int day = get_day_from_timestamp(kline.timestamp[0]);
    const uint last_idx = kline.nb - 1;
    const int last_year = get_year_from_timestamp(kline.timestamp[last_idx]);
    const int last_month = get_month_from_timestamp(kline.timestamp[last_idx]);
    const int last_day = get_day_from_timestamp(kline.timestamp[last_idx]);

    std::time_t difference = std::abs(int(kline.timestamp[last_idx]) - int(kline.timestamp[0]));
    const int days = difference / (24 * 60 * 60);

    // Display info
    std::cout << "\n-------------------------------------" << std::endl;
    std::cout << "TIME RANGE: " << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << " Begin day : " << year << "/" << month << "/" << day << endl;
    std::cout << " End day   : " << last_year << "/" << last_month << "/" << last_day << endl;
    std::cout << " Duration  : " << days << " days" << endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "BEST PARAMETER SET FOUND: " << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << " Time             : " << GREY << GET_CURRENT_TIME_STR() << RESET << std::endl;
    std::cout << " Strategy         : " << BLUE << STRAT_NAME << RESET << std::endl;
    std::cout << " PAIR             : " << PAIR << std::endl;
    std::cout << " EMAs             : " << best.ema1 << " " << best.ema2 << std::endl;
    std::cout << " Gain             : " << best.gain_pc << "%" << std::endl;
    std::cout << " Win rate         : " << best.win_rate << "%" << std::endl;
    std::cout << " max DD           : " << best.max_DD << "%" << std::endl;
    std::cout << " Gain over DDC    : " << best.gain_over_DDC << std::endl;
    std::cout << " Score            : " << GREEN << best.score << RESET << std::endl;
    std::cout << " Number of trades : " << best.nb_posi_entered << std::endl;
    std::cout << "-------------------------------------"
              << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void INITIALIZE_DATA0(const KLINEf &kline)
{
    std::cout << "Initializing data (calculating indicators)...\n"
              << std::endl;

    std::vector<int> list_ema = {};

    for (int i = MIN_EMA; i <= range2.size() + 5; i++)
    {
        list_ema.push_back(i);
    }

    for (const int i : list_ema)
    {
        EMA_LISTS["EMA" + std::to_string(i)] = TALIB_EMA(kline.close, i);
    }
    cout << "Done calculating EMAs." << endl;

    for (uint ii = 0; ii < kline.nb; ii++)
    {
        year.push_back(get_year_from_timestamp(kline.timestamp[ii]));
        hour.push_back(get_hour_from_timestamp(kline.timestamp[ii]));
        month.push_back(get_month_from_timestamp(kline.timestamp[ii]));
        day.push_back(get_day_from_timestamp(kline.timestamp[ii]));
    }

    cout << "Initialized calculations." << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

RUN_RESULTf PROCESS(const KLINEf &KLINEf, const int ema1_v, const int ema2_v)
{
    std::vector<float> EMA1 = EMA_LISTS["EMA" + std::to_string(ema1_v)];
    std::vector<float> EMA2 = EMA_LISTS["EMA" + std::to_string(ema2_v)];

    bool LAST_ITERATION = false, OPEN_LONG_CONDI = false, OPEN_SHORT_CONDI = false, CLOSE_LONG_CONDI = false, CLOSE_SHORT_CONDI = false;
    bool IN_POSITION = false, IN_LONG = false, IN_SHORT = false;
    int nb_profit = 0, nb_loss = 0, NB_POSI_ENTERED = 0;
    int nb_max = KLINEf.nb;
    float current_price = 0, starting_price = 0, pc_change = 0, pc_change_with_max = 0, max_drawdown = 0, price_position_open = 0;
    float USDT_amount = 1000.0;
    float MAX_USDT_AMOUNT = USDT_amount;
    float amount_b = 0;

    std::vector<float> close = KLINEf.close;

    int ii_begin = find_max(range1) + 2;

    for (size_t ii = ii_begin; ii < KLINEf.nb; ii++)
    {
        if (ii == nb_max - 1)
        {
            LAST_ITERATION = true;
        }

        current_price = close[ii];

        if (ii == ii_begin)
        {
            starting_price = close[ii];
        }

        pc_change = (USDT_amount - 1000.0) / 1000.0 * 100.0;

        if (USDT_amount > MAX_USDT_AMOUNT)
        {
            MAX_USDT_AMOUNT = USDT_amount;
        }

        pc_change_with_max = (USDT_amount - MAX_USDT_AMOUNT) / MAX_USDT_AMOUNT * 100.0;

        if (pc_change_with_max < max_drawdown)
        {
            max_drawdown = pc_change_with_max;
        }

        double USDT_amount_check = 100.0;

        if (IN_SHORT && IN_POSITION)
        {
            USDT_amount_check = USDT_amount - (current_price - price_position_open) / price_position_open * FRACTION_PER_POSI * USDT_amount * LEV;
        }
        if (IN_LONG && IN_POSITION)
        {
            USDT_amount_check = USDT_amount + (current_price - price_position_open) / price_position_open * FRACTION_PER_POSI * USDT_amount * LEV;
        }

        if (USDT_amount_check <= 0)
        {
            RUN_RESULTf result;

            result.AMOUNT_USDT = 0;
            result.gain_over_DDC = 0;
            result.gain_pc = 0;
            result.max_DD = -100.0;
            result.nb_posi_entered = NB_POSI_ENTERED;
            result.win_rate = 0;
            result.score = 0;
            result.ema1 = ema1_v;
            result.ema2 = ema2_v;

            return result;
        }

        // Funding fees

        if (IN_POSITION)
        {
            if (((hour[ii] >= 2) && (hour[ii - 1] < 2)) || ((hour[ii] >= 10) && (hour[ii - 1] < 10)) || ((hour[ii] >= 18) && (hour[ii - 1] < 18)))
            {
                float to_rm = current_price / price_position_open * FRACTION_PER_POSI * USDT_amount * LEV * FUNDING_FEE / 100.0;
                USDT_amount = USDT_amount - to_rm;
            }
        }

        // check if should go in position

        OPEN_LONG_CONDI = (EMA2[ii] >= EMA1[ii]) && (EMA2[ii - 1] <= EMA1[ii - 1]);
        OPEN_SHORT_CONDI = (EMA2[ii] <= EMA1[ii]) && (EMA2[ii - 1] >= EMA1[ii - 1]);

        CLOSE_LONG_CONDI = (EMA2[ii] <= EMA1[ii]) && (EMA2[ii - 1] >= EMA1[ii - 1]);
        CLOSE_SHORT_CONDI = (EMA2[ii] >= EMA1[ii]) && (EMA2[ii - 1] <= EMA1[ii - 1]);

        // IT IS IMPORTANT TO CHECK FIRST FOR CLOSING POSITION AND THEN FOR OPENING POSITION

        // CLOSE SHORT
        if ((IN_POSITION && IN_SHORT) && (CLOSE_SHORT_CONDI || LAST_ITERATION))
        {
            amount_b = USDT_amount;

            USDT_amount = USDT_amount - (current_price - price_position_open) / price_position_open * FRACTION_PER_POSI * USDT_amount * LEV;

            // apply FEEs
            if (FEE > 0)
            {
                USDT_amount = USDT_amount - current_price / price_position_open * FRACTION_PER_POSI * amount_b * LEV * FEE / 100.0;
            }
            else
            {
                USDT_amount = USDT_amount - current_price / price_position_open * FRACTION_PER_POSI * amount_b * FEE / 100.0;
            }
            //

            IN_POSITION = false;
            IN_SHORT = false;
            IN_LONG = false;

            pc_change = -(current_price - price_position_open) / price_position_open * 100.0;

            if (current_price < price_position_open)
            {
                nb_profit = nb_profit + 1;
            }
            else
            {
                nb_loss = nb_loss + 1;
            }
        }
        // CLOSE LONG
        if ((IN_POSITION && IN_LONG) && (CLOSE_LONG_CONDI || LAST_ITERATION))
        {
            amount_b = USDT_amount;

            USDT_amount = USDT_amount + (current_price - price_position_open) / price_position_open * FRACTION_PER_POSI * USDT_amount * LEV;

            // apply FEEs
            if (FEE > 0)
            {
                USDT_amount = USDT_amount - current_price / price_position_open * FRACTION_PER_POSI * amount_b * LEV * FEE / 100.0;
            }
            else
            {
                USDT_amount = USDT_amount - current_price / price_position_open * FRACTION_PER_POSI * amount_b * FEE / 100.0;
            }
            //

            IN_POSITION = false;
            IN_SHORT = false;
            IN_LONG = false;

            pc_change = (current_price - price_position_open) / price_position_open * 100.0;

            if (current_price > price_position_open)
            {
                nb_profit = nb_profit + 1;
            }
            else
            {
                nb_loss = nb_loss + 1;
            }
        }

        // Check to open position (should always be after check of closing)

        // OPEN SHORT
        if ((IN_POSITION == false) && (OPEN_SHORT_CONDI && CAN_SHORT))
        {

            price_position_open = current_price;

            // apply FEEs
            if (FEE > 0.0)
            {
                USDT_amount = USDT_amount - FRACTION_PER_POSI * USDT_amount * LEV * FEE / 100.0;
            }
            else
            {
                USDT_amount = USDT_amount - FRACTION_PER_POSI * USDT_amount * FEE / 100.0;
            }

            IN_POSITION = true;
            IN_LONG = false;
            IN_SHORT = true;

            NB_POSI_ENTERED = NB_POSI_ENTERED + 1;
        }
        // OPEN LONG
        if ((IN_POSITION == false) && (OPEN_LONG_CONDI && CAN_LONG))
        {

            price_position_open = current_price;

            // apply FEEs
            if (FEE > 0.0)
            {
                USDT_amount = USDT_amount - FRACTION_PER_POSI * USDT_amount * LEV * FEE / 100.0;
            }
            else
            {
                USDT_amount = USDT_amount - FRACTION_PER_POSI * USDT_amount * FEE / 100.0;
            }

            IN_POSITION = true;
            IN_LONG = true;
            IN_SHORT = false;

            NB_POSI_ENTERED = NB_POSI_ENTERED + 1;
        }
    }

    float gain = (USDT_amount - 1000.0) / 1000.0 * 100.0;
    float WR = float(nb_profit) / float(NB_POSI_ENTERED) * 100.0;

    float DDC = (1.0 / (1.0 + max_drawdown / 100.0) - 1.0) * 100.0;

    float score = gain / DDC * WR;

    i_print++;

    RUN_RESULTf result;

    result.AMOUNT_USDT = USDT_amount;
    result.gain_over_DDC = gain / DDC;
    result.gain_pc = gain;
    result.max_DD = max_drawdown;
    result.nb_posi_entered = NB_POSI_ENTERED;
    result.win_rate = WR;
    result.score = score;
    result.ema1 = ema1_v;
    result.ema2 = ema2_v;
    result.param_str = "PAIR: " + PAIR + "; EMA1: " + std::to_string(ema1_v) + "; EMA2: " + std::to_string(ema2_v);
    result.max_open_trades = 1;
    result.total_fees_paid = NAN;

    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

KLINEf read_input_data0(const std::string &input_file_path)
{
    KLINEf kline;

    ifstream myfile(input_file_path);
    std::string line;
    long ts;
    long ts_prev = 0;
    float op, hi, lo, cl, vol;

    bool skipped_first = false;

    std::cout << "Reading data file..." << std::endl;

    if (myfile.is_open())
    {
        while (getline(myfile, line))
        {
            if (!skipped_first)
            {
                skipped_first = true;
                continue;
            }
            sscanf(line.c_str(), "%ld,%f,%f,%f,%f,%f", &ts, &op, &hi, &lo, &cl, &vol);
            if (ts_prev == ts)
            {
                std::cout << "Found duplciate. Skipping." << std::endl;
                continue;
            }
            ts_prev = ts;
            kline.timestamp.push_back(ts / long(1000));
            kline.open.push_back(op);
            kline.high.push_back(hi);
            kline.low.push_back(lo);
            kline.close.push_back(cl);
        }
        myfile.close();

        kline.nb = int(kline.close.size());
    }
    std::cout << "Done." << std::endl;
    return kline;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    auto startTime = std::chrono::high_resolution_clock::now();

    std::cout << "\n-------------------------------------" << std::endl;
    std::cout << "Strategy to test: EMA crossover simple" << std::endl;
    std::cout << "DATA FILE TO PROCESS: " << DATAFILE << std::endl;

    TA_RetCode retCode;
    retCode = TA_Initialize();
    if (retCode != TA_SUCCESS)
    {
        std::cout << "Cannot initialize TA-Lib !\n"
                  << retCode << "\n";
    }
    else
    {
        std::cout << "Initialized TA-Lib !\n";
    }

    kline = read_input_data0(DATAFILE);

    INITIALIZE_DATA0(kline);

    RUN_RESULTf best{};
    best.gain_over_DDC = -100.0;

    int i_print2 = 0;
    int i_print3 = 0;

    std::cout << "Begining day : " << year[0] << "/" << month[0] << "/" << day[0] << std::endl;
    std::cout << "End day      : " << year.back() << "/" << month.back() << "/" << day.back() << std::endl;
    std::cout << "CAN LONG     : " << CAN_LONG << std::endl; 
    std::cout << "CAN SHORT    : " << CAN_SHORT << std::endl;

    std::vector<doubleEMA_params> param_list{};
    param_list.reserve(range1.size() * range2.size());

    for (const int ema1 : range1)
    {
        for (const int ema2 : range2)
        {
            if (std::abs(ema1 - ema2) < 7 || ema1 < ema2)
                continue;

            const doubleEMA_params to_add{ema1, ema2};

            param_list.push_back(to_add);
        }
    }

    random_shuffle_vector(param_list);

    // MAIN LOOP

    uint nb_done = 0;

    for (const auto para : param_list)
    {
        RUN_RESULTf res = PROCESS(kline, para.ema1, para.ema2);
        nb_done++;

        if (res.score > best.score && res.gain_pc < 10000.0 && res.nb_posi_entered >= 100 && res.max_DD > MAX_ALLOWED_DD)
        {
            best = res;
        }

        i_print3++;
        if (i_print3 == 10000)
        {
            print_best_res(best);
            WRITE_OR_UPDATE_BEST_SCORE_FILE(STRAT_NAME, out_filename, best);
            i_print3 = 0;
            const double pc_done = std::round(double(nb_done) / double(param_list.size()) * 100.0 * 100.0) / 100.0;
            std::cout << "DONE " << nb_done << " / " << param_list.size() << "   = " << pc_done << "%" << std::endl;
        }
    }

    print_best_res(best);
    WRITE_OR_UPDATE_BEST_SCORE_FILE(STRAT_NAME, out_filename, best);

    // End the timer
    auto endTime = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration<double>(endTime - startTime).count();

    std::cout << "CPU wall time: " << duration << " seconds" << std::endl;
    std::cout << "-------------------------------------\n"
              << std::endl;
}
