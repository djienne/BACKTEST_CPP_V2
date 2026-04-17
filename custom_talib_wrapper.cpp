#include "custom_talib_wrapper.hh"
using namespace std;

namespace
{
using SingleInputTalibFn = TA_RetCode (*)(int, int, const float *, int, TA_Integer *, TA_Integer *, TA_Real *);
using OhlcInputTalibFn = TA_RetCode (*)(int, int, const float *, const float *, const float *, int, TA_Integer *, TA_Integer *, TA_Real *);

std::vector<float> build_talib_output(const size_t input_size, const TA_Integer outBeg, const TA_Integer outNbElement, const std::vector<TA_Real> &raw_output, const char *indicator_name)
{
    std::vector<float> output(input_size, 0.0f);
    for (TA_Integer ii = 0; ii < outNbElement; ++ii)
    {
        output[static_cast<size_t>(outBeg + ii)] = raw_output[static_cast<size_t>(ii)];
    }

    if (output.size() != input_size)
    {
        std::cout << "error in " << indicator_name << std::endl;
        std::abort();
    }

    return output;
}

std::vector<float> run_single_input_indicator(const std::vector<float> &vals, const int period, SingleInputTalibFn indicator, const char *indicator_name)
{
    if (vals.empty())
    {
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw_output(vals.size(), 0.0f);

    const TA_RetCode retCode = indicator(0, static_cast<int>(vals.size()) - 1,
                                         vals.data(),
                                         period,
                                         &outBeg,
                                         &outNbElement,
                                         raw_output.data());

    if (retCode != TA_SUCCESS)
    {
        std::cout << "error in " << indicator_name << std::endl;
        std::abort();
    }

    return build_talib_output(vals.size(), outBeg, outNbElement, raw_output, indicator_name);
}

std::vector<float> run_ohlc_indicator(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, OhlcInputTalibFn indicator, const char *indicator_name)
{
    if (close.empty())
    {
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw_output(close.size(), 0.0f);

    const TA_RetCode retCode = indicator(0, static_cast<int>(close.size()) - 1,
                                         high.data(), low.data(), close.data(),
                                         period,
                                         &outBeg,
                                         &outNbElement,
                                         raw_output.data());

    if (retCode != TA_SUCCESS)
    {
        std::cout << "error in " << indicator_name << std::endl;
        std::abort();
    }

    const std::vector<float> output = build_talib_output(close.size(), outBeg, outNbElement, raw_output, indicator_name);
    if (output.size() != high.size() || output.size() != low.size() || output.size() != close.size())
    {
        std::cout << "error in " << indicator_name << std::endl;
        std::cout << output.size() << " " << low.size() << " " << close.size() << " " << high.size() << std::endl;
        std::abort();
    }

    return output;
}
} // namespace
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_MIN(const std::vector<float> &vals, const int period)
{
    return run_single_input_indicator(vals, period, TA_S_MIN, "TALIB_MIN");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_MAX(const std::vector<float> &vals, const int period)
{
    return run_single_input_indicator(vals, period, TA_S_MAX, "TALIB_MAX");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_RSI(const std::vector<float> &vals, const int period)
{
    return run_single_input_indicator(vals, period, TA_S_RSI, "TALIB_RSI");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_EMA(const std::vector<float> &vals, const int period)
{
    return run_single_input_indicator(vals, period, TA_S_EMA, "TALIB_EMA");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_SMA(const std::vector<float> &vals, const int period)
{
    return run_single_input_indicator(vals, period, TA_S_SMA, "TALIB_SMA");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<float> TALIB_ATR(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period)
{
    return run_ohlc_indicator(high, low, close, period, TA_S_ATR, "TALIB_ATR");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Shared StochRSI prologue: RSI -> rolling min/max -> rounded normalized series. The
// three public entries (K, D, not-averaged) differ only in how many trailing SMAs they
// apply.
static std::vector<float> stoch_rsi_raw(const std::vector<float> &vals, const int nb_period_stoch, const int nb_period_rsi)
{
    const std::vector<float> rsi = TALIB_RSI(vals, nb_period_rsi);
    const std::vector<float> highest_rsi = TALIB_MAX(rsi, nb_period_stoch);
    const std::vector<float> lowest_rsi = TALIB_MIN(rsi, nb_period_stoch);

    std::vector<float> stochrsi;
    stochrsi.reserve(rsi.size());
    for (uint i = 0; i < rsi.size(); i++)
    {
        float val = (rsi[i] - lowest_rsi[i]) / (highest_rsi[i] - lowest_rsi[i]);
        if (std::isnan(val) || std::isinf(val))
        {
            val = 0.0f;
        }
        stochrsi.push_back(std::round(val * 1000.0f) / 1000.0f);
    }
    return stochrsi;
}

std::vector<float> TALIB_STOCHRSI_not_averaged(const std::vector<float> &vals, const int nb_period_stoch, const int nb_period_rsi)
{
    return stoch_rsi_raw(vals, nb_period_stoch, nb_period_rsi);
}

std::vector<float> TALIB_STOCHRSI_K(const std::vector<float> &vals, const int nb_period_stoch, const int nb_period_rsi, const int nb_period_k, const int /*nb_period_d*/)
{
    return TALIB_SMA(stoch_rsi_raw(vals, nb_period_stoch, nb_period_rsi), nb_period_k);
}

std::vector<float> TALIB_STOCHRSI_D(const std::vector<float> &vals, const int nb_period_stoch, const int nb_period_rsi, const int nb_period_k, const int nb_period_d)
{
    return TALIB_SMA(TALIB_SMA(stoch_rsi_raw(vals, nb_period_stoch, nb_period_rsi), nb_period_k), nb_period_d);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_TRIX(const std::vector<float> &vals, const int trixLength, const int trixSignal)
{
    std::vector<float> TRIX{};
    std::vector<float> TRIX_PCT{};
    std::vector<float> TRIX_SIGNAL{};
    std::vector<float> TRIX_HISTO{};

    TRIX = TALIB_EMA(vals, trixLength);
    TRIX = TALIB_EMA(TRIX, trixLength);
    TRIX = TALIB_EMA(TRIX, trixLength);

    TRIX_PCT.reserve(vals.size());
    TRIX_HISTO.reserve(vals.size());

    TRIX_PCT.push_back(0.0);
    for (uint i = 1; i < TRIX.size(); i++)
    {
        float val = (TRIX[i] - TRIX[i - 1]) / TRIX[i - 1] * 100.0;
        if (std::isinf(val) || std::isnan(val))
        {
            val = 0.0;
        }
        TRIX_PCT.push_back(val);
    }

    if (TRIX_PCT.size() != TRIX.size())
    {
        std::cout << "ERROR TRIX_PCT.size()!=TRIX.size()" << std::endl;
        std::abort();
    }

    TRIX_SIGNAL = TALIB_SMA(TRIX_PCT, trixSignal);

    for (uint i = 0; i < TRIX_PCT.size(); i++)
    {
        TRIX_HISTO.push_back(TRIX_PCT[i] - TRIX_SIGNAL[i]);
    }

    return TRIX_HISTO;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
SuperTrend TALIB_SuperTrend(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                            const int atr_window, const int atr_multi)
{
    vector<float> hl2{};
    vector<float> upperband{};
    vector<float> lowerband{};
    vector<int> dir{};
    vector<float> trend{};
    vector<float> longg{};
    vector<float> shortt{};
    vector<float> matr{};
    hl2.reserve(high.size());
    upperband.reserve(high.size());
    lowerband.reserve(high.size());
    dir.reserve(high.size());
    trend.reserve(high.size());
    longg.reserve(high.size());
    shortt.reserve(high.size());
    matr.reserve(high.size());

    const vector<float> ATR = TALIB_ATR(high, low, close, atr_window);

    const uint m = close.size();

    for (uint ii = 0; ii < m; ii++)
    {
        dir.push_back(1);
        trend.push_back(0.0);
        longg.push_back(NAN);
        shortt.push_back(NAN);
    }

    for (uint ii = 0; ii < m; ii++)
    {
        // HL2 is the average of high and low prices
        hl2.push_back((high[ii] + low[ii]) / 2.0f);
        matr.push_back(float(atr_multi) * ATR[ii]);
        upperband.push_back(hl2.back() + matr.back());
        lowerband.push_back(hl2.back() - matr.back());
    }

    for (uint ii = 1; ii < m; ii++)
    {
        if (close[ii] > upperband[ii - 1])
        {
            dir[ii] = 1;
        }
        else if (close[ii] < lowerband[ii - 1])
        {
            dir[ii] = -1;
        }
        else
        {
            dir[ii] = dir[ii - 1];
            if (dir[ii] > 0 && lowerband[ii] < lowerband[ii - 1])
            {
                lowerband[ii] = lowerband[ii - 1];
            }
            if (dir[ii] < 0 && upperband[ii] > upperband[ii - 1])
            {
                upperband[ii] = upperband[ii - 1];
            }
        }

        if (dir[ii] > 0)
        {
            trend[ii] = lowerband[ii];
            longg[ii] = lowerband[ii];
        }
        else
        {
            trend[ii] = upperband[ii];
            shortt[ii] = upperband[ii];
        }
    }

    return {dir, lowerband, upperband};
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<float> TALIB_SuperTrend_dir_only(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                                             const int atr_window, const int atr_multi)
{
    vector<float> hl2{};
    vector<float> upperband{};
    vector<float> lowerband{};
    vector<float> dir{};
    vector<float> trend{};
    vector<float> longg{};
    vector<float> shortt{};
    vector<float> matr{};
    hl2.reserve(high.size());
    upperband.reserve(high.size());
    lowerband.reserve(high.size());
    dir.reserve(high.size());
    trend.reserve(high.size());
    longg.reserve(high.size());
    shortt.reserve(high.size());
    matr.reserve(high.size());

    const vector<float> ATR = TALIB_ATR(high, low, close, atr_window);

    const uint m = close.size();

    for (uint ii = 0; ii < m; ii++)
    {
        dir.push_back(1.0);
        trend.push_back(0.0);
        longg.push_back(NAN);
        shortt.push_back(NAN);
    }

    for (uint ii = 0; ii < m; ii++)
    {
        // HL2 is the average of high and low prices
        hl2.push_back((high[ii] + low[ii]) / 2.0f);
        matr.push_back(float(atr_multi) * ATR[ii]);
        upperband.push_back(hl2.back() + matr.back());
        lowerband.push_back(hl2.back() - matr.back());
    }

    for (uint ii = 1; ii < m; ii++)
    {
        if (close[ii] > upperband[ii - 1])
        {
            dir[ii] = 1;
        }
        else if (close[ii] < lowerband[ii - 1])
        {
            dir[ii] = -1;
        }
        else
        {
            dir[ii] = dir[ii - 1];
            if (dir[ii] > 0 && lowerband[ii] < lowerband[ii - 1])
            {
                lowerband[ii] = lowerband[ii - 1];
            }
            if (dir[ii] < 0 && upperband[ii] > upperband[ii - 1])
            {
                upperband[ii] = upperband[ii - 1];
            }
        }

        if (dir[ii] > 0)
        {
            trend[ii] = lowerband[ii];
            longg[ii] = lowerband[ii];
        }
        else
        {
            trend[ii] = upperband[ii];
            shortt[ii] = upperband[ii];
        }
    }

    return dir;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_AO(const std::vector<float> &high, const std::vector<float> &low,
                            const int fast, const int slow)
{
    int fastt = fast;
    int sloww = slow;
    if (sloww < fastt)
    {
        std::swap(fastt, sloww);
    }

    const uint m = high.size();

    vector<float> median_price{};
    vector<float> fast_sma{};
    vector<float> slow_sma{};
    vector<float> ao{};
    median_price.reserve(m);
    ao.reserve(m);

    for (uint ii = 0; ii < m; ii++)
    {
        median_price.push_back(0.5f * (high[ii] + low[ii]));
    }

    fast_sma = TALIB_SMA(median_price, fastt);
    slow_sma = TALIB_SMA(median_price, sloww);
    for (uint ii = 0; ii < m; ii++)
    {
        ao.push_back(fast_sma[ii] - slow_sma[ii]);
    }

    return ao;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<float> TALIB_WILLR(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                               const int length)
{
    return run_ohlc_indicator(high, low, close, length, TA_S_WILLR, "TALIB_WILLR");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void TALIB_BBANDS(const std::vector<float> &close,
                  const float &optInNbDevUp, const float &optInNbDevDn, const int &length,
                  std::vector<float> &OUT_u, std::vector<float> &OUT_m, std::vector<float> &OUT_l)
{
    if (close.empty())
    {
        OUT_u.clear();
        OUT_m.clear();
        OUT_l.clear();
        return;
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> upper_band(close.size(), 0.0f);
    std::vector<TA_Real> middle_band(close.size(), 0.0f);
    std::vector<TA_Real> lower_band(close.size(), 0.0f);

    const TA_RetCode retCode = TA_S_BBANDS(0, static_cast<int>(close.size()) - 1,
                                           close.data(), length,
                                           optInNbDevUp, optInNbDevDn, TA_MAType_SMA,
                                           &outBeg,
                                           &outNbElement,
                                           upper_band.data(), middle_band.data(), lower_band.data());

    if (retCode != TA_SUCCESS)
    {
        std::cout << "error in TALIB_BBANDS" << std::endl;
        std::abort();
    }

    OUT_u = build_talib_output(close.size(), outBeg, outNbElement, upper_band, "TALIB_BBANDS");
    OUT_m = build_talib_output(close.size(), outBeg, outNbElement, middle_band, "TALIB_BBANDS");
    OUT_l = build_talib_output(close.size(), outBeg, outNbElement, lower_band, "TALIB_BBANDS");

    if (OUT_u.size() != close.size() || OUT_m.size() != close.size() || OUT_l.size() != close.size())
    {
        std::cout << "error in TALIB_BBANDS" << std::endl;
        std::abort();
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void print_first_elements(const std::vector<uint> &vec, const int nb)
{
    std::cout << vec.size() << std::endl;
    for (uint i = 0; i < nb; i++)
    {
        std::cout << vec[i] << std::endl;
    }
    std::cout << std::endl;
}

void print_last_elements(const std::vector<uint> &vec, const int nb)
{
    std::cout << vec.size() << std::endl;
    for (uint i = vec.size() - 1; i > vec.size() - nb - 1; i--)
    {
        std::cout << vec[i] << std::endl;
    }
    std::cout << std::endl;
}
