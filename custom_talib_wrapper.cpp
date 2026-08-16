#include "custom_talib_wrapper.hh"
using namespace std;

namespace
{
using SingleInputTalibFn = TA_RetCode (*)(int, int, const float *, int, TA_Integer *, TA_Integer *, TA_Real *);
using OhlcInputTalibFn = TA_RetCode (*)(int, int, const float *, const float *, const float *, int, TA_Integer *, TA_Integer *, TA_Real *);

// Scatters TA-Lib's compacted output back into a full-length series, left-padding the
// warmup region with zeros. The bounds check matters: TA-Lib reports where its output
// starts (outBeg) and how many values it produced (outNbElement), and a wrong pair
// would otherwise write straight past the end of `output`.
std::vector<float> build_talib_output(const size_t input_size, const TA_Integer outBeg, const TA_Integer outNbElement, const std::vector<TA_Real> &raw_output, const char *indicator_name)
{
    if (outBeg < 0 || outNbElement < 0 ||
        static_cast<size_t>(outBeg) + static_cast<size_t>(outNbElement) > input_size ||
        static_cast<size_t>(outNbElement) > raw_output.size())
    {
        std::cerr << "FATAL " << indicator_name << ": TA-Lib returned out-of-range output window"
                  << " (outBeg=" << outBeg << ", outNbElement=" << outNbElement
                  << ", input_size=" << input_size << ")" << std::endl;
        std::abort();
    }

    std::vector<float> output(input_size, 0.0f);
    for (TA_Integer ii = 0; ii < outNbElement; ++ii)
    {
        output[static_cast<size_t>(outBeg + ii)] = raw_output[static_cast<size_t>(ii)];
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
    // Checked before the call, not after: TA-Lib reads high/low/close up to close.size(),
    // so a short high or low array would already have been read out of bounds by the
    // time a post-hoc size comparison could notice.
    if (high.size() != close.size() || low.size() != close.size())
    {
        std::cerr << "FATAL " << indicator_name << ": OHLC input size mismatch"
                  << " (high=" << high.size() << ", low=" << low.size() << ", close=" << close.size() << ")" << std::endl;
        std::abort();
    }

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

    return build_talib_output(close.size(), outBeg, outNbElement, raw_output, indicator_name);
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
// Single implementation behind both public entry points. The `_dir_only` variant used
// to be a verbatim copy of this function differing only in its return type, and both
// copies also built three series (trend/long/short) that neither of them returned.
SuperTrend TALIB_SuperTrend(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                            const int atr_window, const float atr_multi)
{
    const uint m = close.size();

    vector<float> upperband{};
    vector<float> lowerband{};
    vector<int> dir(m, 1);
    upperband.reserve(m);
    lowerband.reserve(m);

    const vector<float> ATR = TALIB_ATR(high, low, close, atr_window);

    for (uint ii = 0; ii < m; ii++)
    {
        // HL2, the midpoint of the bar, is the band centreline.
        const float hl2 = (high[ii] + low[ii]) / 2.0f;
        const float matr = atr_multi * ATR[ii];
        upperband.push_back(hl2 + matr);
        lowerband.push_back(hl2 - matr);
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
            // Inside the bands: hold the previous direction and ratchet the active band
            // so it only ever moves in the trend's favour.
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
    }

    return {dir, lowerband, upperband};
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<float> TALIB_SuperTrend_dir_only(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                                             const int atr_window, const float atr_multi)
{
    const std::vector<int> dir = TALIB_SuperTrend(high, low, close, atr_window, atr_multi).supertrend;
    return std::vector<float>(dir.begin(), dir.end());
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
