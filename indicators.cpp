#include "indicators.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>

using namespace std;

// ---------------------------------------------------------------------------
//  Generic TA-Lib adapters.
//
//  TA-Lib's C API is uniform in shape but verbose: every call needs an output
//  buffer, two out-parameters, a return-code check, and a scatter back to full
//  length. These adapters do that once per input shape, so each public wrapper
//  below is a single line naming the TA_S_* function it forwards to.
// ---------------------------------------------------------------------------
namespace
{
using SingleInputTalibFn = TA_RetCode (*)(int, int, const float *, int, TA_Integer *, TA_Integer *, TA_Real *);
using OhlcInputTalibFn = TA_RetCode (*)(int, int, const float *, const float *, const float *, int, TA_Integer *, TA_Integer *, TA_Real *);
using OhlcNoPeriodTalibFn = TA_RetCode (*)(int, int, const float *, const float *, const float *, TA_Integer *, TA_Integer *, TA_Real *);
using OhlcvTalibFn = TA_RetCode (*)(int, int, const float *, const float *, const float *, const float *, int, TA_Integer *, TA_Integer *, TA_Real *);

[[noreturn]] void indicator_fatal(const char *indicator_name, const std::string &detail)
{
    std::cerr << "FATAL " << indicator_name << ": " << detail << std::endl;
    std::abort();
}

void check_ret(const TA_RetCode retCode, const char *indicator_name)
{
    if (retCode != TA_SUCCESS)
    {
        indicator_fatal(indicator_name, "TA-Lib returned error code " + std::to_string(static_cast<int>(retCode)));
    }
}

// Requires all named series to be the same length, checked *before* the TA-Lib call:
// TA-Lib indexes every input up to the range end, so a short array would already have
// been read out of bounds by the time a post-hoc comparison could notice.
void check_same_size(const char *indicator_name, std::initializer_list<std::pair<const char *, size_t>> series)
{
    if (series.size() == 0)
    {
        return;
    }
    const size_t expected = series.begin()->second;
    for (const auto &s : series)
    {
        if (s.second != expected)
            {
            std::string detail = "input size mismatch:";
            for (const auto &t : series)
            {
                detail += " " + std::string(t.first) + "=" + std::to_string(t.second);
            }
            indicator_fatal(indicator_name, detail);
        }
    }
}

// Scatters TA-Lib's compacted output back into a full-length series, left-padding the
// warmup region with zeros. The bounds check matters: a wrong (outBeg, outNbElement)
// pair would otherwise write straight past the end of `output`.
std::vector<float> build_talib_output(const size_t input_size, const TA_Integer outBeg, const TA_Integer outNbElement,
                                      const std::vector<TA_Real> &raw_output, const char *indicator_name, size_t *warmup)
{
    if (outBeg < 0 || outNbElement < 0 ||
        static_cast<size_t>(outBeg) + static_cast<size_t>(outNbElement) > input_size ||
        static_cast<size_t>(outNbElement) > raw_output.size())
    {
        indicator_fatal(indicator_name, "TA-Lib returned an out-of-range output window (outBeg=" +
                                            std::to_string(outBeg) + ", outNbElement=" + std::to_string(outNbElement) +
                                            ", input_size=" + std::to_string(input_size) + ")");
    }

    std::vector<float> output(input_size, 0.0f);
    for (TA_Integer ii = 0; ii < outNbElement; ++ii)
    {
        output[static_cast<size_t>(outBeg + ii)] = static_cast<float>(raw_output[static_cast<size_t>(ii)]);
    }

    if (warmup != nullptr)
    {
        *warmup = static_cast<size_t>(outBeg);
    }

    return output;
}

void set_warmup(size_t *warmup, const size_t value)
{
    if (warmup != nullptr)
    {
        *warmup = value;
    }
}

std::vector<float> run_single_input(const std::vector<float> &vals, const int period, SingleInputTalibFn indicator,
                                    const char *indicator_name, size_t *warmup)
{
    if (vals.empty())
    {
        set_warmup(warmup, 0);
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw_output(vals.size(), 0.0);

    check_ret(indicator(0, static_cast<int>(vals.size()) - 1, vals.data(), period, &outBeg, &outNbElement, raw_output.data()),
              indicator_name);

    return build_talib_output(vals.size(), outBeg, outNbElement, raw_output, indicator_name, warmup);
}

std::vector<float> run_ohlc(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                            const int period, OhlcInputTalibFn indicator, const char *indicator_name, size_t *warmup)
{
    check_same_size(indicator_name, {{"high", high.size()}, {"low", low.size()}, {"close", close.size()}});
    if (close.empty())
    {
        set_warmup(warmup, 0);
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw_output(close.size(), 0.0);

    check_ret(indicator(0, static_cast<int>(close.size()) - 1, high.data(), low.data(), close.data(), period,
                        &outBeg, &outNbElement, raw_output.data()),
              indicator_name);

    return build_talib_output(close.size(), outBeg, outNbElement, raw_output, indicator_name, warmup);
}

std::vector<float> run_ohlc_no_period(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                                      OhlcNoPeriodTalibFn indicator, const char *indicator_name, size_t *warmup)
{
    check_same_size(indicator_name, {{"high", high.size()}, {"low", low.size()}, {"close", close.size()}});
    if (close.empty())
    {
        set_warmup(warmup, 0);
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw_output(close.size(), 0.0);

    check_ret(indicator(0, static_cast<int>(close.size()) - 1, high.data(), low.data(), close.data(),
                        &outBeg, &outNbElement, raw_output.data()),
              indicator_name);

    return build_talib_output(close.size(), outBeg, outNbElement, raw_output, indicator_name, warmup);
}

std::vector<float> run_ohlcv(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                             const std::vector<float> &volume, const int period, OhlcvTalibFn indicator,
                             const char *indicator_name, size_t *warmup)
{
    check_same_size(indicator_name, {{"high", high.size()}, {"low", low.size()}, {"close", close.size()}, {"volume", volume.size()}});
    if (close.empty())
    {
        set_warmup(warmup, 0);
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw_output(close.size(), 0.0);

    check_ret(indicator(0, static_cast<int>(close.size()) - 1, high.data(), low.data(), close.data(), volume.data(),
                        period, &outBeg, &outNbElement, raw_output.data()),
              indicator_name);

    return build_talib_output(close.size(), outBeg, outNbElement, raw_output, indicator_name, warmup);
}
} // namespace

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Price transforms
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> PRICE_HL2(const std::vector<float> &high, const std::vector<float> &low)
{
    check_same_size("PRICE_HL2", {{"high", high.size()}, {"low", low.size()}});
    std::vector<float> out;
    out.reserve(high.size());
    for (size_t i = 0; i < high.size(); ++i)
    {
        out.push_back(0.5f * (high[i] + low[i]));
    }
    return out;
}

std::vector<float> PRICE_HLC3(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close)
{
    check_same_size("PRICE_HLC3", {{"high", high.size()}, {"low", low.size()}, {"close", close.size()}});
    std::vector<float> out;
    out.reserve(high.size());
    for (size_t i = 0; i < high.size(); ++i)
    {
        out.push_back((high[i] + low[i] + close[i]) / 3.0f);
    }
    return out;
}

std::vector<float> PRICE_OHLC4(const std::vector<float> &open, const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close)
{
    check_same_size("PRICE_OHLC4", {{"open", open.size()}, {"high", high.size()}, {"low", low.size()}, {"close", close.size()}});
    std::vector<float> out;
    out.reserve(open.size());
    for (size_t i = 0; i < open.size(); ++i)
    {
        out.push_back(0.25f * (open[i] + high[i] + low[i] + close[i]));
    }
    return out;
}

HeikinAshi HEIKIN_ASHI(const std::vector<float> &open, const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close)
{
    check_same_size("HEIKIN_ASHI", {{"open", open.size()}, {"high", high.size()}, {"low", low.size()}, {"close", close.size()}});

    HeikinAshi ha{};
    const size_t n = close.size();
    if (n == 0)
    {
        return ha;
    }

    ha.open.reserve(n);
    ha.high.reserve(n);
    ha.low.reserve(n);
    ha.close.reserve(n);

    for (size_t i = 0; i < n; ++i)
    {
        const float ha_close = 0.25f * (open[i] + high[i] + low[i] + close[i]);
        // The first HA open seeds from the raw candle; afterwards it is the running
        // average of the previous HA candle, which is what smooths the series.
        const float ha_open = (i == 0) ? 0.5f * (open[0] + close[0])
                                       : 0.5f * (ha.open[i - 1] + ha.close[i - 1]);
        ha.open.push_back(ha_open);
        ha.close.push_back(ha_close);
        ha.high.push_back(std::max({high[i], ha_open, ha_close}));
        ha.low.push_back(std::min({low[i], ha_open, ha_close}));
    }

    // Bar 0 is seeded rather than derived, so treat it as warmup.
    ha.warmup = 1;
    return ha;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Rolling extremes
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_MIN(const std::vector<float> &vals, const int period, size_t *warmup)
{
    return run_single_input(vals, period, TA_S_MIN, "TALIB_MIN", warmup);
}

std::vector<float> TALIB_MAX(const std::vector<float> &vals, const int period, size_t *warmup)
{
    return run_single_input(vals, period, TA_S_MAX, "TALIB_MAX", warmup);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Moving averages
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_SMA(const std::vector<float> &vals, const int period, size_t *warmup)
{
    return run_single_input(vals, period, TA_S_SMA, "TALIB_SMA", warmup);
}

std::vector<float> TALIB_EMA(const std::vector<float> &vals, const int period, size_t *warmup)
{
    return run_single_input(vals, period, TA_S_EMA, "TALIB_EMA", warmup);
}

std::vector<float> TALIB_WMA(const std::vector<float> &vals, const int period, size_t *warmup)
{
    return run_single_input(vals, period, TA_S_WMA, "TALIB_WMA", warmup);
}

std::vector<float> TALIB_DEMA(const std::vector<float> &vals, const int period, size_t *warmup)
{
    return run_single_input(vals, period, TA_S_DEMA, "TALIB_DEMA", warmup);
}

std::vector<float> TALIB_TEMA(const std::vector<float> &vals, const int period, size_t *warmup)
{
    return run_single_input(vals, period, TA_S_TEMA, "TALIB_TEMA", warmup);
}

std::vector<float> TALIB_KAMA(const std::vector<float> &vals, const int period, size_t *warmup)
{
    return run_single_input(vals, period, TA_S_KAMA, "TALIB_KAMA", warmup);
}

// Hull MA = WMA(2*WMA(n/2) - WMA(n), sqrt(n)). TA-Lib has no HMA, so it is composed
// here; the intermediate difference series is full length, so warmup accumulates.
std::vector<float> TALIB_HMA(const std::vector<float> &vals, const int period, size_t *warmup)
{
    if (vals.empty() || period < 2)
    {
        set_warmup(warmup, 0);
        return vals.empty() ? std::vector<float>{} : TALIB_WMA(vals, std::max(period, 1), warmup);
    }

    const int half = std::max(1, period / 2);
    const int sqrt_period = std::max(1, static_cast<int>(std::lround(std::sqrt(static_cast<double>(period)))));

    size_t warm_half = 0;
    size_t warm_full = 0;
    const std::vector<float> wma_half = TALIB_WMA(vals, half, &warm_half);
    const std::vector<float> wma_full = TALIB_WMA(vals, period, &warm_full);

    std::vector<float> diff;
    diff.reserve(vals.size());
    for (size_t i = 0; i < vals.size(); ++i)
    {
        diff.push_back(2.0f * wma_half[i] - wma_full[i]);
    }

    size_t warm_final = 0;
    std::vector<float> out = TALIB_WMA(diff, sqrt_period, &warm_final);

    // The difference series is only meaningful once the slower WMA has warmed up, so
    // the real warmup is the slower input's plus the final smoothing's.
    const size_t total_warmup = std::min(vals.size(), std::max(warm_half, warm_full) + warm_final);
    for (size_t i = 0; i < total_warmup; ++i)
    {
        out[i] = 0.0f;
    }
    set_warmup(warmup, total_warmup);
    return out;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Momentum / oscillators
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_RSI(const std::vector<float> &vals, const int period, size_t *warmup)
{
    return run_single_input(vals, period, TA_S_RSI, "TALIB_RSI", warmup);
}

std::vector<float> TALIB_ROC(const std::vector<float> &vals, const int period, size_t *warmup)
{
    return run_single_input(vals, period, TA_S_ROC, "TALIB_ROC", warmup);
}

std::vector<float> TALIB_MOM(const std::vector<float> &vals, const int period, size_t *warmup)
{
    return run_single_input(vals, period, TA_S_MOM, "TALIB_MOM", warmup);
}

std::vector<float> TALIB_CCI(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup)
{
    return run_ohlc(high, low, close, period, TA_S_CCI, "TALIB_CCI", warmup);
}

std::vector<float> TALIB_WILLR(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int length, size_t *warmup)
{
    return run_ohlc(high, low, close, length, TA_S_WILLR, "TALIB_WILLR", warmup);
}

std::vector<float> TALIB_ULTOSC(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                                const int period1, const int period2, const int period3, size_t *warmup)
{
    check_same_size("TALIB_ULTOSC", {{"high", high.size()}, {"low", low.size()}, {"close", close.size()}});
    if (close.empty())
    {
        set_warmup(warmup, 0);
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw(close.size(), 0.0);
    check_ret(TA_S_ULTOSC(0, static_cast<int>(close.size()) - 1, high.data(), low.data(), close.data(),
                          period1, period2, period3, &outBeg, &outNbElement, raw.data()),
              "TALIB_ULTOSC");
    return build_talib_output(close.size(), outBeg, outNbElement, raw, "TALIB_ULTOSC", warmup);
}

MACDResult TALIB_MACD(const std::vector<float> &vals, const int fast_period, const int slow_period, const int signal_period)
{
    MACDResult out{};
    if (vals.empty())
    {
        return out;
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> macd(vals.size(), 0.0);
    std::vector<TA_Real> signal(vals.size(), 0.0);
    std::vector<TA_Real> hist(vals.size(), 0.0);

    check_ret(TA_S_MACD(0, static_cast<int>(vals.size()) - 1, vals.data(), fast_period, slow_period, signal_period,
                        &outBeg, &outNbElement, macd.data(), signal.data(), hist.data()),
              "TALIB_MACD");

    out.macd = build_talib_output(vals.size(), outBeg, outNbElement, macd, "TALIB_MACD", &out.warmup);
    out.signal = build_talib_output(vals.size(), outBeg, outNbElement, signal, "TALIB_MACD", nullptr);
    out.histogram = build_talib_output(vals.size(), outBeg, outNbElement, hist, "TALIB_MACD", nullptr);
    return out;
}

StochResult TALIB_STOCH(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                        const int fastk_period, const int slowk_period, const int slowd_period)
{
    check_same_size("TALIB_STOCH", {{"high", high.size()}, {"low", low.size()}, {"close", close.size()}});
    StochResult out{};
    if (close.empty())
    {
        return out;
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> k(close.size(), 0.0);
    std::vector<TA_Real> d(close.size(), 0.0);

    check_ret(TA_S_STOCH(0, static_cast<int>(close.size()) - 1, high.data(), low.data(), close.data(),
                         fastk_period, slowk_period, TA_MAType_SMA, slowd_period, TA_MAType_SMA,
                         &outBeg, &outNbElement, k.data(), d.data()),
              "TALIB_STOCH");

    out.k = build_talib_output(close.size(), outBeg, outNbElement, k, "TALIB_STOCH", &out.warmup);
    out.d = build_talib_output(close.size(), outBeg, outNbElement, d, "TALIB_STOCH", nullptr);
    return out;
}

AroonResult TALIB_AROON(const std::vector<float> &high, const std::vector<float> &low, const int period)
{
    check_same_size("TALIB_AROON", {{"high", high.size()}, {"low", low.size()}});
    AroonResult out{};
    if (high.empty())
    {
        return out;
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> down(high.size(), 0.0);
    std::vector<TA_Real> up(high.size(), 0.0);

    check_ret(TA_S_AROON(0, static_cast<int>(high.size()) - 1, high.data(), low.data(), period,
                         &outBeg, &outNbElement, down.data(), up.data()),
              "TALIB_AROON");

    out.down = build_talib_output(high.size(), outBeg, outNbElement, down, "TALIB_AROON", &out.warmup);
    out.up = build_talib_output(high.size(), outBeg, outNbElement, up, "TALIB_AROON", nullptr);
    return out;
}

std::vector<float> TALIB_AO(const std::vector<float> &high, const std::vector<float> &low, const int fast, const int slow)
{
    int fastt = fast;
    int sloww = slow;
    if (sloww < fastt)
    {
        std::swap(fastt, sloww);
    }

    const std::vector<float> median_price = PRICE_HL2(high, low);
    const std::vector<float> fast_sma = TALIB_SMA(median_price, fastt);
    const std::vector<float> slow_sma = TALIB_SMA(median_price, sloww);

    std::vector<float> ao;
    ao.reserve(high.size());
    for (size_t ii = 0; ii < high.size(); ii++)
    {
        ao.push_back(fast_sma[ii] - slow_sma[ii]);
    }
    return ao;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  StochRSI
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Shared prologue: RSI -> rolling min/max -> normalized series. The three public
// entries differ only in how many trailing SMAs they apply.
//
// The round-to-3-decimals is inherited behaviour, not a numerical necessity; it is
// documented in the header because it makes these differ from TA-Lib's own STOCHRSI.
static std::vector<float> stoch_rsi_raw(const std::vector<float> &vals, const int nb_period_stoch, const int nb_period_rsi)
{
    const std::vector<float> rsi = TALIB_RSI(vals, nb_period_rsi);
    const std::vector<float> highest_rsi = TALIB_MAX(rsi, nb_period_stoch);
    const std::vector<float> lowest_rsi = TALIB_MIN(rsi, nb_period_stoch);

    std::vector<float> stochrsi;
    stochrsi.reserve(rsi.size());
    for (size_t i = 0; i < rsi.size(); i++)
    {
        // Flat RSI over the window makes this 0/0; the warmup region makes it 0/0 too.
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
//  Trend strength
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_ADX(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup)
{
    return run_ohlc(high, low, close, period, TA_S_ADX, "TALIB_ADX", warmup);
}

std::vector<float> TALIB_PLUS_DI(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup)
{
    return run_ohlc(high, low, close, period, TA_S_PLUS_DI, "TALIB_PLUS_DI", warmup);
}

std::vector<float> TALIB_MINUS_DI(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup)
{
    return run_ohlc(high, low, close, period, TA_S_MINUS_DI, "TALIB_MINUS_DI", warmup);
}

DirectionalResult TALIB_DMI(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period)
{
    DirectionalResult out{};
    size_t warm_adx = 0;
    out.adx = TALIB_ADX(high, low, close, period, &warm_adx);
    out.plus_di = TALIB_PLUS_DI(high, low, close, period, nullptr);
    out.minus_di = TALIB_MINUS_DI(high, low, close, period, nullptr);
    // ADX smooths DX on top of the DI warmup, so its lookback dominates the group.
    out.warmup = warm_adx;
    return out;
}

std::vector<float> TALIB_SAR(const std::vector<float> &high, const std::vector<float> &low, const double acceleration, const double maximum, size_t *warmup)
{
    check_same_size("TALIB_SAR", {{"high", high.size()}, {"low", low.size()}});
    if (high.empty())
    {
        set_warmup(warmup, 0);
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw(high.size(), 0.0);
    check_ret(TA_S_SAR(0, static_cast<int>(high.size()) - 1, high.data(), low.data(), acceleration, maximum,
                       &outBeg, &outNbElement, raw.data()),
              "TALIB_SAR");
    return build_talib_output(high.size(), outBeg, outNbElement, raw, "TALIB_SAR", warmup);
}

std::vector<float> TALIB_TRIX(const std::vector<float> &vals, const int trixLength, const int trixSignal)
{
    // Triple-smoothed EMA, then its percent change, then that minus its own SMA. This
    // is the "TRIX histogram" convention, not TA-Lib's TA_TRIX (which returns the
    // percent-change line alone).
    std::vector<float> TRIX = TALIB_EMA(vals, trixLength);
    TRIX = TALIB_EMA(TRIX, trixLength);
    TRIX = TALIB_EMA(TRIX, trixLength);

    std::vector<float> TRIX_PCT;
    TRIX_PCT.reserve(vals.size());
    TRIX_PCT.push_back(0.0f);
    for (size_t i = 1; i < TRIX.size(); i++)
    {
        float val = (TRIX[i] - TRIX[i - 1]) / TRIX[i - 1] * 100.0f;
        if (std::isinf(val) || std::isnan(val))
        {
            val = 0.0f;
        }
        TRIX_PCT.push_back(val);
    }

    const std::vector<float> TRIX_SIGNAL = TALIB_SMA(TRIX_PCT, trixSignal);

    std::vector<float> TRIX_HISTO;
    TRIX_HISTO.reserve(TRIX_PCT.size());
    for (size_t i = 0; i < TRIX_PCT.size(); i++)
    {
        TRIX_HISTO.push_back(TRIX_PCT[i] - TRIX_SIGNAL[i]);
    }

    return TRIX_HISTO;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Volatility
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_ATR(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup)
{
    return run_ohlc(high, low, close, period, TA_S_ATR, "TALIB_ATR", warmup);
}

std::vector<float> TALIB_NATR(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup)
{
    return run_ohlc(high, low, close, period, TA_S_NATR, "TALIB_NATR", warmup);
}

std::vector<float> TALIB_TRANGE(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, size_t *warmup)
{
    return run_ohlc_no_period(high, low, close, TA_S_TRANGE, "TALIB_TRANGE", warmup);
}

std::vector<float> TALIB_STDDEV(const std::vector<float> &vals, const int period, const double nb_dev, size_t *warmup)
{
    if (vals.empty())
    {
        set_warmup(warmup, 0);
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw(vals.size(), 0.0);
    check_ret(TA_S_STDDEV(0, static_cast<int>(vals.size()) - 1, vals.data(), period, nb_dev,
                          &outBeg, &outNbElement, raw.data()),
              "TALIB_STDDEV");
    return build_talib_output(vals.size(), outBeg, outNbElement, raw, "TALIB_STDDEV", warmup);
}

BandsResult TALIB_BBANDS_R(const std::vector<float> &close, const float dev_up, const float dev_dn, const int length)
{
    BandsResult out{};
    if (close.empty())
    {
        return out;
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> upper(close.size(), 0.0);
    std::vector<TA_Real> middle(close.size(), 0.0);
    std::vector<TA_Real> lower(close.size(), 0.0);

    check_ret(TA_S_BBANDS(0, static_cast<int>(close.size()) - 1, close.data(), length, dev_up, dev_dn, TA_MAType_SMA,
                          &outBeg, &outNbElement, upper.data(), middle.data(), lower.data()),
              "TALIB_BBANDS");

    out.upper = build_talib_output(close.size(), outBeg, outNbElement, upper, "TALIB_BBANDS", &out.warmup);
    out.middle = build_talib_output(close.size(), outBeg, outNbElement, middle, "TALIB_BBANDS", nullptr);
    out.lower = build_talib_output(close.size(), outBeg, outNbElement, lower, "TALIB_BBANDS", nullptr);
    return out;
}

void TALIB_BBANDS(const std::vector<float> &close,
                  const float &optInNbDevUp, const float &optInNbDevDn, const int &length,
                  std::vector<float> &OUT_u, std::vector<float> &OUT_m, std::vector<float> &OUT_l)
{
    BandsResult r = TALIB_BBANDS_R(close, optInNbDevUp, optInNbDevDn, length);
    OUT_u = std::move(r.upper);
    OUT_m = std::move(r.middle);
    OUT_l = std::move(r.lower);
}

BandsResult KELTNER_CHANNELS(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                             const int ema_period, const int atr_period, const float atr_multiplier)
{
    check_same_size("KELTNER_CHANNELS", {{"high", high.size()}, {"low", low.size()}, {"close", close.size()}});
    BandsResult out{};
    if (close.empty())
    {
        return out;
    }

    size_t warm_ema = 0;
    size_t warm_atr = 0;
    const std::vector<float> ema = TALIB_EMA(close, ema_period, &warm_ema);
    const std::vector<float> atr = TALIB_ATR(high, low, close, atr_period, &warm_atr);

    out.warmup = std::max(warm_ema, warm_atr);
    out.middle = ema;
    out.upper.reserve(close.size());
    out.lower.reserve(close.size());
    for (size_t i = 0; i < close.size(); ++i)
    {
        out.upper.push_back(ema[i] + atr_multiplier * atr[i]);
        out.lower.push_back(ema[i] - atr_multiplier * atr[i]);
    }
    return out;
}

BandsResult DONCHIAN_CHANNELS(const std::vector<float> &high, const std::vector<float> &low, const int period)
{
    check_same_size("DONCHIAN_CHANNELS", {{"high", high.size()}, {"low", low.size()}});
    BandsResult out{};
    if (high.empty())
    {
        return out;
    }

    size_t warm_hi = 0;
    size_t warm_lo = 0;
    out.upper = TALIB_MAX(high, period, &warm_hi);
    out.lower = TALIB_MIN(low, period, &warm_lo);
    out.warmup = std::max(warm_hi, warm_lo);

    out.middle.reserve(high.size());
    for (size_t i = 0; i < high.size(); ++i)
    {
        out.middle.push_back(0.5f * (out.upper[i] + out.lower[i]));
    }
    return out;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Volume
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_OBV(const std::vector<float> &close, const std::vector<float> &volume, size_t *warmup)
{
    check_same_size("TALIB_OBV", {{"close", close.size()}, {"volume", volume.size()}});
    if (close.empty())
    {
        set_warmup(warmup, 0);
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw(close.size(), 0.0);
    check_ret(TA_S_OBV(0, static_cast<int>(close.size()) - 1, close.data(), volume.data(),
                       &outBeg, &outNbElement, raw.data()),
              "TALIB_OBV");
    return build_talib_output(close.size(), outBeg, outNbElement, raw, "TALIB_OBV", warmup);
}

std::vector<float> TALIB_MFI(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                             const std::vector<float> &volume, const int period, size_t *warmup)
{
    return run_ohlcv(high, low, close, volume, period, TA_S_MFI, "TALIB_MFI", warmup);
}

std::vector<float> TALIB_AD(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                            const std::vector<float> &volume, size_t *warmup)
{
    check_same_size("TALIB_AD", {{"high", high.size()}, {"low", low.size()}, {"close", close.size()}, {"volume", volume.size()}});
    if (close.empty())
    {
        set_warmup(warmup, 0);
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw(close.size(), 0.0);
    check_ret(TA_S_AD(0, static_cast<int>(close.size()) - 1, high.data(), low.data(), close.data(), volume.data(),
                      &outBeg, &outNbElement, raw.data()),
              "TALIB_AD");
    return build_talib_output(close.size(), outBeg, outNbElement, raw, "TALIB_AD", warmup);
}

std::vector<float> TALIB_ADOSC(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                               const std::vector<float> &volume, const int fast_period, const int slow_period, size_t *warmup)
{
    check_same_size("TALIB_ADOSC", {{"high", high.size()}, {"low", low.size()}, {"close", close.size()}, {"volume", volume.size()}});
    if (close.empty())
    {
        set_warmup(warmup, 0);
        return {};
    }

    TA_Integer outBeg = 0;
    TA_Integer outNbElement = 0;
    std::vector<TA_Real> raw(close.size(), 0.0);
    check_ret(TA_S_ADOSC(0, static_cast<int>(close.size()) - 1, high.data(), low.data(), close.data(), volume.data(),
                         fast_period, slow_period, &outBeg, &outNbElement, raw.data()),
              "TALIB_ADOSC");
    return build_talib_output(close.size(), outBeg, outNbElement, raw, "TALIB_ADOSC", warmup);
}

std::vector<float> VWAP_ROLLING(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                                const std::vector<float> &volume, const int period, size_t *warmup)
{
    check_same_size("VWAP_ROLLING", {{"high", high.size()}, {"low", low.size()}, {"close", close.size()}, {"volume", volume.size()}});
    const size_t n = close.size();
    std::vector<float> out(n, 0.0f);
    if (n == 0 || period < 1)
    {
        set_warmup(warmup, 0);
        return out;
    }

    const std::vector<float> typical = PRICE_HLC3(high, low, close);
    const size_t window = static_cast<size_t>(period);

    // Running sums rather than a per-bar re-scan: this is O(n), and it is accumulated
    // in double because a long window of large volumes overflows float precision.
    double pv_sum = 0.0;
    double v_sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        pv_sum += static_cast<double>(typical[i]) * static_cast<double>(volume[i]);
        v_sum += static_cast<double>(volume[i]);

        if (i >= window)
        {
            pv_sum -= static_cast<double>(typical[i - window]) * static_cast<double>(volume[i - window]);
            v_sum -= static_cast<double>(volume[i - window]);
        }

        if (i + 1 >= window)
        {
            // Zero-volume windows exist in thin data; fall back to typical price.
            out[i] = (v_sum > 0.0) ? static_cast<float>(pv_sum / v_sum) : typical[i];
        }
    }

    set_warmup(warmup, std::min(n, window - 1));
    return out;
}

std::vector<float> RELATIVE_VOLUME(const std::vector<float> &volume, const int period, size_t *warmup)
{
    size_t warm = 0;
    const std::vector<float> avg = TALIB_SMA(volume, period, &warm);

    std::vector<float> out(volume.size(), 0.0f);
    for (size_t i = warm; i < volume.size(); ++i)
    {
        out[i] = (avg[i] > 0.0f) ? volume[i] / avg[i] : 0.0f;
    }
    set_warmup(warmup, warm);
    return out;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  SuperTrend
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Single implementation behind both public entry points. The `_dir_only` variant used
// to be a verbatim copy differing only in its return type, and both copies also built
// three series (trend/long/short) that neither of them returned.
SuperTrend TALIB_SuperTrend(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                            const int atr_window, const float atr_multi)
{
    const size_t m = close.size();

    SuperTrend out{};
    out.supertrend.assign(m, 1);
    out.final_lowerband.reserve(m);
    out.final_upperband.reserve(m);

    size_t warm_atr = 0;
    const std::vector<float> ATR = TALIB_ATR(high, low, close, atr_window, &warm_atr);
    out.warmup = warm_atr;

    std::vector<float> &upperband = out.final_upperband;
    std::vector<float> &lowerband = out.final_lowerband;

    for (size_t ii = 0; ii < m; ii++)
    {
        // HL2, the midpoint of the bar, is the band centreline.
        const float hl2 = (high[ii] + low[ii]) / 2.0f;
        const float matr = atr_multi * ATR[ii];
        upperband.push_back(hl2 + matr);
        lowerband.push_back(hl2 - matr);
    }

    for (size_t ii = 1; ii < m; ii++)
    {
        if (close[ii] > upperband[ii - 1])
        {
            out.supertrend[ii] = 1;
        }
        else if (close[ii] < lowerband[ii - 1])
        {
            out.supertrend[ii] = -1;
        }
        else
        {
            // Inside the bands: hold the previous direction and ratchet the active band
            // so it only ever moves in the trend's favour.
            out.supertrend[ii] = out.supertrend[ii - 1];
            if (out.supertrend[ii] > 0 && lowerband[ii] < lowerband[ii - 1])
            {
                lowerband[ii] = lowerband[ii - 1];
            }
            if (out.supertrend[ii] < 0 && upperband[ii] > upperband[ii - 1])
            {
                upperband[ii] = upperband[ii - 1];
            }
        }
    }

    return out;
}

std::vector<float> TALIB_SuperTrend_dir_only(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                                             const int atr_window, const float atr_multi)
{
    const std::vector<int> dir = TALIB_SuperTrend(high, low, close, atr_window, atr_multi).supertrend;
    return std::vector<float>(dir.begin(), dir.end());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Multi-timeframe
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Resampled RESAMPLE_TIMEFRAME(const KLINEf &kline_in, const int bars_per_group, const int ltf_minutes, const int htf_minutes)
{
    if (bars_per_group < 1)
    {
        indicator_fatal("RESAMPLE_TIMEFRAME", "bars_per_group must be >= 1, got " + std::to_string(bars_per_group));
    }
    if (ltf_minutes < 1 || htf_minutes % ltf_minutes != 0 || htf_minutes / ltf_minutes != bars_per_group)
    {
        indicator_fatal("RESAMPLE_TIMEFRAME", "bars_per_group (" + std::to_string(bars_per_group) +
                                                  ") does not match htf/ltf (" + std::to_string(htf_minutes) + "/" +
                                                  std::to_string(ltf_minutes) + ")");
    }

    Resampled result{};
    result.bars_per_group = bars_per_group;
    result.kline.name = kline_in.name;

    KLINEf &out = result.kline;
    const size_t n = kline_in.close.size();
    if (n == 0)
    {
        out.nb = 0;
        return result;
    }

    // Find the first bar that actually opens a higher-timeframe period. Slicing from
    // index 0 regardless -- what both mtf strategies did -- turns a series starting at
    // 16:55 into "hourly" candles spanning 16:55 to 17:50.
    const int64_t htf_seconds = static_cast<int64_t>(htf_minutes) * 60;
    size_t offset = 0;
    while (offset < n && kline_in.timestamp[offset] % htf_seconds != 0)
    {
        ++offset;
    }
    if (offset >= n)
    {
        indicator_fatal("RESAMPLE_TIMEFRAME", "no bar in " + kline_in.name + " starts a " +
                                                  std::to_string(htf_minutes) + "m period");
    }
    result.ltf_offset = offset;

    // A trailing partial group is dropped: an unfinished candle is not a candle.
    const size_t nb_groups = (n - offset) / static_cast<size_t>(bars_per_group);
    out.timestamp.reserve(nb_groups);
    out.open.reserve(nb_groups);
    out.high.reserve(nb_groups);
    out.low.reserve(nb_groups);
    out.close.reserve(nb_groups);
    out.volume.reserve(nb_groups);

    const bool has_volume = kline_in.volume.size() == n;

    for (size_t g = 0; g < nb_groups; ++g)
    {
        const size_t begin = offset + g * static_cast<size_t>(bars_per_group);
        const size_t end = begin + static_cast<size_t>(bars_per_group);

        float hi = kline_in.high[begin];
        float lo = kline_in.low[begin];
        double vol = 0.0;
        for (size_t i = begin; i < end; ++i)
        {
            hi = std::max(hi, kline_in.high[i]);
            lo = std::min(lo, kline_in.low[i]);
            if (has_volume)
            {
                vol += static_cast<double>(kline_in.volume[i]);
            }
        }

        out.timestamp.push_back(kline_in.timestamp[begin]);
        out.open.push_back(kline_in.open[begin]);
        out.high.push_back(hi);
        out.low.push_back(lo);
        out.close.push_back(kline_in.close[end - 1]);
        out.volume.push_back(static_cast<float>(vol));
    }

    out.nb = static_cast<uint>(out.close.size());
    return result;
}

std::vector<float> PROJECT_HTF_TO_LTF(const std::vector<float> &htf_series, const int bars_per_group,
                                      const size_t ltf_size, const size_t ltf_offset, const float fill)
{
    if (bars_per_group < 1)
    {
        indicator_fatal("PROJECT_HTF_TO_LTF", "bars_per_group must be >= 1, got " + std::to_string(bars_per_group));
    }

    const size_t group = static_cast<size_t>(bars_per_group);
    std::vector<float> out(ltf_size, fill);

    // Shift right by one full group: candle g's value is only known once candle g has
    // closed, so it becomes visible on the first bar of group g+1. Without this the
    // strategy would read a higher-timeframe close while that candle was still forming.
    for (size_t g = 0; g < htf_series.size(); ++g)
    {
        const size_t begin = ltf_offset + (g + 1) * group;
        if (begin >= ltf_size)
        {
            break;
        }
        const size_t end = std::min(begin + group, ltf_size);
        std::fill(out.begin() + static_cast<long>(begin), out.begin() + static_cast<long>(end), htf_series[g]);
    }

    return out;
}
