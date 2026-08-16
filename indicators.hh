#pragma once
// ---------------------------------------------------------------------------
//  Technical indicator library.
//
//  Conventions every indicator here follows:
//
//   * Output length always equals input length. TA-Lib returns a compacted array
//     starting at its lookback; the adapters scatter it back to full length and
//     left-pad the warmup region with zeros, so series[i] lines up with candle i.
//
//   * The warmup length is reported, not hidden. Single-output wrappers take an
//     optional `size_t *warmup` out-parameter; multi-output results carry a
//     `warmup` field. series[i] is only a real value for i >= warmup -- reading
//     earlier than that means comparing against padding, which is how a strategy
//     silently trades on zeros. Use strategy_runner::first_tradable_index().
//
//   * Empty input gives empty output rather than aborting.
//
//  ADDING AN INDICATOR
//  -------------------
//  For anything TA-Lib already implements, it is a one-line wrapper over one of
//  the run_* adapters in indicators.cpp -- see TALIB_CCI or TALIB_ROC. Pick the
//  adapter matching the input shape (single series / OHLC / OHLCV / multi-output),
//  declare it here, and add a test in tests.cpp asserting output length, warmup
//  position, and one value with a known analytic answer.
//
//  For anything TA-Lib lacks (Keltner, Donchian, Heikin-Ashi, VWAP), compose it
//  from the wrappers above and keep the same conventions.
// ---------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <ta-lib/ta_libc.h>
#include <unordered_map>
#include <map>
#include "Klinef.hh"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Multi-output results. Each carries the warmup shared by its component series.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct SuperTrend
{
    std::vector<int> supertrend;      // +1 uptrend, -1 downtrend
    std::vector<float> final_lowerband;
    std::vector<float> final_upperband;
    size_t warmup = 0;
};

struct MACDResult
{
    std::vector<float> macd;      // fast EMA - slow EMA
    std::vector<float> signal;    // EMA of macd
    std::vector<float> histogram; // macd - signal
    size_t warmup = 0;
};

struct StochResult
{
    std::vector<float> k;
    std::vector<float> d;
    size_t warmup = 0;
};

struct BandsResult
{
    std::vector<float> upper;
    std::vector<float> middle;
    std::vector<float> lower;
    size_t warmup = 0;
};

struct AroonResult
{
    std::vector<float> down;
    std::vector<float> up;
    size_t warmup = 0;
};

struct DirectionalResult
{
    std::vector<float> adx;
    std::vector<float> plus_di;
    std::vector<float> minus_di;
    size_t warmup = 0;
};

struct HeikinAshi
{
    std::vector<float> open;
    std::vector<float> high;
    std::vector<float> low;
    std::vector<float> close;
    size_t warmup = 0;
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Price transforms. Cheap helpers the indicators below share.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> PRICE_HL2(const std::vector<float> &high, const std::vector<float> &low);
std::vector<float> PRICE_HLC3(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close);
std::vector<float> PRICE_OHLC4(const std::vector<float> &open, const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close);
HeikinAshi HEIKIN_ASHI(const std::vector<float> &open, const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Rolling extremes
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_MIN(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_MAX(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Moving averages
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_SMA(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_EMA(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_WMA(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_DEMA(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_TEMA(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_KAMA(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);
// Hull MA: WMA(2*WMA(n/2) - WMA(n), sqrt(n)). Not in TA-Lib, composed from WMA.
std::vector<float> TALIB_HMA(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Momentum / oscillators
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_RSI(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_ROC(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_MOM(const std::vector<float> &vals, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_CCI(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_WILLR(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int length, size_t *warmup = nullptr);
std::vector<float> TALIB_ULTOSC(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                                const int period1, const int period2, const int period3, size_t *warmup = nullptr);
MACDResult TALIB_MACD(const std::vector<float> &vals, const int fast_period, const int slow_period, const int signal_period);
// Classic stochastic oscillator (%K / %D), distinct from the StochRSI family below.
StochResult TALIB_STOCH(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                        const int fastk_period, const int slowk_period, const int slowd_period);
AroonResult TALIB_AROON(const std::vector<float> &high, const std::vector<float> &low, const int period);
// Awesome Oscillator: SMA(hl2, fast) - SMA(hl2, slow).
std::vector<float> TALIB_AO(const std::vector<float> &high, const std::vector<float> &low, const int fast, const int slow);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  StochRSI
//
//  Note: this family rounds the normalized RSI to 3 decimals before smoothing, a
//  quirk inherited from the original implementation. It is kept because every
//  tracked result and tuned parameter set depends on it, but it means these do not
//  match TA-Lib's own STOCHRSI exactly.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_STOCHRSI_K(const std::vector<float> &vals, const int nb_period_stoch, const int nb_period_rsi, const int k_period, const int d_period);
std::vector<float> TALIB_STOCHRSI_D(const std::vector<float> &vals, const int nb_period_stoch, const int nb_period_rsi, const int k_period, const int d_period);
std::vector<float> TALIB_STOCHRSI_not_averaged(const std::vector<float> &vals, const int nb_period_stoch, const int nb_period_rsi);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Trend strength
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_ADX(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_PLUS_DI(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_MINUS_DI(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup = nullptr);
DirectionalResult TALIB_DMI(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period);
// Parabolic SAR. acceleration/maximum are the usual 0.02 / 0.2.
std::vector<float> TALIB_SAR(const std::vector<float> &high, const std::vector<float> &low, const double acceleration, const double maximum, size_t *warmup = nullptr);
std::vector<float> TALIB_TRIX(const std::vector<float> &vals, const int trixLength, const int trixSignal);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Volatility
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_ATR(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_NATR(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_TRANGE(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close, size_t *warmup = nullptr);
std::vector<float> TALIB_STDDEV(const std::vector<float> &vals, const int period, const double nb_dev, size_t *warmup = nullptr);
BandsResult TALIB_BBANDS_R(const std::vector<float> &close, const float dev_up, const float dev_dn, const int length);
// Keltner: EMA centreline with ATR-scaled bands. Not in TA-Lib.
BandsResult KELTNER_CHANNELS(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                             const int ema_period, const int atr_period, const float atr_multiplier);
// Donchian: rolling high/low envelope with their midline. Not in TA-Lib.
BandsResult DONCHIAN_CHANNELS(const std::vector<float> &high, const std::vector<float> &low, const int period);

// Legacy out-parameter form of Bollinger Bands, kept because several strategies call
// it in a hot loop with pre-allocated destinations. Prefer TALIB_BBANDS_R in new code.
void TALIB_BBANDS(const std::vector<float> &close,
                  const float &optInNbDevUp, const float &optInNbDevDn, const int &length,
                  std::vector<float> &OUT_u, std::vector<float> &OUT_m, std::vector<float> &OUT_l);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Volume
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<float> TALIB_OBV(const std::vector<float> &close, const std::vector<float> &volume, size_t *warmup = nullptr);
std::vector<float> TALIB_MFI(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                             const std::vector<float> &volume, const int period, size_t *warmup = nullptr);
std::vector<float> TALIB_AD(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                            const std::vector<float> &volume, size_t *warmup = nullptr);
std::vector<float> TALIB_ADOSC(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                               const std::vector<float> &volume, const int fast_period, const int slow_period, size_t *warmup = nullptr);
// Rolling VWAP over `period` bars of typical price. Not in TA-Lib. This is a moving
// window, not the session-anchored VWAP -- a backtest over years has no session anchor.
std::vector<float> VWAP_ROLLING(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                                const std::vector<float> &volume, const int period, size_t *warmup = nullptr);
// Volume relative to its own moving average: 1.0 means average volume.
std::vector<float> RELATIVE_VOLUME(const std::vector<float> &volume, const int period, size_t *warmup = nullptr);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  SuperTrend
//
//  atr_multi is a float so the usual fractional multipliers (1.5, 2.5, 3.0) can be
//  swept; it used to be an int, which excluded most real settings.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

SuperTrend TALIB_SuperTrend(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                            const int atr_window, const float atr_multi);
std::vector<float> TALIB_SuperTrend_dir_only(const std::vector<float> &high, const std::vector<float> &low, const std::vector<float> &close,
                                             const int atr_window, const float atr_multi);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Multi-timeframe
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Resampled
{
    KLINEf kline;           // the aggregated higher-timeframe candles
    size_t ltf_offset = 0;  // source index where the first whole higher candle begins
    int bars_per_group = 1;
};

// Aggregates OHLCV candles up to a higher timeframe. `bars_per_group` low-timeframe
// candles become one high-timeframe candle (12 for 5m -> 1h).
//
// Aggregation starts at the first bar that actually opens a higher-timeframe period,
// and that index is reported as `ltf_offset`. This matters: the bundled 5m futures data
// starts at 16:55, so slicing blindly from index 0 -- which both mtf strategies did --
// produced "hourly" candles spanning 16:55 to 17:50. A trailing partial group is
// dropped, since an unfinished candle is not a candle.
Resampled RESAMPLE_TIMEFRAME(const KLINEf &kline_in, const int bars_per_group, const int ltf_minutes, const int htf_minutes);

// Projects a high-timeframe series back onto the low-timeframe index, repeating each
// value `bars_per_group` times and shifting right by one full group.
//
// The shift is what prevents lookahead bias: a high-timeframe candle's value is only
// known once that candle closes, so it must not be visible to the bars that formed it.
// `ltf_offset` must be the value RESAMPLE_TIMEFRAME reported, so the projection lands
// on the same bars the aggregation came from. `fill` is written into the leading bars
// that have no completed higher candle yet.
//
// Neither trailing parameter is defaulted, deliberately: with a default on ltf_offset,
// a caller writing PROJECT_HTF_TO_LTF(s, 12, n, -777.0f) silently binds the fill value
// to the size_t offset and compiles clean.
std::vector<float> PROJECT_HTF_TO_LTF(const std::vector<float> &htf_series, const int bars_per_group,
                                      const size_t ltf_size, const size_t ltf_offset, const float fill);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
