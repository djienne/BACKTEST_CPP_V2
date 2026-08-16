#pragma once
// Compatibility shim. The indicator library moved to indicators.hh, which now covers
// far more than TA-Lib wrappers (Keltner, Donchian, Heikin-Ashi, VWAP, timeframe
// resampling). Existing includes keep working; new code should include indicators.hh.
#include "indicators.hh"
