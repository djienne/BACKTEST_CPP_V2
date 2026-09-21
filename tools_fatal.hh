#pragma once

#include <cstdlib>
#include <iostream>
#include <string>

// Loud abort with file/line context for invalid indicator calls and missing cache
// entries. Data/configuration errors use exceptions so the run reports a clear
// failure and exits without writing a successful result.
//
// This lives in its own header rather than tools.hh so that low-level headers
// (indicator_cache.hh, included from Klinef.hh, itself included by tools.hh) can use
// it without a circular include.
#define BACKTEST_FATAL(msg)                                                                \
    do                                                                                     \
    {                                                                                      \
        std::cerr << "FATAL " << __FILE__ << ":" << __LINE__ << " " << (msg) << std::endl; \
        std::abort();                                                                      \
    } while (0)
