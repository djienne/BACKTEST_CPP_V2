#pragma once

#include <cstdlib>
#include <iostream>
#include <string>

// Loud abort with file/line context. Use at unrecoverable failure points in loaders,
// alignment and cache lookups -- it preserves the "crash hard" semantics the
// backtester relies on, but surfaces the reason instead of a bare exit code.
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
