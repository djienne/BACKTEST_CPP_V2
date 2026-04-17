CXX ?= g++

TALIB_INC := -I./talib/talib_install/include/
TALIB_LIB := -L./talib/talib_install/lib -lta_lib -lpthread

# Using -O3 -fno-trapping-math instead of -Ofast: the backtester is numerics-sensitive
# and -Ofast enables -ffast-math which silently re-associates float operations, which
# invalidates the regression fixtures.
CXXFLAGS_OPT := -O3 -fno-trapping-math -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare

COMMON_SRC := ./custom_talib_wrapper.cpp ./tools.cpp ./trade_core.cpp
COMMON_HDR := tools.hh trade_core.hh custom_talib_wrapper.hh Klinef.hh strategy_runner.hh

STRATEGIES := \
    backtest_double_EMA_float \
    backtest_double_EMA_StochRSI_float_muti_pair \
    BigWill F_BigWill \
    BBTREND F_BBTREND \
    SuperReversal_mtf F_SuperReversal_mtf \
    SuperTrend_EMA_ATR \
    3EMA_SRSI_ATR F_3EMA_SRSI_ATR \
    backtest_TRIX_multi_pair

REGRESSION := verification_regression strategy_regression

.PHONY: default all tests clean
default: backtest_double_EMA_float.exe

all: $(addsuffix .exe,$(STRATEGIES) $(REGRESSION)) tests.exe

%.exe: %.cpp $(COMMON_SRC) $(COMMON_HDR)
	$(CXX) $(CXXFLAGS_OPT) $(TALIB_INC) $(COMMON_SRC) $< $(TALIB_LIB) -o $@

tests: tests.exe
tests.exe: tests.cpp $(COMMON_SRC) $(COMMON_HDR)
	$(CXX) $(CXXFLAGS_OPT) $(TALIB_INC) $(COMMON_SRC) tests.cpp $(TALIB_LIB) -o tests.exe

# Debug build of the default strategy.
debug: backtest_double_EMA_float.cpp $(COMMON_SRC) $(COMMON_HDR)
	$(CXX) -g -O0 $(TALIB_INC) $(COMMON_SRC) ./backtest_double_EMA_float.cpp $(TALIB_LIB) -o ./backtest_double_EMA_float.exe

# Friendly single-strategy aliases kept so existing contributor muscle memory still works.
SR_mtf: SuperReversal_mtf.exe
F_SR_mtf: F_SuperReversal_mtf.exe
trix_multi: backtest_TRIX_multi_pair.exe
BigWill: BigWill.exe
F_BigWill: F_BigWill.exe
EMA2SOTCHRSIMULTI: backtest_double_EMA_StochRSI_float_muti_pair.exe
3EMASRSIATR: 3EMA_SRSI_ATR.exe
F_3EMASRSIATR: F_3EMA_SRSI_ATR.exe
STEMAATR: SuperTrend_EMA_ATR.exe
BBTREND: BBTREND.exe
F_BBTREND: F_BBTREND.exe
verification: verification_regression.exe
strategy_regression: strategy_regression.exe

clean:
	rm -f *.exe
