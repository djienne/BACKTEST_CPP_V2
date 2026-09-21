# ---------------------------------------------------------------------------
#  Trading Strategy Backtester -- build system
#
#    make                    default strategy (build/release/backtest_double_EMA_float.exe)
#    make all                every strategy + regression drivers + unit tests
#    make <Name>             one strategy, e.g. `make BigWill`
#    make tests              unit suite
#    make BUILD=debug all    -O0 -g
#    make BUILD=asan tests   AddressSanitizer + UndefinedBehaviorSanitizer
#    make BUILD=tsan all     ThreadSanitizer
#    make format             clang-format every source in place
#    make clean              remove binaries and objects
#
#  Strategies are discovered from *.cpp automatically -- adding a strategy needs
#  no edit to this file. Objects live under build/$(BUILD)/ with -MMD header
#  dependency tracking, so a header edit triggers exactly the right rebuilds and
#  the four build modes never clobber each other's objects.
# ---------------------------------------------------------------------------

CXX ?= g++

TALIB_DIR := talib/talib_install
TALIB_INC := -I$(TALIB_DIR)/include
# $ORIGIN-relative rpath: binaries find the in-tree TA-Lib without LD_LIBRARY_PATH,
# and keep working if the whole repo is moved.
TALIB_LIB := -L$(TALIB_DIR)/lib -Wl,-rpath,'$$ORIGIN/../../$(TALIB_DIR)/lib' -lta_lib

BUILD ?= release

# -O3 -fno-trapping-math rather than -Ofast: this backtester is numerics-sensitive and
# -Ofast implies -ffast-math, which silently re-associates float operations and
# can invalidate finite-value checks and numerical comparisons.
CXXFLAGS_release := -O3 -fno-trapping-math -DNDEBUG
CXXFLAGS_debug   := -O0 -g
CXXFLAGS_asan    := -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined
CXXFLAGS_tsan    := -O1 -g -fno-omit-frame-pointer -fsanitize=thread
LDFLAGS_release  :=
LDFLAGS_debug    :=
LDFLAGS_asan     := -fsanitize=address,undefined
LDFLAGS_tsan     := -fsanitize=thread

ifeq ($(origin CXXFLAGS_$(BUILD)),undefined)
$(error Unknown BUILD='$(BUILD)'. Use one of: release debug asan tsan)
endif

WARNINGS := -Wall -Wextra -Wshadow

REVISION := $(shell git rev-parse --short HEAD 2>/dev/null || echo unknown)$(shell git diff --quiet -- . ":(exclude)data" || echo -dirty)
CXXFLAGS := -DBACKTEST_REVISION=\"$(REVISION)\" -std=gnu++17 $(CXXFLAGS_$(BUILD)) $(WARNINGS) $(TALIB_INC) -MMD -MP
LDFLAGS  := $(LDFLAGS_$(BUILD))
LDLIBS   := $(TALIB_LIB) -lpthread

OBJDIR := build/$(BUILD)

# Shared library layer, linked into every binary.
COMMON_SRC := indicators.cpp tools.cpp trade_core.cpp config.cpp data_io.cpp
# Non-strategy executables (unit tests + numeric regression harnesses).
DRIVER_SRC := tests.cpp execution_tests.cpp verification_regression.cpp strategy_regression.cpp
# Everything else in the root is a strategy.
STRATEGY_SRC := $(filter-out $(COMMON_SRC) $(DRIVER_SRC),$(wildcard *.cpp))

COMMON_OBJ := $(addprefix $(OBJDIR)/,$(COMMON_SRC:.cpp=.o))
STRATEGIES := $(STRATEGY_SRC:.cpp=)
DRIVERS    := $(DRIVER_SRC:.cpp=)
ALL_EXE    := $(addprefix $(OBJDIR)/,$(addsuffix .exe,$(STRATEGIES) $(DRIVERS)))
ALL_OBJ    := $(addprefix $(OBJDIR)/,$(COMMON_SRC:.cpp=.o) $(STRATEGY_SRC:.cpp=.o) $(DRIVER_SRC:.cpp=.o))
.SECONDARY: $(ALL_OBJ)

.PHONY: default all clean format help $(STRATEGIES) $(DRIVERS) verification

default: $(OBJDIR)/backtest_double_EMA_float.exe

all: $(ALL_EXE)

# Named aliases build the executable in the selected build-mode directory.
$(STRATEGIES) $(DRIVERS): %: $(OBJDIR)/%.exe

# Loader verification driver.
verification: $(OBJDIR)/verification_regression.exe

$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.exe: $(OBJDIR)/%.o $(COMMON_OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(OBJDIR):
	@mkdir -p $(OBJDIR)

format:
	clang-format -i $(wildcard *.cpp) $(wildcard *.hh)

clean:
	rm -f *.exe
	rm -rf build

help:
	@echo "Strategies: $(STRATEGIES)"
	@echo "Drivers   : $(DRIVERS)"
	@echo "Build mode: BUILD=$(BUILD)  (release|debug|asan|tsan)"

-include $(wildcard $(OBJDIR)/*.d)
