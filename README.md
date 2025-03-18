# Trading Strategy Backtester C++

## Purpose
A high-performance backtesting framework for running fast and memory-efficient backtests on various trading strategies (EMA crossover, SuperTrend, Bollinger Bands, etc.).

**Notes:** 
- Strategies with `F_` prefix are designed for futures markets, while those without the prefix are for spot markets.
- Uses `float` (32-bit) variables instead of `double` (64-bit), most of the time, that goes about 4 times faster and uses 4 times less memory for equivalent results.

## Features
- Leverages TA-lib (Technical Analysis Library) for technical indicators
- Extremely fast performance compared to Python implementations (~3000x faster)
- Low memory footprint
- Supports multiple trading pairs
- Includes built-in statistics for strategy evaluation (win rate, drawdown, Calmar ratio, etc.)
- Can be used as boilerplate code to test custom strategies

## Available Strategies

### EMA-Based Strategies
- `backtest_double_EMA_float.cpp`: Basic 2-EMA crossover strategy optimizer
- `backtest_double_EMA_StochRSI_float.cpp`: 2-EMA crossover with Stochastic RSI filter
- `backtest_double_EMA_StochRSI_float_muti_pair.cpp`: Multi-pair version of EMA with StochRSI strategy

### Multiple Indicator Strategies
- `BigWill.cpp`: Spot market strategy using Williams %R, EMAs, and Awesome Oscillator
- `F_BigWill.cpp`: Futures trading strategy combining Awesome Oscillator, Williams %R, and EMA crossovers
- `3EMA_SRSI_ATR.cpp`: Spot market strategy using triple EMA, Stochastic RSI, and ATR
- `F_3EMA_SRSI_ATR.cpp`: Futures trading strategy with three EMAs, Stochastic RSI, and dynamic ATR-based exit levels

### TRIX-Based Strategies
- `backtest_TRIX.cpp`: Optimizer for TRIX indicator strategies (inspired by Crypto Robot)
- `backtest_TRIX_multi_pair.cpp`: Multi-pair TRIX strategy (BTC, ETH, BNB, XRP) with position management

### Bollinger Bands Strategies
- `BBTREND.cpp`: Spot market strategy using Bollinger Bands breakouts
- `F_BBTREND.cpp`: Futures trading strategy for trend breakouts with Bollinger Bands and EMA filter

### SuperTrend Strategies
- `SuperTrend_EMA_ATR.cpp`: Combines SuperTrend indicator with EMA and ATR-based position management
- `SuperReversal_mtf.cpp`: Spot market multi-timeframe strategy using SuperTrend and EMAs
- `F_SuperReversal_mtf.cpp`: Futures multi-timeframe strategy with SuperTrend and EMA crossovers

## Installation

### Required Packages (Ubuntu/Debian)
```bash
sudo apt-get update
sudo apt-get install -y build-essential libjsoncpp-dev nlohmann-json3-dev gcc g++
```

### Setup Options

#### Option 1: Direct installation on Linux/Ubuntu
1. Install a C/C++ compiler (recommended: g++ and gcc)
2. Compile TA-lib:
   ```bash
   cd talib
   ./configure --prefix=/$(pwd)/talib_install/
   make
   make install
   make clean
   ```
3. Compile the strategies:
   ```bash
   cd ..  # return to root folder
   make
   ```

#### Option 2: Windows using WSL (Windows Subsystem for Linux)
1. Install WSL (Windows Subsystem for Linux) on your Windows machine
2. Follow the Linux/Ubuntu instructions above within your WSL environment

#### Using the install script should work too
```bash
sh install.sh
```

#### Option 3: Docker
Use the docker-compose project in the `docker_easy` directory (see `README.md` inside for details).

## Running Backtests

After compilation, execute a strategy backtest:
```bash
./backtest_double_EMA_float.exe  # or any other compiled strategy
```

## Data Format
Data is provided in the `data` directory for BTC and ETH 1h data.

If you want more data you can see on the folder `data_downloader_freqtrade` a project to download data with Freqtrade and Docker as `.json` format (See `README.md` inside).

## Performance Comparison

C++ vs Python implementation of `backtest_double_EMA_StochRSI_float`:

```
C++ Implementation:
-------------------------------------
Number of backtests performed : 9120
Time taken                    : 2 seconds 
RAM usage                     : 14.3 MB
-------------------------------------

Python Implementation:
-------------------------------------
Number of backtests performed :  9120
Time taken                    :  5955 seconds 
RAM usage                     :  86.0 MB
-------------------------------------

The C++ code is ~3000 times faster than the Python version.
```

## Example Results

Sample output from 2-EMA crossover strategy:

```
-------------------------------------
Strategy to test: 2-EMA crossover simple
DATA FILE TO PROCESS: ./data/BTCUSDT_1HOUR.txt
Initialized TA-Lib !
Initialized calculations.
-------------------------------------
Begin day      : 2017/8/17
End day        : 2022/7/6
OPEN/CLOSE FEE : 0.05 %
FUNDING FEE    : 0 %
LEVERAGE       : 1
Minimum number of trades required: 100
CAN LONG : 1 ; CAN SHORT : 0
-------------------------------------
BEST PARAMETER SET FOUND: 
-------------------------------------
Strategy : 2-EMA crossover simple
EMA      : 72 7
Best Gain: 2803.38%
Win rate : 28.0778%
max DD   : -45.9658%
Gain/DDC : 32.9545
Score    : 925.289
Number of trades: 463
-------------------------------------
Number of backtests performed : 39402
Time taken                    : 14 seconds 
RAM usage                     : 43.6 MB
-------------------------------------
```
