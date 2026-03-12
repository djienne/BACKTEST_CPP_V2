default: backtest_double_EMA_float.cpp tools.cpp custom_talib_wrapper.cpp custom_talib_wrapper.hh tools.hh 
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./backtest_double_EMA_float.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./backtest_double_EMA_float.exe

debug: backtest_double_EMA_float.cpp tools.cpp custom_talib_wrapper.cpp custom_talib_wrapper.hh tools.hh  
	${CXX} -g -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./backtest_double_EMA_float.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./backtest_double_EMA_float.exe

SR_mtf :  tools.cpp custom_talib_wrapper.cpp SuperReversal_mtf.cpp custom_talib_wrapper.hh tools.hh  
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./SuperReversal_mtf.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./SuperReversal_mtf.exe

F_SR_mtf :  tools.cpp custom_talib_wrapper.cpp F_SuperReversal_mtf.cpp custom_talib_wrapper.hh tools.hh  
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./F_SuperReversal_mtf.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./F_SuperReversal_mtf.exe

trix_multi: tools.cpp custom_talib_wrapper.cpp backtest_TRIX_multi_pair.cpp custom_talib_wrapper.hh tools.hh Klinef.hh
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./backtest_TRIX_multi_pair.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./backtest_TRIX_multi_pair.exe

BigWill: tools.cpp trade_core.cpp custom_talib_wrapper.cpp BigWill.cpp custom_talib_wrapper.hh tools.hh trade_core.hh
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./trade_core.cpp ./BigWill.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./BigWill.exe

F_BigWill: tools.cpp trade_core.cpp custom_talib_wrapper.cpp F_BigWill.cpp custom_talib_wrapper.hh tools.hh trade_core.hh
	${CXX} -Ofast -ipo -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./trade_core.cpp ./F_BigWill.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./F_BigWill.exe

EMA2SOTCHRSIMULTI: tools.cpp trade_core.cpp custom_talib_wrapper.cpp backtest_double_EMA_StochRSI_float_muti_pair.cpp custom_talib_wrapper.hh tools.hh trade_core.hh
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./trade_core.cpp ./backtest_double_EMA_StochRSI_float_muti_pair.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./backtest_double_EMA_StochRSI_float_muti_pair.exe

3EMASRSIATR: tools.cpp custom_talib_wrapper.cpp ./3EMA_SRSI_ATR.cpp tools.hh custom_talib_wrapper.hh
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./3EMA_SRSI_ATR.cpp -L./talib/talib_install/lib -lta_lib -o ./3EMA_SRSI_ATR.exe

F_3EMASRSIATR: tools.cpp custom_talib_wrapper.cpp ./F_3EMA_SRSI_ATR.cpp tools.hh custom_talib_wrapper.hh
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./F_3EMA_SRSI_ATR.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./F_3EMA_SRSI_ATR.exe

STEMAATR: tools.cpp custom_talib_wrapper.cpp SuperTrend_EMA_ATR.cpp custom_talib_wrapper.hh tools.hh  
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./SuperTrend_EMA_ATR.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./SuperTrend_EMA_ATR.exe

BBTREND: tools.cpp custom_talib_wrapper.cpp BBTREND.cpp custom_talib_wrapper.hh tools.hh  
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./BBTREND.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./BBTREND.exe

F_BBTREND: tools.cpp custom_talib_wrapper.cpp F_BBTREND.cpp custom_talib_wrapper.hh tools.hh  
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./F_BBTREND.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./F_BBTREND.exe

verification: tools.cpp custom_talib_wrapper.cpp verification_regression.cpp custom_talib_wrapper.hh tools.hh Klinef.hh
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./verification_regression.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./verification_regression.exe

strategy_regression: tools.cpp trade_core.cpp custom_talib_wrapper.cpp strategy_regression.cpp custom_talib_wrapper.hh tools.hh trade_core.hh Klinef.hh
	${CXX} -Ofast -I./talib/talib_install/include/ ./custom_talib_wrapper.cpp ./tools.cpp ./trade_core.cpp ./strategy_regression.cpp -L./talib/talib_install/lib -lta_lib -lpthread -o ./strategy_regression.exe
