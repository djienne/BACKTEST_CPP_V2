"""Compare actual C++ EMA trades with an independent numerical ledger."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/"python"))
from double_ema_common import load_candles,evaluate

class ReferenceTests(unittest.TestCase):
    def test_production_ema_matches_independent_trade_ledger(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"candles.csv"
            with path.open("w",newline="") as stream:
                writer=csv.writer(stream);writer.writerow(["date","open","high","low","close","volume"])
                previous=100
                for i in range(2000):
                    close=100+12*math.sin(i*0.12)+0.01*i
                    opening=previous+0.4*math.sin(i*0.31)
                    writer.writerow([1704067200000+i*3600000,opening,max(opening,close)+1,min(opening,close)-1,close,1000])
                    previous=close
            times,ohlc=load_candles(path)
            ref=evaluate(times,ohlc,40,10)
            compiled=evaluate(times,ohlc,40,10,jit=True)
            actual=json.loads(subprocess.check_output([str(ROOT/"build/release/strategy_regression.exe"),str(path),"40","10"],text=True))
            for key in ("wallet","commissions","trades","win_rate_percent","max_drawdown_percent","score"):
                self.assertAlmostEqual(actual[key],ref[key],places=7,msg=key)
                self.assertAlmostEqual(compiled[key],ref[key],places=7,msg=key)
            code={"entry":1,"signal_exit":2,"final_exit":3}
            fills=np.array([[f["timestamp"],code[f["action"]],f["price"],f["quantity"],f["commission"],f["net_pnl"]] for f in actual["fills"]])
            self.assertGreater(len(fills),20)
            np.testing.assert_allclose(fills,ref["fills"],rtol=1e-12,atol=1e-8)
            np.testing.assert_allclose(actual["equity"],ref["equity"],rtol=1e-12,atol=1e-8)

if __name__=="__main__":
    unittest.main()
