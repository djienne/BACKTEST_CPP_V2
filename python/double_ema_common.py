"""Independent spot EMA reference: completed signals, next-open fills, net P&L.

Prices and published EMA values are rounded to float32 like the C++ interface;
cash, fees and the EMA recurrence accumulate in float64. The optional Numba
version compiles this numerical reference, independently of the C++ executor.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np

def load_candles(path):
    rows=np.loadtxt(path,delimiter=",",skiprows=1)
    return rows[:,0].astype(np.int64)//1000, rows[:,1:5].astype(np.float32).astype(np.float64)

def compute_ema(close,period):
    out=np.zeros(len(close),dtype=np.float64)
    if period<2 or period>len(close):
        return out
    value=sum(float(x) for x in close[:period])/period
    out[period-1]=np.float32(value)
    alpha=2/(period+1)
    for i in range(period,len(close)):
        value+=(float(close[i])-value)*alpha
        out[i]=np.float32(value)
    return out

def ledger(times,ohlc,slow,fast,begin,end,fee_percent=0.1):
    cash=1000.0;quantity=0.0;cost=0.0;fees=0.0;peak=1000.0;drawdown=0.0;wins=0;trades=0
    records=np.zeros((2*(end-begin),6),dtype=np.float64);count=0
    equity=np.zeros(end-begin+1,dtype=np.float64);equity[0]=1000
    for i in range(begin,end):
        j=i-1
        enter=fast[j]>=slow[j] and fast[j-1]<=slow[j-1]
        leave=fast[j]<=slow[j] and fast[j-1]>=slow[j-1]
        if quantity>0 and leave:
            price=ohlc[i,0];fee=quantity*price*fee_percent/100
            proceeds=quantity*price-fee;pnl=proceeds-cost;cash+=proceeds;fees+=fee
            if pnl>0:wins+=1
            records[count]=np.array([times[i],2,price,quantity,fee,pnl]);count+=1;quantity=0
        if quantity==0 and enter and i+1<end:
            price=ohlc[i,0];cost=cash;fee=cash*fee_percent/100
            quantity=(cash-fee)/price;cash=0;fees+=fee;trades+=1
            records[count]=np.array([times[i],1,price,quantity,fee,0]);count+=1
        if quantity>0 and i+1==end:
            price=ohlc[i,3];fee=quantity*price*fee_percent/100
            proceeds=quantity*price-fee;pnl=proceeds-cost;cash+=proceeds;fees+=fee
            if pnl>0:wins+=1
            step=times[1]-times[0]
            records[count]=np.array([times[i]+step,3,price,quantity,fee,pnl]);count+=1;quantity=0
        value=cash+quantity*ohlc[i,3]
        peak=max(peak,value);drawdown=min(drawdown,100*(value/peak-1));equity[i-begin+1]=value
    return cash,fees,wins,trades,drawdown,records[:count],equity

_compiled=None
def evaluate(times,ohlc,slow_period,fast_period,begin=None,end=None,jit=False):
    global _compiled
    begin=max(slow_period,fast_period)+1 if begin is None else begin
    end=len(times) if end is None else end
    slow=compute_ema(ohlc[:,3],slow_period);fast=compute_ema(ohlc[:,3],fast_period)
    if begin<max(slow_period,fast_period)+1 or begin>=end:
        raise ValueError("Evaluation requires indicator prehistory and at least one bar")
    kernel=ledger
    if jit:
        if _compiled is None:
            from numba import njit
            _compiled=njit(cache=True)(ledger)
        kernel=_compiled
    cash,fees,wins,trades,dd,fills,equity=kernel(times,ohlc,slow,fast,begin,end)
    gain=100*(cash/1000-1);win=100*wins/trades if trades else 0
    ddc=100*(1/(1+dd/100)-1) if -100<dd<0 else 0
    return dict(wallet=cash,commissions=fees,trades=trades,win_rate_percent=win,max_drawdown_percent=dd,
                gain_percent=gain,score=gain/ddc*win if ddc else 0,fills=fills,equity=equity)

def main(jit=False):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datafile",type=Path,required=True)
    parser.add_argument("--ema1",type=int,required=True,help="slow period")
    parser.add_argument("--ema2",type=int,required=True,help="fast period")
    args=parser.parse_args()
    times,ohlc=load_candles(args.datafile);begin=max(args.ema1,args.ema2)+1
    split=begin+int((len(times)-begin)*0.8)
    started=time.perf_counter()
    training=evaluate(times,ohlc,args.ema1,args.ema2,begin,split,jit)
    holdout=evaluate(times,ohlc,args.ema1,args.ema2,split,len(times),jit)
    for result in (training,holdout):
        result["fills"]=result["fills"].tolist();result["equity"]=result["equity"].tolist()
    print(json.dumps(dict(training=training,holdout=holdout,elapsed_seconds=time.perf_counter()-started)))
