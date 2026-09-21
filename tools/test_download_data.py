"""Repair tests exercise assembly and persistence with supplied upstream responses and archives."""
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import zipfile
import download_data as d

BASE = 1704067200000
STEP = 3600000
ROWS = [[BASE+i*STEP,100+i,102+i,99+i,101+i,10] for i in range(6)]

def archive(rows):
    text=io.StringIO();csv.writer(text).writerows(rows)
    out=io.BytesIO()
    with zipfile.ZipFile(out,"w") as z:z.writestr(zipfile.ZipInfo("candles.csv"),text.getvalue())
    return out.getvalue()

class RepairTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.cache=patch.object(d,"CACHE_DIR",self.root/"cache");self.cache.start();self.addCleanup(self.cache.stop)

    def source(self,url,timeout=60):
        if url.endswith(".CHECKSUM"):
            return (hashlib.sha256(archive(ROWS)).hexdigest()+"  candles.zip").encode()
        if url.endswith(".zip"):
            return archive(ROWS)
        raise AssertionError("Unexpected network request: "+url)

    def test_missing_middle_is_downloaded_and_cpp_loader_reads_it(self):
        path=d.kline_path(self.root,"spot","BTC","1h")
        d.save_rows(path,{r[0]:r for r in ROWS if r[0]!=BASE+2*STEP})
        with patch.object(d,"http_get",side_effect=self.source) as network:
            result,rows=d.ensure_candles(self.root,"spot","BTC","1h",BASE,BASE+6*STEP)
        self.assertTrue(network.called);self.assertEqual(rows,ROWS)
        binary=d.ROOT/"build"/"release"/"verification_regression.exe"
        output=json.loads(subprocess.check_output([str(binary),str(result)],text=True))
        self.assertEqual(output["rows"],6);self.assertEqual(output["close_sum"],621)

    def test_conflicting_downloaded_candles_fail_before_replacement(self):
        path=d.kline_path(self.root,"spot","BTC","1h")
        d.save_rows(path,{ROWS[0][0]:ROWS[0]})
        original=path.read_bytes()
        conflicting=list(ROWS[1]);conflicting[4]+=0.5
        blob=archive(ROWS+[conflicting])
        def source(url,timeout=60):
            return (hashlib.sha256(blob).hexdigest()+" x").encode() if url.endswith(".CHECKSUM") else blob
        with patch.object(d,"http_get",side_effect=source):
            with self.assertRaisesRegex(d.DownloadError,"Conflicting duplicate"):
                d.ensure_candles(self.root,"spot","BTC","1h",BASE,BASE+6*STEP)
        self.assertEqual(path.read_bytes(),original)

    def test_rest_repairs_archive_hole(self):
        def source(url,timeout=60):
            blob=archive([r for i,r in enumerate(ROWS) if i!=3])
            if url.endswith(".CHECKSUM"):return (hashlib.sha256(blob).hexdigest()+" x").encode()
            if url.endswith(".zip"):return blob
            if "/klines?" in url:return json.dumps([ROWS[3]]).encode()
            raise AssertionError(url)
        with patch.object(d,"http_get",side_effect=source):
            _,rows=d.ensure_candles(self.root,"futures","BTC","1h",BASE,BASE+6*STEP)
        self.assertEqual(rows,ROWS)

    def test_unavailable_hole_fails_after_attempt(self):
        def source(url,timeout=60):
            if "/klines?" in url:return b"[]"
            raise FileNotFoundError(url)
        with patch.object(d,"http_get",side_effect=source) as network:
            with self.assertRaisesRegex(d.DownloadError,"unresolved gap"):
                d.ensure_candles(self.root,"spot","BTC","1h",BASE,BASE+STEP)
        self.assertTrue(network.called)

    def test_complete_data_uses_no_network_and_offline_holes_fail(self):
        path=d.kline_path(self.root,"spot","BTC","1h");d.save_rows(path,{r[0]:r for r in ROWS})
        with patch.object(d,"http_get",side_effect=AssertionError("network")):
            _,rows=d.ensure_candles(self.root,"spot","BTC","1h",BASE,BASE+6*STEP,True)
            self.assertEqual(rows,ROWS)
            with self.assertRaises(d.DownloadError):
                d.ensure_candles(self.root,"spot","BTC","1h",BASE,BASE+7*STEP,True)

    def test_corrupt_zip_is_replaced_not_trusted(self):
        path=self.root/"cache.zip";path.write_bytes(b"broken");path.with_suffix(".sha256").write_text("0"*64)
        with patch.object(d,"http_get",side_effect=self.source):
            blob=d.fetch_zip("https://archive/test.zip",path)
        self.assertEqual(d.zip_rows(blob),[[str(x) for x in r] for r in ROWS])

    def test_atomic_failure_preserves_existing_file(self):
        path=self.root/"existing";path.write_text("good")
        with patch.object(d.os,"replace",side_effect=OSError("interrupted")):
            with self.assertRaises(OSError):d.write_atomic(path,"new")
        self.assertEqual(path.read_text(),"good");self.assertFalse(list(self.root.glob("*.tmp")))

    def test_invalid_rows_and_duplicate_conflicts_rejected(self):
        bad=list(ROWS[0]);bad[2]=90
        with self.assertRaises(d.DownloadError):d.validate_row(bad,STEP)
        bad=list(ROWS[0]);bad[4]=float("nan")
        with self.assertRaises(d.DownloadError):d.validate_row(bad,STEP)
        path=self.root/"test.csv";path.write_text("date,open,high,low,close,volume\n"+
            ",".join(map(str,ROWS[0]))+"\n"+",".join(map(str,[BASE,200,202,199,201,10]))+"\n")
        with self.assertRaises(d.DownloadError):d.read_rows(path,STEP)
        micro=list(ROWS[0]);micro[0]*=1000
        self.assertEqual(d.validate_row(micro,STEP),ROWS[0])

    def test_funding_variable_intervals_and_deleted_record_repair(self):
        source=[dict(fundingTime=BASE+i*STEP,fundingRate="0.001",markPrice="110") for i in (0,4,6)]
        with patch.object(d,"fetch_zip",return_value=None), patch.object(d,"http_get",side_effect=[json.dumps(source).encode(),b"[]"]):
            path,events=d.ensure_funding(self.root,"BTC",BASE,BASE+8*STEP)
        self.assertEqual([e["timestamp_ms"] for e in events],[BASE,BASE+4*STEP,BASE+6*STEP])
        ledger=json.loads(path.read_text());ledger["events"].pop(1);path.write_text(json.dumps(ledger))
        with self.assertRaises(d.DownloadError):d.ensure_funding(self.root,"BTC",BASE,BASE+8*STEP,True)
        with patch.object(d,"fetch_zip",return_value=None), patch.object(d,"http_get",side_effect=[json.dumps(source).encode(),b"[]"]) as network:
            _,events=d.ensure_funding(self.root,"BTC",BASE,BASE+8*STEP)
        self.assertTrue(network.called);self.assertEqual(len(events),3)

    def test_funding_archive_repairs_missing_rest_event(self):
        blob=archive([[BASE,4,0.001],[BASE+4*STEP,4,0.002]])
        page=[dict(fundingTime=BASE,fundingRate="0.001",markPrice="100")]
        with patch.object(d,"fetch_zip",return_value=blob), patch.object(d,"http_get",
                side_effect=[json.dumps(page).encode(),b"[]",json.dumps([[BASE+4*STEP,110]]).encode()]):
            _,events=d.ensure_funding(self.root,"BTC",BASE,BASE+8*STEP)
        self.assertEqual(events[1],dict(timestamp_ms=BASE+4*STEP,rate=0.002,mark_price=110,interval_hours=4))

    def test_off_boundary_funding_fetches_finer_execution_candles(self):
        cfg=dict(data_dir=str(self.root), coins=["BTC"],
            strategies={"study":dict(market="futures",timeframe="4h")},
            run=dict(start="2024-01-01",end="2024-01-01T08:00:00"))
        path=self.root/"config.json";path.write_text(json.dumps(cfg))
        def source(url,timeout=60):
            if "/fundingRate/" in url:
                raise FileNotFoundError(url)
            if "/fundingRate?" in url:
                cursor=int(d.urllib.parse.parse_qs(d.urllib.parse.urlparse(url).query)["startTime"][0])
                return json.dumps([dict(fundingTime=BASE+i*STEP,fundingRate="0.001",markPrice="100")
                    for i in (1,5) if BASE+i*STEP>=cursor]).encode()
            tf=url.split("/")[-2]
            step=d.TIMEFRAME_MS[tf]
            rows=[[ts,100,101,99,100,10] for ts in range(BASE,BASE+8*STEP,step)]
            blob=archive(rows)
            return (hashlib.sha256(blob).hexdigest()+" candles.zip").encode() if url.endswith(".CHECKSUM") else blob
        with patch.object(d,"http_get",side_effect=source):
            manifest=d.ensure_strategy(path,"study")
        self.assertEqual(manifest["signal_seconds"],14400)
        self.assertEqual(manifest["execution_seconds"],3600)
        rows=d.read_rows(manifest["execution_files"][0],STEP)
        self.assertEqual(len(rows),8)
        self.assertEqual(sum(r[4] for r in rows.values()),800)

    def test_funding_missing_period_fails(self):
        with patch.object(d,"fetch_zip",return_value=None), patch.object(d,"http_get",return_value=b"[]"):
            with self.assertRaisesRegex(d.DownloadError,"unresolved interval"):
                d.ensure_funding(self.root,"BTC",BASE,BASE+24*STEP)

if __name__=="__main__":
    unittest.main()
