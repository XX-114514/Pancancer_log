"""Bounded public Range transport, complete official MD5; no analysis selection."""
import argparse,hashlib,json,time,urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone

URL='https://zenodo.org/records/11198494/files/tx_merfish_Pt43.csv?download=1'
SIZE=263338496
MD5='e7e82a1fd863df190ef271441b2c48bd'
BLOCK=4*1024**2

def fetch(start):
 stop=min(SIZE,start+BLOCK)-1
 for attempt in range(3):
  try:
   req=urllib.request.Request(URL,headers={'Range':f'bytes={start}-{stop}','Accept-Encoding':'identity','User-Agent':'CM02-spatial-public-qualification/1'})
   with urllib.request.urlopen(req,timeout=60) as r:
    assert r.status==206 and r.headers['Content-Range']==f'bytes {start}-{stop}/{SIZE}'
    data=r.read(BLOCK+1);assert len(data)==stop-start+1
   return start,data
  except Exception:
   if attempt==2:raise
   time.sleep(1)

def main(a):
 assert not a.output.exists();a.output.mkdir(mode=0o700)
 f=a.output/'tx_merfish_Pt43.csv';completed=0
 with f.open('xb') as sink,ThreadPoolExecutor(max_workers=8) as pool:
  sink.truncate(SIZE)
  for future in as_completed([pool.submit(fetch,start) for start in range(0,SIZE,BLOCK)]):
   start,data=future.result();sink.seek(start);sink.write(data);completed+=len(data)
   if completed//(32*1024**2)!=(completed-len(data))//(32*1024**2):print(json.dumps(dict(phase='range_chunks_written',bytes=completed,total=SIZE)),flush=True)
 md5=hashlib.md5();sha=hashlib.sha256()
 with f.open('rb') as h:
  for block in iter(lambda:h.read(BLOCK),b''):md5.update(block);sha.update(block)
 assert completed==SIZE and md5.hexdigest()==MD5
 receipt=dict(status='FULL_SELECTED_PUBLIC_MERFISH_TRANSCRIPT_FILE_OFFICIAL_MD5_VERIFIED',finished_utc=datetime.now(timezone.utc).isoformat(),url=URL,source_record='https://zenodo.org/records/11198494',bytes=SIZE,md5=MD5,sha256=sha.hexdigest(),selection='Smallest of all4 published patient transcript CSVs before any expression or spatial effects',qualified_for_biological_analysis=False,reason='Complete numeric/gene/segmentation/geometry qualification remains separate')
 (a.output/'download_receipt_public.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args())
