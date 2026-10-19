"""Bounded public GEO IPv4 ranges, exact Content-Range/length/full SHA."""
import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from prepare_CM02_panel_bridge_reference_scope_v1 import budget

def main(a):
    budget(a.root,256*1024**2)
    records=json.loads(a.manifest.read_text())
    out=a.root/'GeoMx_public_server_ranges_ipv4_v3';assert not out.exists();out.mkdir(mode=0o700)
    parts=out/'transport_parts';parts.mkdir(mode=0o700);jobs=[];handles={}
    chunk=1024**2
    for i,r in enumerate(records):
        fd=os.open(out/r['file'],os.O_CREAT|os.O_RDWR,0o600);os.ftruncate(fd,r['bytes']);handles[i]=fd
        jobs.extend((i,s,min(s+chunk,r['bytes'])-1) for s in range(0,r['bytes'],chunk))
    def fetch(job):
        i,start,end=job;r=records[i];key=str(i)+'_'+str(start);body=parts/(key+'.part');header=parts/(key+'.headers')
        for attempt in range(3):
            process=subprocess.run(['curl','-4','--noproxy','*','--fail','--silent','--show-error','--connect-timeout','10','--max-time','90','--range',str(start)+'-'+str(end),'--dump-header',str(header),'--output',str(body),r['url']],capture_output=True)
            if process.returncode==0:
                h=header.read_text();matches=re.findall(r'Content-Range:\s*bytes\s+(\d+)-(\d+)/(\d+)',h,re.I)
                if matches and tuple(map(int,matches[-1]))==(start,end,r['bytes']) and body.stat().st_size==end-start+1:
                    data=body.read_bytes();assert os.pwrite(handles[i],data,start)==len(data);return len(data)
            time.sleep(attempt+1)
        raise RuntimeError('Exact public range failed: '+key)
    total=0
    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:
        for n in pool.map(fetch,jobs):
            total+=n;print(json.dumps(dict(phase='range_verified',bytes=total,total=sum(r['bytes'] for r in records))),flush=True)
    for i,r in enumerate(records):
        os.fsync(handles[i]);os.close(handles[i]);h=hashlib.sha256()
        with (out/r['file']).open('rb') as f:
            for b in iter(lambda:f.read(4194304),b''):h.update(b)
        assert h.hexdigest()==r['sha256']
    proof=dict(status='ALL5_PUBLIC_GEOMX_SOURCES_EXACT_RANGES_LENGTH_AND_ACCEPTED_LOCAL_SHA_VERIFIED',finished_utc=datetime.now(timezone.utc).isoformat(),files=records,
               transport_workers=16,producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),no_clinical_inference=True)
    (out/'download_technical_public.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps({k:v for k,v in proof.items() if k!='files'}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);main(p.parse_args())
