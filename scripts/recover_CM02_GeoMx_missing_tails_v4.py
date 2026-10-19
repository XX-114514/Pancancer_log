"""Reuse validated public range prefixes; fetch only missing small tails."""
import argparse
import concurrent.futures
import hashlib
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from prepare_CM02_panel_bridge_reference_scope_v1 import budget

def main(a):
    budget(a.root,256*1024**2)
    records=json.loads(a.manifest.read_text());old=a.root/'GeoMx_public_server_ranges_ipv4_v3/transport_parts'
    out=a.root/'GeoMx_public_server_tail_recovery_v4';assert not out.exists();out.mkdir(mode=0o700);parts=out/'transport_parts';parts.mkdir(mode=0o700)
    jobs=[];fds={};reused=0
    for i,r in enumerate(records):
        fd=os.open(out/r['file'],os.O_CREAT|os.O_RDWR,0o600);os.ftruncate(fd,r['bytes']);fds[i]=fd
        for start in range(0,r['bytes'],1024**2):
            end=min(start+1024**2,r['bytes'])-1;key=str(i)+'_'+str(start);n=0;body=old/(key+'.part');header=old/(key+'.headers')
            if body.exists() and header.exists():
                ranges=re.findall(r'Content-Range:\s*bytes\s+(\d+)-(\d+)/(\d+)',header.read_text(),re.I)
                if ranges and tuple(map(int,ranges[-1]))==(start,end,r['bytes']):
                    data=body.read_bytes();assert len(data)<=end-start+1;n=len(data)
                    assert os.pwrite(fd,data,start)==n;reused+=n
            jobs.extend((i,s,min(s+32768,end+1)-1) for s in range(start+n,end+1,32768))
    print(json.dumps(dict(phase='valid_range_prefixes_reused',bytes=reused,missing_small_ranges=len(jobs))),flush=True)
    def fetch(job):
        i,start,end=job;r=records[i];key=str(i)+'_'+str(start);body=parts/(key+'.part');header=parts/(key+'.headers')
        for attempt in range(3):
            p=subprocess.run(['curl','-4','--noproxy','*','--fail','--silent','--show-error','--connect-timeout','8','--max-time','25','--range',str(start)+'-'+str(end),'--dump-header',str(header),'--output',str(body),r['url']],capture_output=True)
            if p.returncode==0:
                ranges=re.findall(r'Content-Range:\s*bytes\s+(\d+)-(\d+)/(\d+)',header.read_text(),re.I)
                if ranges and tuple(map(int,ranges[-1]))==(start,end,r['bytes']) and body.stat().st_size==end-start+1:
                    data=body.read_bytes();assert os.pwrite(fds[i],data,start)==len(data);return len(data)
        raise RuntimeError('Missing public tail failed: '+key)
    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:added=sum(pool.map(fetch,jobs))
    for i,r in enumerate(records):
        os.fsync(fds[i]);os.close(fds[i]);h=hashlib.sha256()
        with (out/r['file']).open('rb') as f:
            for b in iter(lambda:f.read(4194304),b''):h.update(b)
        assert h.hexdigest()==r['sha256']
    proof=dict(status='ALL5_PUBLIC_GEOMX_SOURCES_EXACT_RANGES_LENGTH_AND_ACCEPTED_LOCAL_SHA_VERIFIED',finished_utc=datetime.now(timezone.utc).isoformat(),files=records,
               reused_exact_range_prefix_bytes=reused,new_missing_tail_bytes=added,
               original_failed_transport_preserved=True,producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'download_technical_public.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps({k:v for k,v in proof.items() if k!='files'}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);main(p.parse_args())
