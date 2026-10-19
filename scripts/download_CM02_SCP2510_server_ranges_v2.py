"""Bounded parallel GCS ranges; every complete file must match accepted MD5/SHA."""
import argparse
import hashlib
import json
import os
import re
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

from prepare_CM02_panel_bridge_reference_scope_v1 import budget


def main(a):
    budget(a.root,512*1024**2)
    metadata=a.root/'server_public_source_metadata_v1'
    page=(metadata/'page.html').read_text();match=re.search(r'window\.SCP\.readOnlyTokenObject\s*=\s*',page)
    assert match
    token=json.JSONDecoder().raw_decode(page[match.end():])[0]['access_token']
    explore=json.loads((metadata/'explore.json').read_text())
    accepted=json.loads(a.accepted.read_text());inventory=json.loads((metadata/'inventory.json').read_text())
    official={f['bucket_location']:f for f in inventory['study_files']}
    assert inventory['public'] and explore['isPublic']
    out=a.root/'SCP2510_author_processed_server_ranges_v2';assert not out.exists();out.mkdir(mode=0o700)
    files=sorted(accepted['files'],key=lambda r:r['bytes']);handles={};jobs=[];block=1024**2
    for row in files:
        name=row['relative_path'];assert name in official and official[name]['upload_file_size']==row['bytes']
        dest=out/name;dest.parent.mkdir(parents=True,exist_ok=True,mode=0o700)
        fd=os.open(dest,os.O_CREAT|os.O_EXCL|os.O_RDWR,0o600);os.ftruncate(fd,row['bytes']);handles[name]=fd
        url='https://storage.googleapis.com/download/storage/v1/b/'+urllib.parse.quote(explore['bucketId'],safe='')+'/o/'+urllib.parse.quote(name,safe='')+'?alt=media'
        for start in range(0,row['bytes'],block):jobs.append((name,url,start,min(start+block,row['bytes'])-1,row['bytes']))

    def fetch(job):
        name,url,start,stop,size=job
        for attempt in range(3):
            try:
                request=urllib.request.Request(url,headers={'Authorization':'Bearer '+token,'Range':f'bytes={start}-{stop}',
                                                           'Accept-Encoding':'identity','User-Agent':'CM02-public-bounded-range/2'})
                with urllib.request.urlopen(request,timeout=60) as response:
                    assert response.status==206 or (response.status==200 and start==0 and stop+1==size)
                    if response.status==206:assert response.headers['Content-Range']==f'bytes {start}-{stop}/{size}'
                    data=response.read(stop-start+2)
                assert len(data)==stop-start+1
                assert os.pwrite(handles[name],data,start)==len(data)
                return len(data)
            except Exception:
                if attempt==2:raise
                time.sleep(1)
    done=0
    try:
        with ThreadPoolExecutor(max_workers=16) as pool:
            for future in as_completed([pool.submit(fetch,j) for j in jobs]):
                before=done;done+=future.result()
                if done//(8*1024**2)!=before//(8*1024**2):print(json.dumps(dict(phase='verified_range_bytes_written',bytes=done,total=sum(f['bytes'] for f in files))),flush=True)
    finally:
        for fd in handles.values():os.close(fd)
    receipts=[]
    for row in files:
        md=hashlib.md5();h=hashlib.sha256()
        with (out/row['relative_path']).open('rb') as f:
            for part in iter(lambda:f.read(4*1024**2),b''):md.update(part);h.update(part)
        assert md.hexdigest()==row['official_md5'] and h.hexdigest()==row['sha256']
        receipts.append(row)
    receipt=dict(status='ALL8_SERVER_PUBLIC_FILES_LENGTH_RANGE_MD5_AND_ACCEPTED_LOCAL_SHA_VERIFIED',
                 finished_utc=datetime.now(timezone.utc).isoformat(),total_bytes=done,files=receipts,
                 anonymous_public_read_token_not_recorded=True,partial_older_transport_not_accepted=True,
                 transport_workers=16,range_chunk_bytes=block)
    (out/'download_technical_public.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='files'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--accepted',type=Path,required=True);main(p.parse_args())
