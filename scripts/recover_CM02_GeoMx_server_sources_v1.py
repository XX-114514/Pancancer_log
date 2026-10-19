"""Replicate five public GEO files; verify against accepted local source hashes."""
import argparse
import hashlib
import json
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from prepare_CM02_panel_bridge_reference_scope_v1 import budget

def main(a):
    budget(a.root,512*1024**2)
    expected=json.loads(a.manifest.read_text())
    out=a.root/'GeoMx_public_raw_recovery_server_v1';assert not out.exists();out.mkdir(mode=0o700)
    opener=urllib.request.build_opener(urllib.request.ProxyHandler({}))
    accepted=[]
    for record in expected:
        path=out/record['file'];h=hashlib.sha256();n=0
        with opener.open(record['url'],timeout=60) as response,path.open('wb') as f:
            assert int(response.headers['Content-Length'])==record['bytes']
            while True:
                b=response.read(1024**2)
                if not b:break
                f.write(b);h.update(b);n+=len(b)
        assert n==record['bytes'] and h.hexdigest()==record['sha256']
        accepted.append(record);print(json.dumps(dict(phase='public_source_verified',file=record['file'],bytes=n)),flush=True)
    proof=dict(status='ALL5_PUBLIC_GEOMX_SOURCES_LENGTH_AND_ACCEPTED_LOCAL_SHA_VERIFIED',
               finished_utc=datetime.now(timezone.utc).isoformat(),files=accepted,
               producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),no_clinical_inference=True)
    (out/'download_technical_public.json').write_text(json.dumps(proof,indent=2)+'\n')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);main(p.parse_args())
