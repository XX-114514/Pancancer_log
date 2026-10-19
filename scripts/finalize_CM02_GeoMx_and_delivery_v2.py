"""Sequential post-download qualification/audit/package, never scientific effects."""
import argparse
import json
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

def main(a):
    source=a.root/'GeoMx_public_server_tail_recovery_v4'
    terminal=source/'download_technical_public.json'
    deadline=time.monotonic()+900
    while not terminal.exists():
        if time.monotonic()>deadline:raise TimeoutError('Public source download has no accepted terminal')
        time.sleep(10)
    assert json.loads(terminal.read_text())['status']=='ALL5_PUBLIC_GEOMX_SOURCES_EXACT_RANGES_LENGTH_AND_ACCEPTED_LOCAL_SHA_VERIFIED'
    code=a.root/'code';producer=a.root/'GeoMx_server_probe_qualification_v1';audit=a.root/'GeoMx_server_independent_audit_v3'
    commands=[['qualify_CM02_GeoMx_public_probe_mapping_v1.py','--source',str(source),'--output',str(producer)],
              ['audit_CM02_GeoMx_probe_workbooks_v3.py','--source',str(source),'--producer',str(producer),'--output',str(audit)]]
    for command in commands:subprocess.run([sys.executable,str(code/command[0]),*command[1:]],check=True)
    target=a.root/'private_final_handoff_v1.json';assert not target.exists();shutil.copyfile(code/target.name,target)
    subprocess.run([sys.executable,str(code/'package_CM02_panel_bridge_delivery_v1.py'),'--root',str(a.root)],check=True)
    receipt=dict(status='PUBLIC_SOURCE_SERVER_REPLICA_COMPLETE_INDEPENDENT_PROBE_AUDIT_AND_PRIVATE_DELIVERY_ACCEPTED',
                 finished_utc=datetime.now(timezone.utc).isoformat(),no_new_spatial_clinical_inference=True)
    (a.root/'finalization_technical_public.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args())
