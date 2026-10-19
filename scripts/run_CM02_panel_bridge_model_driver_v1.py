"""One bounded scientific driver, explicit source hashes and terminal exit."""
import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def main(a):
    worker=a.root/'code/fit_CM02_panel_bridge_LOSO_v1.py'
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    audit=a.root/'reference_native_audit_v2/independent_native_audit_technical_public.json'
    assert audit.exists()
    started=dict(started_utc=datetime.now(timezone.utc).isoformat(),driver_pid=os.getpid(),
                 worker_sha256=sha(worker),native_audit_sha256=sha(audit),
                 depth_contract_sha256=sha(a.root/'target_shared_depth_contract_v2.json'))
    receipt=a.root/'model_driver_start_technical_public.json';assert not receipt.exists()
    receipt.write_text(json.dumps(started,indent=2)+'\n')
    with (a.root/'model_validation_v1.log').open('xb') as f:
        result=subprocess.run([sys.executable,str(worker),'--root',str(a.root),
                               '--depth',str(a.root/'target_shared_depth_contract_v2.json')],stdout=f,stderr=subprocess.STDOUT)
    (a.root/'model_driver_terminal_technical_public.json').write_text(json.dumps(dict(
        finished_utc=datetime.now(timezone.utc).isoformat(),exit_code=result.returncode,
        start_receipt_sha256=sha(receipt),not_independent_model_acceptance=True),indent=2)+'\n')
    return result.returncode


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    sys.exit(main(p.parse_args()))
