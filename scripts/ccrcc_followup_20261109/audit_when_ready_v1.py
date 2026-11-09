"""Finite audit handoff; terminal failure is reported instead of waiting forever."""
from pathlib import Path
from datetime import datetime, timezone
import os
import sys
import time
import json
import subprocess
from resumable_public_zip_v1 import exclusive_lock, atomic_json

r = Path(sys.argv[1]).resolve()
deadline = time.monotonic() + 14 * 3600
env = os.environ.copy()
env['PYTHONPATH'] = str(r / 'code') + os.pathsep + str(r / 'pydeps')
env['PYTHONDONTWRITEBYTECODE'] = '1'
for k in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
    env[k] = '4'
with exclusive_lock(r / 'audit_controller.lock'):
    status = 'INPUT_OR_COMPUTATION_WAIT_EXPIRED'
    while time.monotonic() < deadline:
        if (r / 'followup_v1/complete.json').exists():
            with (r / 'logs/independent_audit_PRIVATE.log').open('ab', buffering=0) as log:
                result = subprocess.run([sys.executable, '-B', str(r / 'code/audit_ccrcc_followup_v1.py'), str(r)],
                                        env=env, stdout=log, stderr=log,
                                        timeout=max(1, deadline - time.monotonic()))
            status = 'AUDIT_PASS' if result.returncode == 0 else 'AUDIT_FAILED_REQUIRES_REVIEW'
            break
        p = r / 'controller_terminal.json'
        if p.exists() and json.loads(p.read_text())['status'].startswith('FAILED'):
            status = 'UPSTREAM_FAILED_REQUIRES_REVIEW'
            break
        time.sleep(30)
    atomic_json(r / 'audit_controller_terminal.json', {'UTC': datetime.now(timezone.utc).isoformat(),
                                                      'status': status})
