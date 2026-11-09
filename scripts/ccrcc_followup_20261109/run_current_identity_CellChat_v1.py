"""Finite handoff: use the unchanged identity input, never wait on reannotation."""
from pathlib import Path
from datetime import datetime, timezone
import json
import os
import sys
import subprocess
import time
from resumable_public_zip_v1 import atomic_json, exclusive_lock

r = Path(sys.argv[1]).resolve(); rb = sys.argv[2]
o = r / 'current_identity_forward_v1'
o.mkdir(exist_ok=True)
deadline = time.monotonic() + 12 * 3600
with exclusive_lock(r / 'current_identity_CellChat_controller.lock'):
    while not (o / 'complete.json').exists() and time.monotonic() < deadline:
        time.sleep(30)
    if not (o / 'complete.json').exists():
        raise TimeoutError('Current-identity input wait expired; no annotation dependency')
    env = os.environ.copy()
    for name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS']:
        env[name] = '4'
    log = r / 'logs/current_identity_CellChat_PRIVATE.log'
    with log.open('ab', buffering=0) as f:
        result = subprocess.run([rb, '--vanilla', str(r / 'code/ccrcc_current_identity_spatial_CellChat_v1.R'),
                                 str(o), 'ccRCC_current_identity'], env=env, stdout=f, stderr=f,
                                timeout=max(1, deadline - time.monotonic()))
    audit = None
    if result.returncode == 0:
        with (r / 'logs/current_identity_CellChat_audit_PRIVATE.log').open('ab', buffering=0) as f:
            audit = subprocess.run([rb, '--vanilla', str(r / 'code/audit_ccrcc_current_identity_CellChat_v1.R'),
                                    str(o), 'ccRCC_current_identity'], env=env, stdout=f, stderr=f,
                                   timeout=max(1, deadline - time.monotonic())).returncode
    atomic_json(o / 'CellChat_controller_terminal.json', {'UTC': datetime.now(timezone.utc).isoformat(),
                'CellChat_exit': result.returncode, 'audit_exit': audit,
                'original_identity_changed': False, 'identity_review_prerequisite': False})
