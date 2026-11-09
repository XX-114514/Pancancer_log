"""Finite server controller. Never stops older jobs or accepts partial inputs."""
from pathlib import Path
from datetime import datetime, timezone
import os
import sys
import json
import time
import subprocess
import traceback
from resumable_public_zip_v1 import atomic_json, exclusive_lock

root = Path(sys.argv[1]).resolve()
deadline = time.monotonic() + 14 * 3600
env = os.environ.copy()
env['PYTHONPATH'] = str(root / 'pydeps') + os.pathsep + str(root / 'code')
env['MPLCONFIGDIR'] = str(root / 'mpl_cache')
env['NUMBA_CACHE_DIR'] = str(root / 'numba_cache')
env['PYTHONDONTWRITEBYTECODE'] = '1'
for name in ['OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
    env[name] = '4'


def wait_for(files):
    while time.monotonic() < deadline:
        if all((root / f).exists() for f in files):
            return
        time.sleep(30)
    raise TimeoutError('Finite 14-hour input/dependency wait ended')


def run(script):
    with (root / 'logs' / (script + '_PRIVATE.log')).open('ab', buffering=0) as log:
        result = subprocess.run([sys.executable, '-B', str(root / 'code' / script), str(root)],
                                env=env, stdout=log, stderr=log,
                                timeout=max(1, deadline - time.monotonic()))
    if result.returncode:
        raise RuntimeError(script + ' returned ' + str(result.returncode))


with exclusive_lock(root / 'controller.lock'):
    try:
        wait_for(['inputs/complete_PRIVATE.json', 'pydeps/pyarrow/__init__.py', 'pydeps/shapely/__init__.py'])
        if not (root / 'spatial_pilot_v1/complete.json').exists():
            run('ccrcc_spatial_pilot_v1.py')
        run('audit_ccrcc_pilot_v1.py')
        old = json.loads((root / 'config/baseline_expected_aggregate.json').read_text())
        new = json.loads((root / 'spatial_pilot_v1/complete.json').read_text())
        checks = ['cells_input', 'QC_pass', 'candidate_counts', 'distance', 'conditional_null']
        for key in checks:
            assert old[key] == new[key], 'Original baseline replay mismatch: ' + key
        atomic_json(root / 'baseline_replay_audit.json', {'pass': True, 'equal_fields': checks,
                                                       'UTC': datetime.now(timezone.utc).isoformat()})
        wait_for(['inputs/boundary_complete_PRIVATE.json'])
        if not (root / 'followup_v1/complete.json').exists():
            run('ccrcc_followup_v1.py')
        status = 'COMPUTED_AWAITING_INDEPENDENT_FOLLOWUP_AUDIT_AND_SCIENTIFIC_REVIEW'
    except Exception as error:
        (root / 'logs/controller_failure_PRIVATE.txt').write_text(traceback.format_exc())
        status = 'FAILED_OR_INPUT_WAIT_EXPIRED_REQUIRES_REVIEW'
    atomic_json(root / 'controller_terminal.json', {'UTC': datetime.now(timezone.utc).isoformat(),
                                                  'status': status, 'old_tasks_preserved': True})
