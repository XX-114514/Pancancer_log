"""Bounded CPU niche queue; no changes to running cell2location models."""
import os, sys, time, json, subprocess, concurrent.futures
from pathlib import Path
from datetime import datetime, timezone

r = Path(os.environ['RUN_ROOT'])
project = Path(os.environ['PROJECT_ROOT'])
envroot = Path(sys.executable).parents[3] / 'envs'
cpu_python = envroot / 'spatial_cellcharter_20260905/bin/python'
env = os.environ.copy()
env.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4',
           CUDA_VISIBLE_DEVICES='', LC_ALL='C', LANG='C', PYTHONNOUSERSITE='1')
jobs = json.loads((r / 'inputs/section_jobs_PRIVATE.json').read_text())

def stage(name, command, timeout=86400):
    begin = datetime.now(timezone.utc).isoformat()
    with (r / 'launch' / (name + '.log')).open('xb') as log:
        try:
            p = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
            exit_code = p.returncode
        except subprocess.TimeoutExpired:
            exit_code = 124
    (r / 'launch' / (name + '_terminal.json')).open('x').write(json.dumps({
        'UTC_start': begin, 'UTC_end': datetime.now(timezone.utc).isoformat(), 'exit_code': exit_code,
    }, indent=2))
    return exit_code

def fit(job):
    key = job['section_code']
    if stage('niche_fit_' + key, [str(cpu_python), '-B', str(r / 'scripts/fit_CellCharter_Visium.py'), str(r), key]):
        return {'section_code': key, 'status': 'CELLCHARTER_FIT_FAILED'}
    if stage('niche_audit_' + key, [sys.executable, '-B', str(r / 'scripts/audit_CellCharter_Visium.py'), str(r), key]):
        return {'section_code': key, 'status': 'CELLCHARTER_INDEPENDENT_AUDIT_FAILED'}
    return {'section_code': key, 'status': 'PASS_CELLCHARTER_FIT_AND_INDEPENDENT_AUDIT'}

with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    results = list(pool.map(fit, jobs))
(r / 'launch/all_niches_terminal.json').open('x').write(json.dumps({
    'UTC': datetime.now(timezone.utc).isoformat(), 'results': results,
}, indent=2))
start = time.monotonic()
while not (r / 'launch/all_scores_terminal.json').exists():
    if time.monotonic() - start > 432000:
        raise TimeoutError('Five-day bounded wait for real CM scoring')
    time.sleep(30)
for script in ['summarize_CellCharter_Visium.py', 'audit_CellCharter_cohort.py']:
    code = "import os,runpy;from pathlib import Path;runpy.run_path(str(Path(os.environ['RUN_ROOT'])/'scripts'/os.environ['STAGE_SCRIPT']),init_globals={'PROJECT_ROOT':Path(os.environ['PROJECT_ROOT'])})"
    env['STAGE_SCRIPT'] = script
    if stage(script, [sys.executable, '-B', '-c', code], timeout=7200):
        sys.exit(1)
