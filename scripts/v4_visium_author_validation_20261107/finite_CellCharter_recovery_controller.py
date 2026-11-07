"""Retry only pre-training zero-HVG failures; reuse accepted models unchanged."""
import os, sys, time, json, subprocess, concurrent.futures
from pathlib import Path
from datetime import datetime, timezone

r = Path(os.environ['RUN_ROOT'])
envroot = Path(sys.executable).parents[3] / 'envs'
cpu_python = envroot / 'spatial_cellcharter_20260905/bin/python'
env = os.environ.copy()
env.update(OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='4', MKL_NUM_THREADS='4',
           CUDA_VISIBLE_DEVICES='', LC_ALL='C', LANG='C', PYTHONNOUSERSITE='1')
start = time.monotonic()

def wait_for(path):
    while not path.exists():
        if time.monotonic() - start > 432000:
            raise TimeoutError('Five-day bounded niche recovery/cohort wait expired')
        time.sleep(30)

def stage(name, command, timeout=86400):
    begin = datetime.now(timezone.utc).isoformat()
    with (r / 'launch' / (name + '.log')).open('xb') as log:
        try:
            result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
            code = result.returncode
        except subprocess.TimeoutExpired:
            code = 124
    (r / 'launch' / (name + '_terminal.json')).open('x').write(json.dumps({
        'UTC_start': begin, 'UTC_end': datetime.now(timezone.utc).isoformat(), 'exit_code': code,
    }, indent=2))
    return code

original = json.loads((r / 'inputs/section_jobs_PRIVATE.json').read_text())
def already_failed(job):
    path = r/'launch'/('niche_fit_' + job['section_code'] + '_terminal.json')
    return 0 if path.exists() and json.loads(path.read_text())['exit_code'] else 1
original.sort(key=already_failed)

def recover(row):
    key = row['section_code']
    wait_for(r/'launch'/('niche_fit_' + key + '_terminal.json'))
    fitted = json.loads((r/'launch'/('niche_fit_' + key + '_terminal.json')).read_text())
    if fitted['exit_code']:
        status = 'CELLCHARTER_FIT_FAILED'
    else:
        wait_for(r/'launch'/('niche_audit_' + key + '_terminal.json'))
        audited = json.loads((r/'launch'/('niche_audit_' + key + '_terminal.json')).read_text())
        status = 'PASS_CELLCHARTER_FIT_AND_INDEPENDENT_AUDIT' if not audited['exit_code'] else 'CELLCHARTER_INDEPENDENT_AUDIT_FAILED'
    result = {'section_code': key, 'status': status, 'niche_branch': 'niches', 'successful_model_reused': not fitted['exit_code']}
    if status != 'CELLCHARTER_FIT_FAILED':
        return result
    log = (r / 'launch' / ('niche_fit_' + key + '.log')).read_text(errors='replace')
    if 'assert np.all(np.asarray(b.X.sum(axis=1)).ravel() > 0)' not in log:
        return result
    result.update(niche_branch='niches_retry2', successful_model_reused=False)
    if stage('niche_retry2_fit_' + key, [str(cpu_python), '-B', str(r / 'scripts/fit_CellCharter_Visium_v2.py'), str(r), key]):
        result['status'] = 'CELLCHARTER_QUALIFIED_VIEW_FIT_FAILED'
    elif stage('niche_retry2_audit_' + key, [sys.executable, '-B', str(r / 'scripts/audit_CellCharter_Visium_v2.py'), str(r), key]):
        result['status'] = 'CELLCHARTER_QUALIFIED_VIEW_AUDIT_FAILED'
    else:
        result['status'] = 'PASS_CELLCHARTER_FIT_AND_INDEPENDENT_AUDIT'
    return result

with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
    results = list(pool.map(recover, original))
(r / 'launch/all_niches_accepted_v2_terminal.json').open('x').write(json.dumps({
    'UTC': datetime.now(timezone.utc).isoformat(), 'results': results,
    'successful_original_models_reused_without_refit': True,
    'zero_HVG_spots_unassessable_only_for_niche_main_model_unchanged': True,
}, indent=2))
wait_for(r / 'launch/all_scores_terminal.json')
for script in ['summarize_CellCharter_Visium_v2.py', 'audit_CellCharter_cohort_v2.py']:
    env['STAGE_SCRIPT'] = script
    code = "import os,runpy;from pathlib import Path;runpy.run_path(str(Path(os.environ['RUN_ROOT'])/'scripts'/os.environ['STAGE_SCRIPT']),init_globals={'PROJECT_ROOT':Path(os.environ['PROJECT_ROOT'])})"
    if stage(script, [sys.executable, '-B', '-c', code], timeout=7200):
        sys.exit(1)
