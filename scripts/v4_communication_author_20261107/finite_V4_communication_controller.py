"""Finite CPU production after audited inputs/pilot; no impact on active GPU fits."""
import os, sys, json, subprocess, concurrent.futures
from pathlib import Path
from datetime import datetime, timezone

r=Path(os.environ['RUN_ROOT']);project=Path(os.environ['PROJECT_ROOT'])
envroot=Path(sys.executable).parents[3]/'envs';cpdb=envroot/'cpdb/bin/python';rbin=envroot/'scrna_r/bin/Rscript'
env=os.environ.copy();env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='',LC_ALL='C',LANG='C')
def stage(name,command,timeout=86400):
    start=datetime.now(timezone.utc).isoformat()
    with (r/'launch'/(name+'.log')).open('xb') as log:
        try:code=subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=timeout).returncode
        except subprocess.TimeoutExpired:code=124
    (r/'launch'/(name+'_terminal.json')).open('x').write(json.dumps({'stage':name,'UTC_start':start,'UTC_end':datetime.now(timezone.utc).isoformat(),'exit_code':code},indent=2))
    return code
def injected(script):
    env['STAGE_SCRIPT']=script
    return [sys.executable,'-B','-c',"import os,runpy;from pathlib import Path;runpy.run_path(str(Path(os.environ['RUN_ROOT'])/'scripts'/os.environ['STAGE_SCRIPT']),init_globals={'PROJECT_ROOT':Path(os.environ['PROJECT_ROOT'])})"]
for script in ['prepare_V4_communication.py','audit_V4_communication_inputs.py']:
    if stage(script,injected(script),21600):sys.exit(1)
for name in ['pilot','global']:
    if stage('CPDB_'+name,[str(cpdb),'-B',str(r/'scripts/run_V4_CPDB.py'),str(r),name]):sys.exit(1)
    if stage('audit_CPDB_'+name,[sys.executable,'-B',str(r/'scripts/audit_V4_CPDB.py'),str(r),name]):sys.exit(1)
# CellChat runs once on a separate CPU while CPDB CMT groups retain all cells.
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    cellchat=pool.submit(stage,'CellChat_global',[str(rbin),'--vanilla',str(r/'scripts/run_V4_CellChat.R'),str(r)])
    response=pool.submit(stage,'cytokine_response',injected('run_V4_cytokine_response.py'),21600)
    groups=[]
    members=__import__('pandas').read_csv(r/'audit/network_members.tsv',sep='\t')
    # Author code CMT-specific inference first; pooled-other comparator is separate.
    tags=['CM%02d_High'%i for i in range(1,15)]+[cm+'_Low' for cm in sorted(members.cm.unique())]
    for tag in tags:
        rc=stage('CPDB_'+tag,[str(cpdb),'-B',str(r/'scripts/run_V4_CPDB.py'),str(r),tag])
        audit_rc=stage('audit_CPDB_'+tag,[sys.executable,'-B',str(r/'scripts/audit_V4_CPDB.py'),str(r),tag]) if rc==0 else None
        groups.append({'tag':tag,'fit_exit_code':rc,'audit_exit_code':audit_rc})
    terminal={'UTC':datetime.now(timezone.utc).isoformat(),'groups':groups,'CellChat_exit_code':cellchat.result(),'cytokine_response_exit_code':response.result()}
(r/'launch/all_communication_terminal.json').open('x').write(json.dumps(terminal,indent=2))
if any(x['fit_exit_code'] or x['audit_exit_code'] for x in groups):sys.exit(1)
for script in ['summarize_V4_communication.py','audit_V4_communication_cohort.py']:
    if stage(script,injected(script),21600):sys.exit(1)
