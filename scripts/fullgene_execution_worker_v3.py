"""Verified execution phases; computational preparation never releases an axis.

PANCANCER_PREFLIGHT_ROOT=... PANCANCER_ANNOTATION_ROOT=... python worker.py production|audit
Requires completed resource phase and a frozen code manifest; no automatic acceptance.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,sys,json,hashlib,shutil,subprocess

def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):
    tmp=Path(str(p)+'.atomic_tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');os.replace(tmp,p)
def proc_present(pid):return Path('/proc/'+str(pid)).exists()

def main():
    phase=sys.argv[1];assert phase in ['production','audit']
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);O=R/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1';C=R/'fullgene_receiver_target_code_v1'
    freeze=read(C/'fullgene_code_freeze_v4.json');assert sha(O/'predeclared_plan.json')==freeze['plan_sha256']
    for name,value in freeze['code_shas'].items():assert sha(C/name)==value
    assert sha(R/'immune_same_name_family_source_review_v1/supported_LR_native_cell_panel_v1/inventory.json')==freeze['accepted_inputs']['accepted_native_inventory_sha256']
    assert sha(C/'fullgene_offline_resource_recovery_v4_spec.json')==freeze['resource_contract_sha256']
    assets=read(R/'family_LR_fullgene_assets_v4.execution.json');assert assets['status']=='TERMINAL_SUCCESS_RESOURCE_ONLY' and assets['exit_code']==0 and not proc_present(assets['pid'])
    assert assets['script_sha256']==freeze['code_shas']['verify_transferred_fullgene_priors_v4.py'] and assets['contract_sha256']==freeze['resource_contract_sha256']
    receipts=read(O/'verified_prior_inputs_v4/download_receipts.json');assert receipts['status']=='ALL3_OFFICIAL_MD5_VERIFIED'
    regpath=R/'CM_program_screen_evidence_register_v1.json';reg=read(regpath);assert reg['schema']=='native_LR_program_evidence_register_v16' and reg['active_jobs']==[]
    assert sha(reg['latest_verified_handoff'])==read(O/'predeclared_plan.json')['parent_handoff_sha256']
    if phase=='production':assert not (O/'inventory.json').exists() and not (O/'private_source_unit_arrays').exists()
    else:
        previous=read(R/'family_LR_fullgene_production_v1.execution.json');assert previous['status']=='TERMINAL_SUCCESS_PENDING_FULLGENE_AUDIT' and previous['exit_code']==0
        assert not proc_present(previous['pid']) and not proc_present(previous['child_pid']);assert (O/'inventory.json').exists() and not (O/'independent_audit.json').exists()
    script=C/('fullgene_receiver_target_readiness_v3.py' if phase=='production' else 'audit_fullgene_receiver_target_readiness_v3.py')
    status=R/('family_LR_fullgene_'+phase+'_v1.execution.json');assert not status.exists();log=R/('family_LR_fullgene_'+phase+'_v1.stdout.log');assert not log.exists()
    rootpath=AR/'control/handoffs/20261011_GSE243013_preflight_v1.json';checkpoint=R/'cross_source_cell_identity_programme_localization_v1/current_execution_checkpoint.json'
    for label,p in [('register',regpath),('root',rootpath),('checkpoint',checkpoint)]:shutil.copy2(p,O/'execution_state_snapshots'/(phase+'_before_'+label+'.json'))
    j=dict(status='RUNNING',phase=phase,started_UTC=now(),pid=os.getpid(),role='fullgene_preparation_or_computational_audit_not_axis_validation',worker_sha256=sha(__file__),script_sha256=sha(script));write(status,j)
    active=[dict(phase='fullgene_'+phase,pid=j['pid'],role='worker',status_manifest=str(status),observed_UTC=now())]
    def state(execution,active_jobs):
        d=read(regpath);d['active_jobs']=active_jobs;write(regpath,d);root=read(rootpath);root.update(execution_status=execution,active_jobs=active_jobs,last_verified_utc=now());write(rootpath,root);write(checkpoint,dict(verified_utc=now(),latest_complete_handoff=d['latest_verified_handoff'],execution_status=execution,active_jobs=active_jobs,verified_execution={str(status):j},scientific_release=False,axis_validation_complete=False))
    state('RUNNING_FULLGENE_'+phase.upper(),active);code=1;child=None
    try:
        subprocess.run(['python3',str(AR/'control/refresh_state.py')],check=True)
        env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',MKL_NUM_THREADS='4')
        with log.open('xb') as output:
            child=subprocess.Popen([sys.executable,str(script)],env=env,stdout=output,stderr=subprocess.STDOUT)
            j['child_pid']=child.pid;write(status,j);active.append(dict(phase='fullgene_'+phase,pid=child.pid,role='scientific_preparation_child' if phase=='production' else 'independent_auditor_child',status_manifest=str(status),observed_UTC=now()));state('RUNNING_FULLGENE_'+phase.upper(),active);code=child.wait()
    except Exception as e:
        j['error_type']=type(e).__name__;code=1
        if child is not None:
            # Do not clear a live child from execution state if publication fails.
            j['pipeline_exit_code']=child.wait()
    finally:
        j.update(status=('TERMINAL_SUCCESS_PENDING_FULLGENE_AUDIT' if phase=='production' else 'TERMINAL_SUCCESS_INDEPENDENT_FULLGENE_AUDIT') if code==0 else 'TERMINAL_FAILED_'+phase.upper(),finished_UTC=now(),exit_code=code);write(status,j)
        current=read(regpath);assert any(a['pid']==j['pid'] for a in current['active_jobs']);state(j['status'],[])
        try:subprocess.run(['python3',str(AR/'control/refresh_state.py')],check=True)
        except Exception as e:
            j.update(status='TERMINAL_STATE_REFRESH_FAILED',pipeline_exit_code=code,exit_code=1,state_error_type=type(e).__name__);write(status,j);state(j['status'],[]);code=1
    print(json.dumps(j),flush=True);raise SystemExit(code)
if __name__=='__main__':main()
