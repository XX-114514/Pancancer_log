"""Official resources recovered by installed aria2, retaining all previous bytes.

Fixed same human3files and wholefile MD5. No statistical models.
Requires the stopped v2 worker be absent. Force-save retains download control.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,shutil,subprocess,sys,time

def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def write(p,d):Path(p).write_text(json.dumps(d,indent=2)+'\n')
def digest(p,k='sha256'):
    h=hashlib.new(k)
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);A=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);O=R/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1';C=R/'fullgene_receiver_target_code_v1'
    manifest=R/'family_LR_fullgene_assets_v3.execution.json';D=O/'verified_prior_inputs_v3';assert not D.exists() and not manifest.exists()
    old=read(R/'family_LR_fullgene_assets_v2.execution.json');assert old['status']=='TERMINAL_STOPPED_REPEAT_INCOMPLETE_RANGE_BODIES' and not Path('/proc/'+str(old['pid'])).exists()
    regpath=R/'CM_program_screen_evidence_register_v1.json';reg=read(regpath);assert reg['schema']=='native_LR_program_evidence_register_v16' and reg['active_jobs']==[]
    assert not (O/'verified_prior_inputs_v2/ligand_target_matrix_nsga2r_final.rds').exists()
    assert shutil.which('aria2c')=='/usr/bin/aria2c';D.mkdir();shutil.copy2(__file__,C/'fetch_server_priors_aria2_v3.py')
    rootpath=A/'control/handoffs/20261011_GSE243013_preflight_v1.json';checkpoint=R/'cross_source_cell_identity_programme_localization_v1/current_execution_checkpoint.json'
    for name,p in [('register',regpath),('root',rootpath),('checkpoint',checkpoint)]:shutil.copy2(p,O/'execution_state_snapshots'/('assets_v3_before_'+name+'.json'))
    j=dict(status='RUNNING',started_UTC=now(),pid=os.getpid(),role='public_prior_resource_fetch_not_scientific_model',script_sha256=digest(__file__),download_cli='aria2_1.34.0',file_deadline_seconds=7200,max_tries=20,connections=8,connect_timeout_seconds=20,transfer_timeout_seconds=60,force_save_control=True);write(manifest,j)
    active=[dict(phase='fullgene_assets_v3',pid=j['pid'],role='resource_worker',status_manifest=str(manifest),observed_UTC=now())]
    def state(active_jobs,status):
        reg=read(regpath);reg['active_jobs']=active_jobs;write(regpath,reg);root=read(rootpath);root.update(execution_status=status,active_jobs=active_jobs,last_verified_utc=now());write(rootpath,root);write(checkpoint,dict(verified_utc=now(),latest_complete_handoff=reg['latest_verified_handoff'],execution_status=status,active_jobs=active_jobs,verified_execution={str(manifest):j},scientific_release=False,axis_validation_complete=False))
    state(active,'RUNNING_FULLGENE_ARIA2_RESOURCE_RECOVERY');code=1;child=None;records=[]
    try:
        subprocess.run(['python3',str(A/'control/refresh_state.py')],check=True)
        for item in read(O/'predeclared_plan.json')['prior']['files']:
            start=now();f=D/item['key'];previous=O/'verified_prior_inputs_v2'/item['key'];seed=0
            if previous.exists() and previous.stat().st_size==item['size'] and digest(previous,'md5')==item['checksum'].split(':')[1]:
                shutil.copy2(previous,f);method='copied_official_MD5_verified_v2'
            else:
                # Old pieces were opened only after exact206/Content-Range check.
                # Seed only a contiguous byte prefix; all old pieces stay unchanged.
                folder=O/'range_chunk_attempts_v2'/item['key'];pieces=[]
                with f.open('xb') as out:
                    for lo in range(0,item['size'],1048576):
                        length=min(1048576,item['size']-lo);candidates=[p for p in folder.glob(str(lo).zfill(12)+'.attempt*.bin') if 0<p.stat().st_size<=length]
                        if not candidates:break
                        p=max(candidates,key=lambda q:q.stat().st_size);b=p.read_bytes();out.write(b);seed+=len(b);pieces.append(dict(offset=lo,bytes=len(b),old_piece=p.name,sha256=digest(p)))
                        if len(b)<length:break
                write(D/'contiguous_seed_receipts.json',dict(bytes=seed,wholefile_not_yet_verified=True,pieces=pieces));j.update(current_file=item['key'],seed_prefix_bytes=seed);write(manifest,j)
                args=['/usr/bin/aria2c','--no-conf=true','--continue=true','--auto-file-renaming=false','--allow-overwrite=false','--file-allocation=none','--force-save=true','--check-integrity=true','--checksum='+item['checksum'].replace(':','='),'--split=8','--max-connection-per-server=8','--min-split-size=1M','--piece-length=1M','--connect-timeout=20','--timeout=60','--max-tries=20','--retry-wait=3','--summary-interval=30','--http-accept-gzip=false','--dir='+str(D),'--out='+item['key'],item['url']]
                env={k:v for k,v in os.environ.items() if k.lower() not in ['http_proxy','https_proxy','all_proxy']}
                with (R/'family_LR_fullgene_assets_v3.aria2.stdout.log').open('xb') as output:
                    child=subprocess.Popen(args,env=env,stdout=output,stderr=subprocess.STDOUT);j['child_pid']=child.pid;write(manifest,j);active.append(dict(phase='fullgene_assets_v3',pid=child.pid,role='aria2_resource_child',status_manifest=str(manifest),observed_UTC=now()));state(active,'RUNNING_FULLGENE_ARIA2_RESOURCE_RECOVERY')
                    try:rc=child.wait(timeout=7200)
                    except subprocess.TimeoutExpired:
                        child.terminate();rc=child.wait();j['deadline_hit']=True;raise
                assert rc==0, 'aria2 terminal error, preserve output';method='aria2_contiguous_seed_resume_official_checksum'
            assert f.stat().st_size==item['size'] and digest(f,'md5')==item['checksum'].split(':')[1]
            row=dict(file=item['key'],bytes=item['size'],official_md5=item['checksum'],observed_md5='md5:'+digest(f,'md5'),sha256=digest(f),started_UTC=start,finished_UTC=now(),method=method,seed_bytes=seed);records.append(row);write(D/'download_receipts.json',dict(status='FETCHING' if len(records)<3 else 'ALL3_OFFICIAL_MD5_VERIFIED',records=records));print(json.dumps(row),flush=True)
        # Frozen biological producer retains its v2 input directory. Only append
        # the previously absent, complete target file; do not overwrite v2 history.
        dest=O/'verified_prior_inputs_v2/ligand_target_matrix_nsga2r_final.rds';assert not dest.exists();shutil.copy2(D/dest.name,dest);assert digest(dest)==records[-1]['sha256'];code=0
    except Exception as e:
        j['error_type']=type(e).__name__
        if child is not None and child.poll() is None:j['child_exit_code_after_error']=child.wait()
    finally:
        j.update(status='TERMINAL_SUCCESS_RESOURCE_ONLY' if code==0 else 'TERMINAL_FAILED_ARIA2_RESOURCE_FETCH',finished_UTC=now(),exit_code=code);write(manifest,j);state([],j['status']);subprocess.run(['python3',str(A/'control/refresh_state.py')],check=True)
    print(json.dumps(j),flush=True);raise SystemExit(code)
if __name__=='__main__':main()
