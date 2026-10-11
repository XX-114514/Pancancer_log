"""Bounded official Range resource recovery; not a scientific model.

Requires PANCANCER_PREFLIGHT_ROOT and PANCANCER_ANNOTATION_ROOT.
Preserves old partials. Eight concurrent 1MiB ranges, three attempts/range,
60second socket timeout, 2400second deadline/file; wholefile official MD5.
"""
from pathlib import Path
from datetime import datetime,timezone
from urllib.request import Request,build_opener,ProxyHandler
from concurrent.futures import ThreadPoolExecutor,as_completed
import os,json,time,hashlib,shutil,subprocess

def now():return datetime.now(timezone.utc).isoformat()
def read(p):return json.loads(Path(p).read_text())
def write(p,x):
    tmp=Path(str(p)+'.atomic_tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');os.replace(tmp,p)
def digest(p,name):
    h=hashlib.new(name)
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()
def get_range(item,lo,hi,folder,deadline):
    for attempt in range(1,4):
        assert time.monotonic()<deadline,'File deadline exceeded'
        # Distinct query guards against a proxy cache ignoring the Range header.
        u=item['url']+'?pancancer_range_v2='+str(lo)+'_'+str(hi)+'_'+str(attempt)
        start=now();path=folder/(str(lo).zfill(12)+'.attempt'+str(attempt)+'.bin')
        try:
            with build_opener(ProxyHandler({})).open(Request(u,headers={'Range':f'bytes={lo}-{hi}','User-Agent':'PanCancer-public-prior-audit'}),timeout=60) as r:
                assert r.status==206 and r.headers.get('Content-Range')==f'bytes {lo}-{hi}/{item["size"]}','Range response mismatch'
                n=0
                with path.open('xb') as f:
                    while True:
                        b=r.read(min(65536,hi-lo+2-n))
                        if not b:break
                        n+=len(b);assert n<=hi-lo+1 and time.monotonic()<deadline
                        f.write(b)
                assert n==hi-lo+1,'Incomplete range'
            return dict(lo=lo,hi=hi,bytes=n,attempt=attempt,chunk_file=path.name,sha256=digest(path,'sha256'),started_UTC=start,finished_UTC=now())
        except Exception as e:
            write(folder/(str(lo).zfill(12)+'.attempt'+str(attempt)+'.failure.json'),dict(error_type=type(e).__name__,started_UTC=start,finished_UTC=now(),bytes=path.stat().st_size if path.exists() else 0))
            if attempt==3:raise
            time.sleep(min(attempt*2,4))

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);O=R/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1';C=R/'fullgene_receiver_target_code_v1'
    status=R/'family_LR_fullgene_assets_v2.execution.json';dest=O/'verified_prior_inputs_v2';chunks=O/'range_chunk_attempts_v2'
    assert not status.exists() and not dest.exists() and not chunks.exists()
    regpath=R/'CM_program_screen_evidence_register_v1.json';reg=read(regpath);assert reg['schema']=='native_LR_program_evidence_register_v16' and reg['active_jobs']==[]
    failed=read(R/'family_LR_fullgene_assets_v1.execution.json');assert failed['status']=='TERMINAL_FAILED_RESOURCE_FETCH' and not Path('/proc/'+str(failed['pid'])).exists()
    assert (O/'predeclared_plan.json').exists() and C.exists();dest.mkdir();chunks.mkdir();shutil.copy2(__file__,C/'fetch_server_priors_v2.py')
    j=dict(status='RUNNING',started_UTC=now(),pid=os.getpid(),role='public_prior_resource_fetch_not_scientific_model',script_sha256=digest(__file__,'sha256'),threads=8,chunk_bytes=1048576,max_attempts=3,timeout_seconds=60,deadline_seconds_per_file=2400,transport='direct_urllib_no_env_proxy_exact_HTTP206_ranges');write(status,j)
    rootpath=AR/'control/handoffs/20261011_GSE243013_preflight_v1.json';checkpoint=R/'cross_source_cell_identity_programme_localization_v1/current_execution_checkpoint.json'
    for label,p in [('register',regpath),('root',rootpath),('checkpoint',checkpoint)]:shutil.copy2(p,O/'execution_state_snapshots'/('assets_v2_before_'+label+'.json'))
    active=[dict(phase='fullgene_assets_v2',pid=j['pid'],role=j['role'],status_manifest=str(status),observed_UTC=now())];reg['active_jobs']=active;write(regpath,reg);root=read(rootpath);root.update(execution_status='RUNNING_FULLGENE_PRIOR_ASSETS_NOT_SCIENTIFIC_MODEL',active_jobs=active,last_verified_utc=now());write(rootpath,root);write(checkpoint,dict(verified_utc=now(),latest_complete_handoff=reg['latest_verified_handoff'],execution_status=root['execution_status'],active_jobs=active,scientific_release=False,axis_validation_complete=False))
    code=1;records=[]
    try:
        subprocess.run(['python3',str(AR/'control/refresh_state.py')],check=True)
        for item in read(O/'predeclared_plan.json')['prior']['files']:
            f=dest/item['key'];start=now();old=O/'prior_inputs'/item['key']
            if old.exists() and old.stat().st_size==item['size'] and digest(old,'md5')==item['checksum'].split(':')[1]:
                shutil.copy2(old,f);mode='copied_already_verified_v1';ranges=[]
            else:
                folder=chunks/item['key'];folder.mkdir();deadline=time.monotonic()+2400;ranges=[]
                with ThreadPoolExecutor(max_workers=8) as pool:
                    jobs=[pool.submit(get_range,item,lo,min(lo+1048575,item['size']-1),folder,deadline) for lo in range(0,item['size'],1048576)]
                    for future in as_completed(jobs):
                        ranges.append(future.result())
                        progress=dict(file=item['key'],complete_chunks=len(ranges),total_chunks=len(jobs),verified_range_bytes=sum(x['bytes'] for x in ranges),observed_UTC=now(),resource_only=True)
                        write(dest/'range_progress.json',progress)
                        if len(ranges)%8==0 or len(ranges)==len(jobs):print(json.dumps(progress),flush=True)
                ranges.sort(key=lambda x:x['lo']);cursor=0
                with f.open('xb') as output:
                    for part in ranges:
                        assert part['lo']==cursor;chunk=folder/part['chunk_file'];assert chunk.stat().st_size==part['bytes'] and digest(chunk,'sha256')==part['sha256'];output.write(chunk.read_bytes());cursor=part['hi']+1
                assert cursor==item['size'];mode='exact_ranges_assembled'
                write(folder/'complete_range_receipts.json',ranges)
            assert f.stat().st_size==item['size'] and digest(f,'md5')==item['checksum'].split(':')[1]
            record=dict(file=item['key'],bytes=item['size'],official_md5=item['checksum'],observed_md5='md5:'+digest(f,'md5'),sha256=digest(f,'sha256'),started_UTC=start,finished_UTC=now(),method=mode,chunks=len(ranges));records.append(record);print(json.dumps(record),flush=True)
            write(dest/'download_receipts.json',dict(status='FETCHING' if len(records)<3 else 'ALL3_OFFICIAL_MD5_VERIFIED',records=records))
        code=0
    except Exception as e:
        j['error_type']=type(e).__name__;raise
    finally:
        j.update(status='TERMINAL_SUCCESS_RESOURCE_ONLY' if code==0 else 'TERMINAL_FAILED_RESOURCE_FETCH',finished_UTC=now(),exit_code=code);write(status,j)
        reg=read(regpath);assert len(reg['active_jobs'])==1 and reg['active_jobs'][0]['pid']==j['pid'];reg['active_jobs']=[];write(regpath,reg);root=read(rootpath);root.update(execution_status=j['status'],active_jobs=[],last_verified_utc=now());write(rootpath,root);write(checkpoint,dict(verified_utc=now(),latest_complete_handoff=reg['latest_verified_handoff'],execution_status=root['execution_status'],active_jobs=[],verified_execution={str(status):j},scientific_release=False,axis_validation_complete=False));subprocess.run(['python3',str(AR/'control/refresh_state.py')],check=True)
    print(json.dumps(j),flush=True)
if __name__=='__main__':main()
