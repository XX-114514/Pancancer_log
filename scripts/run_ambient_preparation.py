#!/usr/bin/env python3
"""Detached bounded orchestration, persistent heartbeats and final data audit."""
import argparse
from datetime import datetime,timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import shutil

NATIVE=['CXG_b252b015_b488_4d5c_b16e_968c13e48a2c','CXG_67b6b9ac_b841_42ce_b2d4_ed2ec12327bd',
        'CXG_843c1459_e852_429c_9339_17b5760e7716','CXG_34f5307e_7b4d_4a48_b68f_2ba844c6414b','GSE256326']
QC=['GSE294109','GSE286827','GSE292046','CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea']

def now():return datetime.now(timezone.utc).isoformat()

def main():
    p=argparse.ArgumentParser();p.add_argument('--run-root',required=True);p.add_argument('--collection',required=True)
    p.add_argument('--rscript',required=True);p.add_argument('--parallel-datasets',type=int,default=4)
    p.add_argument('--library-workers',type=int,default=4);p.add_argument('--qc-only',action='store_true')
    p.add_argument('--code-git-sha',required=True)
    p.add_argument('--storage-plan',required=True,help='Audited effective storage capacity and source peak estimates')
    args=p.parse_args();root=Path(args.run_root);collection=Path(args.collection)
    storage=json.loads(Path(args.storage_plan).read_text())
    reserve=int(storage['reserve_bytes']);peak_by_source=storage['peak_bytes_by_source']
    if shutil.disk_usage(root).free < reserve+sum(int(v) for v in peak_by_source.values()):
        raise RuntimeError('Insufficient output-volume capacity for expected peak plus reserve')
    for name in ['objects','logs','manifests','config']:(root/name).mkdir(exist_ok=True)
    lock=(root/'logs/orchestration.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (root/'manifests/execution_start.json').exists():raise RuntimeError('Existing attempt preserved; use new run root')
    queue=QC.copy()+([] if args.qc_only else NATIVE);active=[];completed=[];sequence=0
    env=os.environ.copy()
    for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']:
        env[k]='1'
    worker=Path(__file__).with_name('preprocess_doublet_ambient.py')
    rworker=Path(__file__).with_name('decontx_library_worker.R')
    audit=Path(__file__).with_name('audit_ambient_preparation.py')
    def write(name,value):
        with (root/'manifests'/name).open('x') as f:json.dump(value,f,indent=2)
    write('execution_start.json',dict(start_time=now(),pid=os.getpid(),queue=queue,
        parallel_datasets=args.parallel_datasets,library_workers=args.library_workers,
        script_sha256={s.name:hashlib.sha256(s.read_bytes()).hexdigest() for s in [worker,rworker,audit,Path(__file__)]},
        python=sys.version.split()[0],code_git_sha=args.code_git_sha,final_source_AC=False,global_integration=False))
    while queue or active:
        free=shutil.disk_usage(root).free
        needed=reserve+sum(int(peak_by_source[n]) for n in queue)+sum(int(peak_by_source[n]) for n,p,f in active)
        while queue and len(active)<args.parallel_datasets:
            needed=reserve+sum(int(peak_by_source[n]) for n in queue)+sum(int(peak_by_source[n]) for n,p,f in active)
            free=shutil.disk_usage(root).free
            if free<needed:
                print(now(),'STORAGE_PAUSE',free,needed,flush=True)
                break
            name=queue.pop(0);command=[sys.executable,'-u',str(worker),'--collection',str(collection),
                 '--output',str(root/'objects'),'--dataset',name,'--rscript',args.rscript,
                 '--rworker',str(rworker),'--workers',str(args.library_workers)]
            if name in NATIVE:
                command+=['--native-object',str(collection/'objects'/(name+'_native_counts_preserved_v1.zarr'))]
            log=(root/'logs'/(name+'.log')).open('x')
            proc=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT)
            active.append((name,proc,log));print(now(),name,'START',proc.pid,flush=True)
        for name,proc,log in active.copy():
            if proc.poll() is not None:
                record=dict(dataset=name,exit_code=proc.returncode,end_time=now())
                completed.append(record);log.close();active.remove((name,proc,log))
                print(now(),name,'END',proc.returncode,flush=True)
        free=shutil.disk_usage(root).free
        needed=reserve+sum(int(peak_by_source[n]) for n in queue)+sum(int(peak_by_source[n]) for n,p,f in active)
        state=dict(utc=now(),active=[dict(dataset=n,pid=p.pid) for n,p,f in active],
                   completed=completed,pending=queue,output_volume_free_bytes=shutil.disk_usage(root).free,
                   reserve_bytes=reserve,storage_pause=bool(queue and free<needed))
        write('heartbeat_%05d.json'%sequence,state);sequence+=1
        if queue or active:time.sleep(20)
    command=[sys.executable,'-u',str(audit),'--output-root',str(root/'objects'),'--collection',str(collection),
             '--report',str(root/'manifests/independent_preparation_audit.json')]
    with (root/'logs/independent_audit.log').open('x') as log:
        proc=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT)
        while proc.poll() is None:
            write('audit_heartbeat_%05d.json'%sequence,dict(utc=now(),pid=proc.pid,phase='full_count_audit'))
            sequence+=1;time.sleep(20)
    write('execution_terminal.json',dict(end_time=now(),datasets=completed,audit_exit_code=proc.returncode,
       status='completed' if proc.returncode==0 and all(c['exit_code']==0 for c in completed) else 'needs_repair',
       final_source_AC=False,global_integration=False))
    return 0 if proc.returncode==0 and all(c['exit_code']==0 for c in completed) else 1

if __name__=='__main__':sys.exit(main())
