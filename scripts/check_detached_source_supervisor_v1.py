"""Artificial-only launch/collect lifecycle checks on the source Linux runtime.

Launch and collection are separate commands. The successful delayed worker must
finish after the first controlling SSH process has already closed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys
import uuid
from launch_external_TDLN_detached_RAM_v1 import BOOTSTRAP


def launch(directory):
    directory=Path(directory)
    sources={n:(directory/n).read_text() for n in (
        'supervise_external_TDLN_RAM_v1.py','accept_external_TDLN_fixed_spatial_v1.py')}
    jobs=[]
    cases={
        'success_after_control_exit':('import time;time.sleep(12)\n'+
            'print(json.dumps(dict(kind="COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED",ARTIFICIAL_ONLY=True)),flush=True)\n',40),
        'actual_worker_exit7':('raise SystemExit(7)\n',40),
        'individual_frame_rejected':('print(json.dumps(dict(kind="COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED",patient_ids=["ARTIFICIAL_CANARY"])),flush=True)\n',40),
        'deadline_terminates_worker':('import time;time.sleep(30)\n',2),
        'orphan_group_cleanup':('import subprocess,pathlib\n'+
            'z=subprocess.Popen([sys.executable,"-c","import time;time.sleep(60)"])\n'+
            'pathlib.Path(sys.argv[1],"ARTIFICIAL_descendant_pid.json").write_text(json.dumps(dict(pid=z.pid)))\n'+
            'print(json.dumps(dict(kind="COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED",ARTIFICIAL_ONLY=True)),flush=True)\n',40),
    }
    for label,(body,budget) in cases.items():
        name='pancancer_TDLN_detached_ARTIFICIAL_'+uuid.uuid4().hex
        worker='import sys,json,os\nprint(json.dumps(dict(kind="SOURCE_FIXED_WORKFLOW_STARTED",ARTIFICIAL_ONLY=True)),flush=True)\n'+body
        files=dict(sources,ARTIFICIAL_worker=worker)
        files['ARTIFICIAL_worker.py']=files.pop('ARTIFICIAL_worker')
        plan=dict(remote_control_name=name,source_files_SHA256={n:hashlib.sha256(s.encode()).hexdigest() for n,s in files.items()},
            source_worker_file='ARTIFICIAL_worker.py',maximum_download_seconds=budget-1,
            maximum_R_seconds=1,supervisor_grace_seconds=0,
            supervisor_address_cap_bytes=256*1024**2,inherited_address_hard_cap_bytes=2*1024**3,
            ARTIFICIAL_ONLY=True)
        bundle=dict(plan_text=json.dumps(plan),files=files)
        payload=(json.dumps(bundle)+'\n').encode()
        first=subprocess.run([sys.executable,'-B','-c',BOOTSTRAP],input=payload,
            stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=20,check=True)
        receipt=json.loads(first.stdout)
        assert receipt['status']=='FINITE_SOURCE_RAM_SUPERVISOR_LAUNCHED'
        # Same exact frozen launch is idempotent, including when already done.
        again=subprocess.run([sys.executable,'-B','-c',BOOTSTRAP],input=payload,
            stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=20,check=True)
        repeated=json.loads(again.stdout)
        assert repeated['status']=='EXISTING_FROZEN_SOURCE_LAUNCH' and repeated['identity']==receipt['identity']
        jobs.append(dict(label=label,root=str(Path('/dev/shm')/name),launch=receipt['identity']))
    output=directory/'ARTIFICIAL_detached_jobs.json'
    with output.open('x') as stream:json.dump(jobs,stream)
    print(json.dumps(dict(status='ARTIFICIAL_DETACHED_JOBS_LAUNCHED_ONCE',
                         registry=str(output),jobs=len(jobs),no_real_native_or_clinical_data=True)))


def collect(registry):
    jobs=json.loads(Path(registry).read_text());checks=[];pending=[]
    for job in jobs:
        root=Path(job['root']);path=root/'source_actual_exit.json'
        if not path.exists():
            pid=job['launch']['source_supervisor_pid'];proc=Path('/proc')/str(pid)/'stat'
            assert proc.exists() and int(proc.read_text().split()[21])==job['launch']['source_supervisor_start_ticks']
            pending.append(job['label']);continue
        result=json.loads(path.read_text());raw=(root/'source_aggregate_frames.jsonl').read_bytes()
        assert hashlib.sha256(raw).hexdigest()==result['aggregate_frames_SHA256']
        assert hashlib.sha256((root/'source_plan.json').read_bytes()).hexdigest()==result['source_plan_SHA256']
        label=job['label']
        if label in ('success_after_control_exit','orphan_group_cleanup'):
            assert result['effective_exit_code']==result['actual_worker_exit_code']==0
        elif label=='actual_worker_exit7':
            assert result['actual_worker_exit_code']==result['effective_exit_code']==7
        elif label=='individual_frame_rejected':
            assert result['effective_exit_code']==2 and b'ARTIFICIAL_CANARY' not in raw
        elif label=='deadline_terminates_worker':
            assert result['effective_exit_code']==2 and result['actual_worker_exit_code']!=0
        if label=='orphan_group_cleanup':
            pid=json.loads((root/'ARTIFICIAL_descendant_pid.json').read_text())['pid']
            proc=Path('/proc')/str(pid)/'stat'
            assert not proc.exists() or proc.read_text().split()[2]=='Z'
        checks.append(dict(label=label,actual_worker_exit=result['actual_worker_exit_code'],
                           effective_exit=result['effective_exit_code'],completed_UTC=result['observed_UTC']))
    answer=dict(status='PASS_ARTIFICIAL_DETACHED_SOURCE_LIFECYCLE' if not pending else 'SAME_ARTIFICIAL_JOBS_PENDING',
        checks=checks,pending=pending,no_real_native_or_clinical_data=True,
        no_fixture_or_history_deleted=True)
    print(json.dumps(answer))
    return answer


if __name__=='__main__':
    parser=argparse.ArgumentParser();group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--launch');group.add_argument('--collect');args=parser.parse_args()
    launch(args.launch) if args.launch else collect(args.collect)
