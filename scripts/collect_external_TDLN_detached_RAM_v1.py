"""Read a single frozen source task; observation failure never restarts it."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import shlex
import subprocess


SNAPSHOT = r'''
import sys,json,pathlib,hashlib,os,re
name=sys.argv[1]
assert re.fullmatch(r'pancancer_TDLN_detached_[a-zA-Z0-9_]+',name)
root=pathlib.Path('/dev/shm')/name
assert root.is_dir() and not root.is_symlink() and root.stat().st_uid==os.getuid()
plan=(root/'source_plan.json').read_bytes()
def read(name):
 p=root/name
 return json.loads(p.read_text()) if p.exists() else None
def exact(identity,prefix='source_'):
 if identity is None:return None
 pid=identity[prefix+'pid'];p=pathlib.Path('/proc')/str(pid)/'stat'
 if not p.exists():return False
 fields=p.read_text().split()
 return int(fields[21])==identity[prefix+'start_ticks'] and fields[2]!='Z'
launch=read('source_launch_identity.json');worker=read('source_worker_identity.json')
exited=read('source_actual_exit.json');frames={};raw=b''
p=root/'source_aggregate_frames.jsonl'
if p.exists():
 assert p.stat().st_size<=2*1024**2
 raw=p.read_bytes()
 complete=raw[:raw.rfind(b'\n')+1]
 for line in complete.splitlines():
  frame=json.loads(line);frames[frame['kind']]=frame
if exited is not None:
 assert hashlib.sha256(raw).hexdigest()==exited['aggregate_frames_SHA256']
 assert hashlib.sha256(plan).hexdigest()==exited['source_plan_SHA256']
print(json.dumps(dict(source_plan_SHA256=hashlib.sha256(plan).hexdigest(),launch_identity=launch,
 exact_supervisor_live=exact(launch,'source_supervisor_'),worker_identity=worker,
 exact_worker_live=exact(worker),actual_exit=exited,latest_frames=frames)))
'''


def save_equal(path,value):
    if path.exists():
        assert json.loads(path.read_text())==value
    else:
        with path.open('x') as stream:
            json.dump(value,stream,indent=2);stream.write('\n')


def main():
    task=Path(os.environ['PANCANCER_EXTERNAL_FIXED_TASK'])
    frozen=task/'private_frozen_fixed_native_workflow_v1.json';plan=json.loads(frozen.read_text())
    args=['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes',
        '-o','ConnectTimeout=15','-o','ProxyCommand=none','-o','ProxyJump=none','-o','IPQoS=none',
        os.environ['PANCANCER_READONLY_SSH_TARGET'],
        shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -c '+shlex.quote(SNAPSHOT)+' '+shlex.quote(plan['remote_control_name'])]
    child=subprocess.run(args,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,timeout=40)
    assert child.returncode==0 and len(child.stdout)<=2*1024**2, 'OBSERVATION_FAILED_DO_NOT_RESTART'
    snapshot=json.loads(child.stdout);assert snapshot['source_plan_SHA256']==hashlib.sha256(frozen.read_bytes()).hexdigest()
    from accept_external_TDLN_fixed_spatial_v1 import reject_individual_payload
    frames=snapshot['latest_frames']
    for frame in frames.values():reject_individual_payload(frame)
    mapping={'SOURCE_FIXED_WORKFLOW_STARTED':'private_actual_source_start_v1.json',
        'COMPLETE_INPUT_IDENTITY_ACCEPTED':'private_complete_input_identity_v1.json',
        'COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED':'private_complete_source_fixed_workflow_v1.json',
        'EXTERNAL_FIXED_WORKFLOW_NOT_ACCEPTED':'private_source_not_accepted_v1.json',
        'EXTERNAL_FIXED_WORKFLOW_NATIVE_GATE_FAILURE':'private_native_gate_failure_v1.json'}
    for kind,name in mapping.items():
        if kind in frames:save_equal(task/name,frames[kind])
    if 'PUBLIC_INPUT_RAM_BYTES_PROGRESS' in frames:
        (task/'private_live_input_progress_v1.json').write_text(json.dumps(frames['PUBLIC_INPUT_RAM_BYTES_PROGRESS'])+'\n')
    observed=datetime.now(timezone.utc).isoformat()
    (task/'private_latest_detached_source_snapshot_v1.json').write_text(json.dumps(dict(observed_UTC=observed,snapshot=snapshot),indent=2)+'\n')
    exited=snapshot['actual_exit']
    if exited is None:
        live=snapshot['exact_worker_live'] is True or snapshot['exact_supervisor_live'] is True
        print(json.dumps(dict(status='SAME_SOURCE_TASK_PENDING' if live else 'SOURCE_PROCESS_ABSENT_WITHOUT_EXIT_REQUIRES_REVIEW',observed_UTC=observed,
            exact_worker_live=snapshot['exact_worker_live'],exact_supervisor_live=snapshot['exact_supervisor_live'],
            complete_input='COMPLETE_INPUT_IDENTITY_ACCEPTED' in frames)))
        return
    save_equal(task/'private_source_supervisor_actual_exit_v1.json',exited)
    if exited['effective_exit_code']!=0 or 'COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED' not in frames:
        target=task/'private_actual_failure_exit_v1.json'
        if not target.exists():save_equal(target,dict(status='FIXED_NATIVE_WORKFLOW_NOT_ACCEPTED',
            observed_UTC=observed,actual_exit_code=exited['effective_exit_code'],
            actual_worker_exit_code=exited['actual_worker_exit_code'],exit_code_scope='SOURCE_WORKER_WAIT_AND_OUTPUT_SUPERVISION'))
        print('{"status":"DETACHED_SOURCE_TERMINAL_FAILURE_REQUIRES_REVIEW"}')
        raise SystemExit(2)
    assert exited['actual_worker_exit_code']==0
    source=task/'private_complete_source_fixed_workflow_v1.json'
    target=task/'private_actual_exit_v1.json'
    if not target.exists():save_equal(target,dict(status='SOURCE_EXIT0_FIXED_WORKFLOW_AND_RAM_RELEASE_COMPLETE',
        observed_UTC=observed,actual_exit_code=0,exit_code_scope='ACTUAL_SOURCE_WORKER_WAIT_RESULT_FROM_SUPERVISOR',
        frozen_contract_SHA256=hashlib.sha256(frozen.read_bytes()).hexdigest(),
        complete_source_SHA256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_supervisor_exit_receipt_SHA256=hashlib.sha256((task/'private_source_supervisor_actual_exit_v1.json').read_bytes()).hexdigest()))
    print('{"status":"DETACHED_SOURCE_COMPLETE_ACTUAL_WORKER_EXIT0_RESULTS_PENDING_CONSUMER"}')


if __name__=='__main__':
    main()
