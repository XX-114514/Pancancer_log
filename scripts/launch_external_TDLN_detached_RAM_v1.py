"""Deploy generic code/control only and launch one finite source RAM task.

An interrupted launch is retried with the same frozen control name. Existing
roots without a matching launch receipt require inspection, never a new job.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import shlex
import subprocess


BOOTSTRAP = r'''
import sys,json,pathlib,hashlib,os,re,subprocess,datetime
assert sys.platform.startswith('linux')
os.umask(0o077)
line=sys.stdin.buffer.readline(1048577)
assert len(line)<=1048576 and line.endswith(b'\n')
bundle=json.loads(line)
plan_text=bundle['plan_text'];plan=json.loads(plan_text)
assert re.fullmatch(r'pancancer_TDLN_detached_[a-zA-Z0-9_]+',plan['remote_control_name'])
root=pathlib.Path('/dev/shm')/plan['remote_control_name']
plan_sha=hashlib.sha256(plan_text.encode()).hexdigest()
if root.exists():
 assert not root.is_symlink() and root.stat().st_uid==os.getuid()
 assert (root/'source_plan.json').is_file() and hashlib.sha256((root/'source_plan.json').read_bytes()).hexdigest()==plan_sha
 receipt=root/'source_launch_identity.json'
 if not receipt.is_file():
  print(json.dumps(dict(status='PARTIAL_LAUNCH_REQUIRES_INSPECTION_NO_RESTART')))
  raise SystemExit(2)
 print(json.dumps(dict(status='EXISTING_FROZEN_SOURCE_LAUNCH',identity=json.loads(receipt.read_text()))))
 raise SystemExit(0)
assert set(bundle['files'])==set(plan['source_files_SHA256'])
for name,text in bundle['files'].items():
 assert pathlib.Path(name).name==name and name not in ('.','..')
 assert hashlib.sha256(text.encode()).hexdigest()==plan['source_files_SHA256'][name]
root.mkdir(mode=0o700)
with (root/'source_plan.json').open('x') as f:f.write(plan_text)
for name,text in bundle['files'].items():
 with (root/name).open('x') as f:f.write(text)
child=subprocess.Popen([sys.executable,'-B',str(root/'supervise_external_TDLN_RAM_v1.py'),str(root)],
 stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL,start_new_session=True,cwd=str(root))
identity=dict(started_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
 source_supervisor_pid=child.pid,
 source_supervisor_start_ticks=int((pathlib.Path('/proc')/str(child.pid)/'stat').read_text().split()[21]),
 source_control_root=str(root),source_plan_SHA256=plan_sha,
 detached_session=True,generic_code_and_source_aggregate_files_only=True)
with (root/'source_launch_identity.json').open('x') as f:json.dump(identity,f)
print(json.dumps(dict(status='FINITE_SOURCE_RAM_SUPERVISOR_LAUNCHED',identity=identity)))
'''


def make_bundle(plan,code,plan_text):
    for name,sha in plan['code_SHA256'].items():
        assert hashlib.sha256((code/name).read_bytes()).hexdigest()==sha
    files={name:(code/name).read_text() for name in plan['source_files_SHA256'] if name!='generic_fixed_spatial.R'}
    files['generic_fixed_spatial.R']='\n'.join((code/n).read_text() for n in plan['R_components'])
    assert hashlib.sha256(files['generic_fixed_spatial.R'].encode()).hexdigest()==plan['R_code_SHA256']
    assert all(hashlib.sha256(v.encode()).hexdigest()==plan['source_files_SHA256'][n] for n,v in files.items())
    return dict(plan_text=plan_text,files=files)


def main():
    task=Path(os.environ['PANCANCER_EXTERNAL_FIXED_TASK']);code=Path(__file__).parent
    frozen=task/'private_frozen_fixed_native_workflow_v1.json';plan_text=frozen.read_text();plan=json.loads(plan_text)
    bundle=make_bundle(plan,code,plan_text)
    payload=(json.dumps(bundle)+'\n').encode();assert len(payload)<=1048576
    intent=dict(remote_control_name=plan['remote_control_name'],
        frozen_contract_SHA256=hashlib.sha256(frozen.read_bytes()).hexdigest(),
        source_bundle_SHA256=hashlib.sha256(payload).hexdigest(),
        bootstrap_SHA256=hashlib.sha256(BOOTSTRAP.encode()).hexdigest())
    target=task/'private_detached_launch_intent_v1.json'
    if target.exists():
        assert json.loads(target.read_text())==intent
    else:
        with target.open('x') as stream:json.dump(intent,stream,indent=2)
    args=['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes',
        '-o','ConnectTimeout=15','-o','ProxyCommand=none','-o','ProxyJump=none','-o','IPQoS=none',
        os.environ['PANCANCER_READONLY_SSH_TARGET'],
        shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -c '+shlex.quote(BOOTSTRAP)]
    child=subprocess.run(args,input=payload,stdout=subprocess.PIPE,stderr=subprocess.DEVNULL,timeout=60)
    assert len(child.stdout)<=65536
    result=json.loads(child.stdout)
    receipt=dict(observed_UTC=datetime.now(timezone.utc).isoformat(),actual_launch_control_exit=child.returncode,
        launch_result=result,scientific_completion_not_implied=True)
    path=task/'private_detached_launch_actual_v1.json'
    if path.exists():
        assert json.loads(path.read_text())['launch_result']['identity']==result['identity']
    else:
        with path.open('x') as stream:json.dump(receipt,stream,indent=2)
    assert child.returncode==0 and result['status'] in ('FINITE_SOURCE_RAM_SUPERVISOR_LAUNCHED','EXISTING_FROZEN_SOURCE_LAUNCH')
    print(json.dumps(dict(status=result['status'],source_started_UTC=result['identity']['started_UTC'],
                         completion_and_native_input_pending=True)))


if __name__=='__main__':
    main()
