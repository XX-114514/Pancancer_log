"""Read-only RAM measurement, then consume a separately frozen statistics contract."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
import os
import shlex
import subprocess
import time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    task=Path(os.environ['PANCANCER_SUBTYPE_STATE_TASK']);code=Path(__file__).parent
    freeze=task/'private_frozen_measurement_v1.json';plan=json.loads(freeze.read_text())
    for name,digest in plan['code_SHA256'].items():assert sha(code/name)==digest
    folder=task/'measurement_v1';folder.mkdir(exist_ok=False)
    save=lambda n,v:(folder/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    remote='import types,json\nbase=types.ModuleType("original_author_measurement")\n'
    remote+='exec('+repr((code/'NSCLC_author_state_counts_RAM_v1.py').read_text())+',base.__dict__)\n'
    remote+=(code/'NSCLC_subtype_marker_measurement_RAM_v1.py').read_text()+'\n'
    remote+='try:\n payload=base.read_source('+repr(os.environ['PANCANCER_PREFLIGHT_ROOT'])+','+repr(plan['base_measurement'])+')\n'
    remote+=' payload=augment(payload,'+repr(os.environ['PANCANCER_PREFLIGHT_ROOT'])+','+repr(plan)+')\n print(json.dumps(dict(kind="RAM_payload",value=payload),allow_nan=False),flush=True)\n'
    remote+='except Exception as error:\n print(json.dumps(dict(kind="failed",error_type=type(error).__name__)),flush=True)\n raise SystemExit(1)\n'
    args=['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','ConnectTimeout=15']
    if os.environ.get('PANCANCER_READONLY_PROXY_COMMAND'):args+=['-o','ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND']]
    args += [os.environ['PANCANCER_READONLY_SSH_TARGET'],'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 '+shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -']
    started=datetime.now(timezone.utc).isoformat();save('private_checkpoint.json',dict(stage='READONLY_SOURCE_INPUT',started_UTC=started))
    result=subprocess.run(args,input=remote,text=True,capture_output=True,timeout=240)
    assert result.returncode==0,'SUBTYPE_MEASUREMENT_SOURCE_NOT_ACCEPTED'
    message=json.loads(result.stdout);assert message['kind']=='RAM_payload';payload=message['value'];result=None;message=None
    save('private_original_parent_proof.json',payload['proof']);save('private_subtype_marker_proof.json',payload['subtype_measurement_proof'])
    save('private_terminal.json',dict(status='COMPLETE_MEASUREMENT_RAM_WAITING_STATISTICS',started_UTC=started,
        completed_UTC=datetime.now(timezone.utc).isoformat(),measurement_contract_SHA256=sha(freeze),
        subtype_proof_SHA256=sha(folder/'private_subtype_marker_proof.json'),all_individual_arrays_RAM_only=True,
        clinical_or_coordination_effects_not_computed=True))
    print(json.dumps(dict(status='SUBTYPE_MEASUREMENT_COMPLETE_RAM_WAITING_CONTRACT')),flush=True)
    path=task/'private_frozen_statistics_v1.json';deadline=time.monotonic()+21600
    while not path.exists():
        if time.monotonic()>deadline:save('private_governance_terminal.json',dict(status='WAIT_EXPIRED_RAM_RELEASED'));return
        time.sleep(2)
    contract=json.loads(path.read_text());assert contract['measurement_terminal_SHA256']==sha(folder/'private_terminal.json')
    for name,digest in contract['code_SHA256'].items():assert sha(code/name)==digest
    spec=importlib.util.spec_from_file_location('frozen_subtype_statistics',code/contract['statistics_module'])
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.run(payload,task,contract)


if __name__=='__main__':main()
