"""Hold source-specific author subtype counts in RAM for a frozen next question."""
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import shlex
import subprocess
import time


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    task=Path(os.environ['PANCANCER_AUTHOR_STATE_TASK']);code=Path(__file__).parent
    freeze=task/'private_frozen_author_state_measurement_v1.json';plan=json.loads(freeze.read_text())
    for name,digest in plan['code_SHA256'].items():assert sha(code/name)==digest
    out=task/'measurement_v1';out.mkdir(exist_ok=False)
    save=lambda name,value:(out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    source=(code/'NSCLC_author_state_counts_RAM_v1.py').read_text()+'\n'
    source+='try:\n payload=read_source('+repr(os.environ['PANCANCER_PREFLIGHT_ROOT'])+','+repr(plan)+')\n print(json.dumps(dict(kind="RAM_payload",value=payload),allow_nan=False),flush=True)\n'
    source+='except Exception as error:\n print(json.dumps(dict(kind="failed",error_type=type(error).__name__)),flush=True)\n raise SystemExit(1)\n'
    args=['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','ConnectTimeout=15']
    if os.environ.get('PANCANCER_READONLY_PROXY_COMMAND'):args+=['-o','ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND']]
    args += [os.environ['PANCANCER_READONLY_SSH_TARGET'],
        'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 '+shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -']
    started=datetime.now(timezone.utc).isoformat()
    save('private_checkpoint.json',dict(stage='READONLY_AUTHOR_STATE_COUNT_SOURCE_STARTED',started_UTC=started))
    result=subprocess.run(args,input=source,text=True,capture_output=True)
    try:
        assert result.returncode==0,'AUTHOR_STATE_SOURCE_NOT_ACCEPTED'
        message=json.loads(result.stdout);assert message['kind']=='RAM_payload'
        payload=message['value'];result=None;message=None
    except Exception:
        save('private_failed_terminal.json',dict(status='AUTHOR_STATE_SOURCE_NOT_ACCEPTED_NO_EFFECTS',
             observed_UTC=datetime.now(timezone.utc).isoformat(),no_individual_stderr_or_payload_saved=True))
        raise RuntimeError('AUTHOR_STATE_SOURCE_NOT_ACCEPTED') from None
    save('private_source_coverage.json',payload['source_coverage']);save('private_measurement_proof.json',payload['proof'])
    save('private_terminal.json',dict(status='COMPLETE_AUTHOR_STATE_MEASUREMENT_CASE_COUNTS_RAM_WAITING_CONTRACT',
        started_UTC=started,completed_UTC=datetime.now(timezone.utc).isoformat(),
        freeze_SHA256=sha(freeze),coverage_SHA256=sha(out/'private_source_coverage.json'),proof_SHA256=sha(out/'private_measurement_proof.json'),
        individual_arrays_RAM_only=True,clinical_and_coordination_effects_not_computed=True))
    print(json.dumps(dict(status='AUTHOR_STATE_COUNTS_READY_WAITING_FROZEN_STATISTICS')),flush=True)
    path=task/'private_frozen_author_state_family_statistics_v1.json';deadline=time.monotonic()+21600
    while not path.exists():
        if time.monotonic()>deadline:
            save('private_governance_terminal.json',dict(status='CONTRACT_WAIT_EXPIRED_COUNTS_RELEASED_NO_EFFECTS'));return
        time.sleep(2)
    contract=json.loads(path.read_text())
    assert contract['measurement_terminal_SHA256']==sha(out/'private_terminal.json')
    if contract.get('do_not_fit'):
        save('private_governance_terminal.json',dict(status='EXPLICIT_MEASUREMENT_ONLY',reason=contract['do_not_fit']));return
    module=code/'NSCLC_author_state_family_statistics_v1.py'
    assert sha(module)==contract['statistics_module_SHA256']
    spec=importlib.util.spec_from_file_location('frozen_author_state_family',module);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    m.run(payload,task,contract)


if __name__=='__main__':main()
