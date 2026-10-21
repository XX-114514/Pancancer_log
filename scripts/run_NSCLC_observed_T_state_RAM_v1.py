"""Hold qualified case data in RAM while a statistics contract is frozen.

Only aggregate coverage/proofs/checkpoints are saved. A checked named module
can consume the RAM payload after the new clinical contract exists. No lock or
checkpoint file by itself is a liveness claim; monitor the running tool handle.
"""
from datetime import datetime,timezone
from pathlib import Path
import ast
import hashlib
import importlib.util
import json
import os
import shlex
import subprocess
import sys
import time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def remote_program(code):
    original=(code/'prepare_NSCLC_full_native_gene_units_v1.py').read_text();tree=ast.parse(original)
    selected=[n for n in tree.body if isinstance(n,(ast.ClassDef,ast.FunctionDef)) and n.name in ['TrackedReader','chunks']]
    assert len(selected)==2
    s='import types\nevidence=types.ModuleType("fixed_evidence")\n'
    s+='exec('+repr((code/'observed_T_marker_evidence_v1.py').read_text())+',evidence.__dict__)\n'
    s+=(code/'NSCLC_observed_T_state_RAM_reader_v1.py').read_text()+'\n'
    s+='\n'.join(ast.get_source_segment(original,n) for n in selected)+'\n'
    s+='def progress(v):\n print(json.dumps({"kind":"progress","value":v}),flush=True)\n'
    s+='data=read_source('+repr(os.environ['PANCANCER_PREFLIGHT_ROOT'])+',progress)\n'
    s+='print(json.dumps({"kind":"RAM_payload","value":data},allow_nan=False),flush=True)\n'
    return s


def main():
    task=Path(os.environ['PANCANCER_NSCLC_STATE_TASK']);code=Path(__file__).parent
    freeze=task/'private_frozen_full_native_measurement_v1.json';contract=json.loads(freeze.read_text())
    for n,d in contract['code_SHA256'].items():assert sha(code/n)==d
    out=task/'native_measurement_v1';out.mkdir(exist_ok=False)
    save=lambda n,v:(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    started=datetime.now(timezone.utc).isoformat();payload=None
    save('private_checkpoint.json',dict(stage='START_READONLY_SOURCE',started_UTC=started))
    args=['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','ConnectTimeout=15']
    if os.environ.get('PANCANCER_READONLY_PROXY_COMMAND'):args+=['-o','ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND']]
    args+=[os.environ['PANCANCER_READONLY_SSH_TARGET'],'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 '+shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -']
    proc=subprocess.Popen(args,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    source=remote_program(code);proc.stdin.write(source);proc.stdin.close()
    for line in proc.stdout:
        obj=json.loads(line)
        if obj['kind']=='progress':
            save('private_checkpoint.json',dict(**obj['value'],observed_UTC=datetime.now(timezone.utc).isoformat(),started_UTC=started) if 'started_UTC' not in obj['value'] else dict(**obj['value'],observed_UTC=datetime.now(timezone.utc).isoformat()))
        else:
            assert obj['kind']=='RAM_payload' and payload is None;payload=obj['value']
    rc=proc.wait();errors=proc.stderr.read()
    if rc or payload is None:
        # Full stderr can contain individual keys; never persist or print it.
        save('private_failed_terminal.json',dict(status='NATIVE_SOURCE_NOT_ACCEPTED',exit_code=rc,
            last_stage=json.loads((out/'private_checkpoint.json').read_text())['stage'],observed_UTC=datetime.now(timezone.utc).isoformat(),
            original_source_and_control_preserved=True,no_source_individual_output_saved=True))
        raise RuntimeError('READONLY_NATIVE_SOURCE_NOT_ACCEPTED')
    save('private_complete_source_coverage.json',payload['coverage'])
    save('private_source_input_proof.json',payload['proof'])
    save('private_measurement_terminal.json',dict(status='COMPLETE_NATIVE_MARKER_CAPTURE_SOURCE_PROOF_ACCEPTED_CASE_ARRAYS_RAM',
        started_UTC=started,completed_UTC=datetime.now(timezone.utc).isoformat(),frozen_SHA256=sha(freeze),
        proof_SHA256=sha(out/'private_source_input_proof.json'),coverage_SHA256=sha(out/'private_complete_source_coverage.json'),
        individual_case_counts_quality_and_clinical_data_RAM_only=True,new_clinical_effects_computed=False))
    save('private_checkpoint.json',dict(stage='WAITING_FROZEN_STATISTICS_CONTRACT_CASE_ARRAYS_RAM_ONLY',observed_UTC=datetime.now(timezone.utc).isoformat()))
    print(json.dumps(dict(status='NATIVE_CAPTURE_COMPLETE_WAITING_NEW_STATISTICS_CONTRACT')),flush=True)
    path=task/'private_frozen_NSCLC_state_statistics_v1.json';deadline=time.monotonic()+21600
    while not path.exists():
        if time.monotonic()>deadline:
            save('private_governance_terminal.json',dict(status='CONTRACT_WAIT_EXPIRED_NO_EFFECTS_CASE_ARRAYS_RELEASED',observed_UTC=datetime.now(timezone.utc).isoformat()))
            return
        time.sleep(2)
    plan=json.loads(path.read_text())
    assert plan['measurement_terminal_SHA256']==sha(out/'private_measurement_terminal.json')
    if plan.get('do_not_fit'):
        save('private_governance_terminal.json',dict(status='MEASUREMENT_ONLY_EXPLICIT_NO_FIT',reason=plan['do_not_fit'],observed_UTC=datetime.now(timezone.utc).isoformat()));return
    module=code/'NSCLC_observed_T_state_statistics_v1.py'
    assert sha(module)==plan['statistics_module_SHA256']
    spec=importlib.util.spec_from_file_location('frozen_NSCLC_state_statistics',module);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    m.run(payload,task,plan)


if __name__=='__main__':main()
