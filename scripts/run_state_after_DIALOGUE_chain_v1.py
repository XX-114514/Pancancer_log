"""Bounded prepare acceptance then official DIALOGUE; no clinical or LR selection."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,subprocess,time,hashlib
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=Path('/tmp/state_after_DIALOGUE_fit_plan_v1.json');plan=json.loads(P.read_text());assert sha(__file__)==plan['driver_sha256'];O=V/'state_after_DIALOGUE_chain_20261016_v1';assert not O.exists();O.mkdir();S=V/'state_after_DIALOGUE_inputs_20261016_v2';start=time.monotonic()
def mark(phase,**kw):(O/'progress.json').write_text(json.dumps(dict(phase=phase,observed_UTC=datetime.now(timezone.utc).isoformat(),**kw))+'\n')
mark('WAITING_ORIGINAL_PREPARE_V2')
while not (S/'inventory.json').exists():
 e=json.loads((V/'state_after_DIALOGUE_input_preparation_20261016_v2.execution.json').read_text());proc=Path('/proc',str(e['child_pid']),'cmdline');assert proc.exists() and b'prepare_state_after_DIALOGUE_v2.py' in proc.read_bytes(),e;assert time.monotonic()-start<8*3600;time.sleep(15)
assert sha('/tmp/audit_state_after_DIALOGUE_inputs_v1.py')==plan['input_auditor_sha256'];assert sha('/tmp/run_state_after_DIALOGUE_v1.R')==plan['runner_sha256'];runtime=R.parent/plan['runtime_directory'];assert sha(runtime/'sources/DIALOGUE-9c146cc.tar.gz')==plan['official_archive_sha256'];U=R/'all_CM_state_native_input_preflight_20261014_v1/private_all_frozen_sample_node_units.tsv.gz';fi=json.loads((U.parent/'inventory.json').read_text());assert sha(U)==fi['private_unit_table_sha256'];env=os.environ.copy();env['PANCANCER_DIALOGUE_RUNTIME_ROOT']=str(runtime)
mark('INDEPENDENT_ALL28_INPUT_READBACK')
with (O/'input_audit.log').open('x') as f:r=subprocess.run([os.sys.executable,'/tmp/audit_state_after_DIALOGUE_inputs_v1.py'],env=env,stdout=f,stderr=subprocess.STDOUT,timeout=4*3600)
assert r.returncode==0,r.returncode;a=V/'state_after_DIALOGUE_input_acceptance_20261016_v1/inventory.json';ai=json.loads(a.read_text());assert ai['status']=='ACCEPTED_ALL14_DIALOGUE_CONTEXTS_ALL28_COUNT_ID_PCA_QC_READBACK' and ai['input_inventory_sha256']==sha(S/'inventory.json')
mark('OFFICIAL_DIALOGUE_ALL14_DISCOVERY_CONTEXTS')
F=V/'state_after_DIALOGUE_models_20261016_v1'
with (O/'official_fit.log').open('x') as f:r=subprocess.run([os.environ['PANCANCER_DIALOGUE_RSCRIPT'],'/tmp/run_state_after_DIALOGUE_v1.R',str(S),str(F),str(runtime),str(U)],env=env,stdout=f,stderr=subprocess.STDOUT,timeout=24*3600)
assert r.returncode==0,r.returncode;assert (F/'execution_complete.txt').exists();files={str(p.relative_to(F)):sha(p) for p in F.rglob('*') if p.is_file()};out=dict(status='ALL14_OFFICIAL_DIALOGUE_CONTEXTS_ATTEMPTED_HOLDOUT_AND_BIOLOGY_PENDING',observed_UTC=datetime.now(timezone.utc).isoformat(),plan_sha256=sha(P),input_acceptance_sha256=sha(a),official_archive_sha256=plan['official_archive_sha256'],files=files,heldout_validated=False,biological_state_or_LR_validation=False);(F/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');mark('DISCOVERY_TERMINAL_REVIEW_PER_CONTEXT_STATUS_NOT_AUTOMATIC_MCP_ACCEPTANCE');print(json.dumps({k:v for k,v in out.items() if k!='files'}))
