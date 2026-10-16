"""Wait for accepted official fits, freeze-project every supported slot, then all-slot internal holdout."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time,subprocess
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);plan_path=Path('/tmp/state_after_DIALOGUE_projection_plan_v1.json');plan=json.loads(plan_path.read_text());assert sha(__file__)==plan['driver_sha256']
O=V/'state_after_DIALOGUE_projection_chain_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
mark=lambda phase:(O/'progress.json').write_text(json.dumps(dict(phase=phase,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
M=V/'state_after_DIALOGUE_models_20261016_v2';start=time.monotonic();mark('WAITING_ALL14_OFFICIAL_ATTEMPTS')
while not (M/'inventory.json').exists():
 e=json.loads((V/'state_after_DIALOGUE_chain_20261016_v2.execution.json').read_text());p=Path('/proc',str(e['child_pid']),'cmdline');assert p.exists() and e['producer'].encode() in p.read_bytes(),e;assert time.monotonic()-start<26*3600;time.sleep(20)
a=V/'state_after_DIALOGUE_input_acceptance_20261016_v2/inventory.json';assert json.loads(a.read_text())['status']=='ACCEPTED_ALL14_DIALOGUE_CONTEXTS_ALL28_COUNT_ID_PCA_QC_READBACK';iv=json.loads((M/'inventory.json').read_text());assert iv['input_acceptance_sha256']==sha(a)
for rel,h in iv['files'].items():assert sha(M/rel)==h,rel
for f,h in plan['scripts'].items():assert sha(Path('/tmp')/f)==h,f
S=V/'state_after_DIALOGUE_inputs_20261016_v3';P=V/'state_after_DIALOGUE_frozen_projection_20261016_v1';T=V/'state_after_DIALOGUE_internal_holdout_20261016_v1';U=R/'all_CM_state_native_input_preflight_20261014_v1/private_all_frozen_sample_node_units.tsv.gz'
assert sha(U)==json.loads((U.parent/'inventory.json').read_text())['private_unit_table_sha256']
mark('FROZEN_SCORE_PARITY_AND_PROJECTION')
with (O/'projection.log').open('x') as f:r=subprocess.run([os.environ['PANCANCER_DIALOGUE_RSCRIPT'],'/tmp/project_state_after_DIALOGUE_frozen_v1.R',str(S),str(M),str(P),str(R.parent/plan['runtime_directory']),str(U)],stdout=f,stderr=subprocess.STDOUT,timeout=4*3600)
assert r.returncode==0,r.returncode;assert (P/'execution_complete.txt').exists()
# Slot parity failures remain explicit missing; no algorithm or covariate fallback.
mark('ALL28_FIXED_INTERNAL_HOLDOUT_SLOTS')
with (O/'holdout.log').open('x') as f:r=subprocess.run([os.sys.executable,'/tmp/test_state_after_DIALOGUE_holdout_v1.py',str(P),str(M),str(T)],stdout=f,stderr=subprocess.STDOUT,timeout=2*3600)
assert r.returncode==0,r.returncode;assert json.loads((T/'inventory.json').read_text())['fixed_slots']==28
files={str(p.relative_to(P)):sha(p) for p in P.rglob('*') if p.is_file()}
(P/'inventory.json').write_text(json.dumps(dict(status='ALL_SUPPORTED_SLOTS_ATTEMPTED_REVIEW_PER_SLOT_PARITY_AND_INTERNAL_HOLDOUT_NOT_EXTERNAL_VALIDATION',observed_UTC=datetime.now(timezone.utc).isoformat(),plan_sha256=sha(plan_path),model_inventory_sha256=sha(M/'inventory.json'),input_acceptance_sha256=sha(a),files=files,holdout_inventory_sha256=sha(T/'inventory.json')),indent=2)+'\n');mark('TERMINAL_PER_SLOT_SCIENTIFIC_REVIEW_REQUIRED')
