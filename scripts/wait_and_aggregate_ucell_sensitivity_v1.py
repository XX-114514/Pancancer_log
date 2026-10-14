"""Wait on the specific original projection chain, then audit once; never restart fits."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,subprocess,time,sys
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);O=V/'all_source_ucell_sensitivity_supervision_20261014_v1';assert not O.exists();O.mkdir();start=time.monotonic()
def mark(phase):(O/'progress.json').write_text(json.dumps(dict(phase=phase,observed_UTC=datetime.now(timezone.utc).isoformat()),indent=2)+'\n')
mark('WAITING_ORIGINAL_FULL_PROJECTION_CHAIN')
while not (V/'all_source_native_UCell_audit_aggregation_20261014_v2/inventory.json').exists():
 e=json.loads((V/'all_source_ucell_audit_supervision_20261014_v2.execution.json').read_text());p=Path('/proc',str(e.get('child_pid')),'cmdline');assert p.exists() and b'wait_and_audit_all_source_ucell_v2.py' in p.read_bytes(),('projection_chain_not_live',e);assert time.monotonic()-start<36*3600,'bounded36hour_wait';time.sleep(30)
mark('INDEPENDENT_NATIVE_SCORE_AUDIT_AND_PRIVATE_UNIT_AGGREGATION')
with (O/'audit.log').open('x') as f:r=subprocess.run([sys.executable,'/tmp/aggregate_ucell_sensitivity_units_v1.py'],stdout=f,stderr=subprocess.STDOUT,timeout=24*3600)
assert r.returncode==0,r.returncode
mark('AUDIT_COMPLETE_PENDING_BIOLOGICAL_STATE_AND_ALL9CM_REVIEW')
