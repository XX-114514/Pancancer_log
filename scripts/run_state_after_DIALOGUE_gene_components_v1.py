"""Freeze gene-wise quality projection for later LR-gene-masked receptor response."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time,subprocess
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();plan_path=Path('/tmp/state_after_DIALOGUE_gene_component_plan_v1.json');plan=json.loads(plan_path.read_text());assert sha(__file__)==plan['driver_sha256'];P=V/'state_after_DIALOGUE_frozen_projection_20261016_v1';start=time.monotonic()
while not (P/'inventory.json').exists():
 e=json.loads((V/'state_after_DIALOGUE_projection_chain_20261016_v1.execution.json').read_text());p=Path('/proc',str(e['child_pid']),'cmdline');assert p.exists() and e['producer'].encode() in p.read_bytes(),e;assert time.monotonic()-start<32*3600;time.sleep(20)
iv=json.loads((P/'inventory.json').read_text())
for f,h in iv['files'].items():assert sha(P/f)==h
script=Path('/tmp/export_state_after_DIALOGUE_gene_components_v1.R');assert sha(script)==plan['R_sha256'];O=V/'state_after_DIALOGUE_gene_components_20261016_v1';assert not O.exists()
r=subprocess.run([os.environ['PANCANCER_DIALOGUE_RSCRIPT'],str(script),str(V/'state_after_DIALOGUE_inputs_20261016_v3'),str(V/'state_after_DIALOGUE_models_20261016_v2'),str(P),str(O),str(R.parent/plan['runtime_directory']),str(R/'all_CM_state_native_input_preflight_20261014_v1/private_all_frozen_sample_node_units.tsv.gz')],timeout=4*3600);assert r.returncode==0;assert (O/'execution_complete.txt').exists();out=dict(status='FROZEN_GENE_CONTRIBUTIONS_FOR_LR_MASKING_NOT_LR_EFFECT_TESTED',observed_UTC=datetime.now(timezone.utc).isoformat(),plan_sha256=sha(plan_path),projection_inventory_sha256=sha(P/'inventory.json'),files={p.name:sha(p) for p in O.iterdir() if p.is_file()},no_LR_or_p_used_for_programme_selection=True);(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
