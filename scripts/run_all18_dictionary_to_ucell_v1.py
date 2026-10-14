"""Continue the already authorized full task after original model audits finish."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,subprocess,time
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);O=V/'all18_dictionary_to_ucell_execution_20261014_v1';assert not O.exists();O.mkdir();start=time.monotonic()
def mark(phase,**extra):(O/'progress.json').write_text(json.dumps(dict(phase=phase,observed_UTC=datetime.now(timezone.utc).isoformat(),**extra),indent=2)+'\n')
mark('WAITING_ORIGINAL24_DIAGNOSTIC_SCOPES')
while True:
 e=json.loads((V/'remaining_state_diagnostic_supervision_20261014_v2.execution.json').read_text());p=V/'remaining_state_diagnostic_supervision_20261014_v2/inventory.json'
 if p.exists():
  q=json.loads(p.read_text());assert q['completed_scopes']==24 and q['status']=='COMPLETE_ALL6_BLOCKS_AND18_MAJORS_FIXED_NUMERIC_DIAGNOSTICS_NOT_STATES';break
 cmd=Path('/proc',str(e.get('child_pid')),'cmdline');assert cmd.exists() and b'supervise_remaining_state_diagnostics_v2.py' in cmd.read_bytes(),('original_diagnostics_not_live',e);assert time.monotonic()-start<8*3600,'eight_hour_wait_limit';time.sleep(30)
scanpy=os.environ['PANCANCER_SCANPY_PYTHON'];ucell=str(R/'pyucell_runtime_20261014_v1/venv/bin/python');env=os.environ.copy();env.update(PANCANCER_ALL_STATE_DICTIONARY_PLAN='/tmp/all_lineage_state_candidate_dictionary_plan_v4.json',PANCANCER_ALL_LINEAGE_SPECTRA_DIAGNOSTIC_PLAN='/tmp/all_lineage_spectra_diagnostic_plan_v2.json',PANCANCER_UCELL_FULL_PROJECTION_PLAN='/tmp/all_source_ucell_full_projection_plan_v1.json')
for phase,python,script,timeout in [('ASSEMBLE_ALL18',scanpy,'assemble_all_lineage_state_candidate_dictionary_v4.py',3600),('AUDIT_ALL_SIGNATURES',scanpy,'audit_all18_candidate_dictionary_v1.py',3600),('PROJECT_ALL4532933',ucell,'project_all_source_ucell_stream_v1.py',24*3600)]:
 mark(phase)
 with (O/(phase+'.log')).open('x') as f:r=subprocess.run([python,'/tmp/'+script],env=env,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
 assert r.returncode==0,(phase,r.returncode)
mark('FULL_PROJECTION_COMPLETE_PENDING_INDEPENDENT_SCORE_AUDIT_NOT_STATES')
