"""Bounded all18 comparison supervision of existing fits only."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,subprocess,time,sys
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);O=V/'all_lineage_cross_method_supervision_20261014_v1';assert not O.exists();O.mkdir();plan=Path('/tmp/all_lineage_method_comparison_plan_v1.json');p=json.loads(plan.read_text());majors=sorted(m for mm in p['blocks'].values() for m in mm);done=[];started=time.monotonic()
def write(name,obj):(O/name).write_text(json.dumps(obj,indent=2)+'\n')
while len(done)<18:
 for i,major in enumerate(majors):
  if major in done:continue
  block=next(b for b,mm in p['blocks'].items() if major in mm)
  ca=R/'all_lineage_cNMF_diagnostic_audit_20261014_v2/major_00/independent_audit.json' if i==0 else V/'all_lineage_cNMF_diagnostic_audit_20261014_v3'/f'major_{i:02}'/'independent_audit.json'
  sa=R/'all_lineage_Spectra_diagnostic_audit_20261014_v1'/block/'independent_audit.json' if block=='lymphoid_T' else V/'all_lineage_Spectra_diagnostic_audit_20261014_v2'/block/'independent_audit.json'
  if not ca.exists() or not sa.exists():continue
  with (O/f'major_{i:02}.log').open('x') as f:
   r=subprocess.run([sys.executable,'/tmp/compare_all_lineage_methods_v1.py','--major',major,'--plan',str(plan)],stdout=f,stderr=subprocess.STDOUT,timeout=1800)
  assert r.returncode==0,('comparison_failed',major,r.returncode)
  q=V/'all_lineage_cross_method_comparison_20261014_v1'/f'major_{i:02}'/'inventory.json';iv=json.loads(q.read_text());assert iv['status']=='COMPLETE_ALL_PERMITTED_METHOD_PAIRS_NOT_STATE_VALIDATION';done.append(major);write('progress.json',dict(status='RUNNING_ALL18',completed_majors=done,observed_UTC=datetime.now(timezone.utc).isoformat()))
 if len(done)==18:break
 assert time.monotonic()-started<8*3600,'bounded_8hour_deadline'
 e=json.loads((V/'remaining_state_diagnostic_supervision_20261014_v2.execution.json').read_text());pid=e.get('child_pid');proc=Path('/proc',str(pid),'cmdline');live=proc.exists() and b'supervise_remaining_state_diagnostics_v2.py' in proc.read_bytes()
 assert live or e.get('exit_code')==0,('parent_diagnostics_not_live_or_success',e)
 if e.get('exit_code')==0:raise RuntimeError('parent_complete_but_required_audits_missing')
 time.sleep(30)
write('inventory.json',dict(status='COMPLETE_ALL18_METHOD_COMPARISONS_NOT_STATES',completed_majors=done,observed_UTC=datetime.now(timezone.utc).isoformat(),formal_states=False))
