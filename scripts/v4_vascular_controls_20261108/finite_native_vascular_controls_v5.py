"""Bounded dependency stage for six original donors; no Git/alerts or state models."""
from pathlib import Path
import os,json,time,subprocess,traceback,hashlib,csv
from datetime import datetime,timezone
from run_native_vascular_controls_v5 import run
from vascular_spatial_controls_v2 import patient_summary
w=Path(os.environ['VASCULAR_CONTROLS_WORK']);b=Path(os.environ['BREAST_PARENT_WORK']);root=b/'local_native_review_PRIVATE';deadline=time.monotonic()+18*3600;out=w/'native_controls_v5_PRIVATE';out.mkdir(exist_ok=False);done={};log=out/'events_PRIVATE.jsonl'
xi=sorted([v for v in json.loads((b/'sources_PRIVATE/anndata_objects.zip_inventory_PRIVATE.json').read_text())['members'] if '/xenium/' in v['member']],key=lambda v:v['member']);sheet=list(csv.DictReader((b/'sources_PRIVATE/author_code/config/sample_sheet_xenium.csv').open()));bridge={v['section_id']:v['donor_id'] for v in sheet};patients={f'S{i:02d}':bridge[Path(v['member']).stem] for i,v in enumerate(xi,1)};probe_genes={v['gene_name'] for v in csv.DictReader((b/'sources_PRIVATE/author_panel.txt').open())};assert len(set(patients.values()))==6
(out/'patient_bridge_PRIVATE.json').write_text(json.dumps(patients,indent=2))
def event(**z):
 with log.open('a') as f:f.write(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),**z})+'\n')
event(stage='all6',status='WAITING_COMPLETE_CRC_SHA_NO_STATE_PREREQUISITE',plan_sha256=hashlib.sha256((w/'vascular_controls_plan_before_native_effects_v5.json').read_bytes()).hexdigest())
while time.monotonic()<deadline and len(done)<6:
 for i in range(1,7):
  u=f'S{i:02d}'
  if u in done:continue
  native=b/'first_native_local_PRIVATE' if i==1 else root/'sources'/u;mapping=root/'sources/mapping_all18'
  if not (native/'complete_PRIVATE.json').exists() or not (mapping/'complete_PRIVATE.json').exists():continue
  try:
   event(stage=u,status='START_NATIVE_IDENTITY_AND_VASCULAR_CONTROLS');z=run(out,u,native,mapping,expected_member=xi[i-1]["member"],panel_gene_names=probe_genes);done[u]='complete';event(stage=u,status='CASE_COMPUTATION_COMPLETE_BIOLOGY_AND_VISUAL_REVIEW_PENDING')
   subprocess.run([os.environ['VASCULAR_PLOTTING_PYTHON'],str(w/'render_vascular_identity_review_v3.py'),str(out/'spatial_controls_v1'/u)],check=True,env={**os.environ,'MPLCONFIGDIR':os.environ['VASCULAR_MPL_CACHE']})
  except Exception as e:done[u]='failed';(out/(u+'_failure_PRIVATE.txt')).write_text(traceback.format_exc());event(stage=u,status='failed',exception_type=type(e).__name__)
 if len(done)<6:time.sleep(30)
accepted=[u for u,status in done.items() if status=='complete']
if accepted:
 patient_summary([out/'spatial_controls_v1'/u for u in accepted],[patients[u] for u in accepted],out/'cohort_private_summary')
(out/'terminal_PRIVATE.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'requested_donors':6,'native_cases_computed':len(accepted),'failed_or_pending':6-len(accepted),'minimum5_patients_required_for_formal_P':True,'status':'FINITE_STAGE_ENDED_REVIEW_ALL66_HYPOTHESES_AND_CASE_QUALIFICATION'},indent=2)+'\n')
