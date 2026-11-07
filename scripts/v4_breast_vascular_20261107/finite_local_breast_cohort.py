"""Finite local fallback for all six frozen native inputs; no automatic Git or alerts."""
from pathlib import Path
import os,json,time,hashlib,traceback
from datetime import datetime,timezone
from selective_public_zip import download_selected
from qualify_breast_xenium_v2 import run
w=Path(os.environ['BREAST_NATIVE_WORK']);r=w/'local_native_review_PRIVATE';s=w/'sources_PRIVATE';deadline=time.monotonic()+18*3600
members=json.loads((s/'cell_type_mapping.zip_inventory_PRIVATE.json').read_text())['members'];allmaps=sorted(z['member'] for z in members if '/xenium/' in z['member']);xi=sorted([z for z in json.loads((s/'anndata_objects.zip_inventory_PRIVATE.json').read_text())['members'] if '/xenium/' in z['member']],key=lambda z:z['member']);assert len(xi)==6 and len(allmaps)==18
log=r/'cohort_events_PRIVATE.jsonl'
def event(**z):
 with log.open('a') as f:f.write(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),**z})+'\n')
def wait_for(p):
 while not p.exists():
  if time.monotonic()>deadline:raise TimeoutError('Finite source dependency deadline')
  time.sleep(20)
try:
 event(stage='cohort',status='STARTED_LOCAL_FALLBACK_SOURCE_AND_SCIENCE_UNCHANGED')
 wait_for(r/'sources/mapping/complete_PRIVATE.json');old=json.loads((r/'sources/mapping/complete_PRIVATE.json').read_text());remaining=[n for n in allmaps if n not in {m['member'] for m in old['members']}];assert len(remaining)==15
 more=download_selected('https://zenodo.org/records/18917100/files/cell_type_mapping.zip?download=1',remaining,r/'sources/mapping_remaining',350*1024**2);event(stage='mapping_remaining',status='15_CRC_SHA_accepted')
 full=r/'sources/mapping_all18';full.mkdir(exist_ok=False);rows=[]
 for folder,receipt in [(r/'sources/mapping',old),(r/'sources/mapping_remaining',more)]:
  for v in receipt['members']:
   src=folder/v['local_name'];assert src.stat().st_size==v['bytes'] and hashlib.sha256(src.read_bytes()).hexdigest()==v['sha256'];name=f'member_{len(rows):03d}_PRIVATE.csv';os.link(src,full/name);rows.append({**v,'local_name':name})
 (full/'complete_PRIVATE.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'status':'18_MEMBERS_CRC_SHA_ACCEPTED_COMBINED_WITHOUT_DATA_CHANGE','members':rows},indent=2)+'\n')
 wait_for(w/'first_native_local_PRIVATE/complete_PRIVATE.json');event(stage='S01',status='native_CRC_SHA_accepted_existing_finite_qualification_handles_first')
 for i,z in enumerate(xi[1:],2):
  if time.monotonic()>deadline:raise TimeoutError('Finite cohort deadline')
  u=f'S{i:02d}';out=r/'sources'/u;event(stage=u,status='start_native_download')
  download_selected('https://zenodo.org/records/18917100/files/anndata_objects.zip?download=1',[z['member']],out,int(z['compressed_bytes']+40*1024**2));event(stage=u,status='native_CRC_SHA_accepted')
  try:run(r,u,out,'mapping_all18');event(stage=u,status='descriptive_qualification_complete_not_final_identity')
  except Exception as e:(r/(u+'_qualification_failure_PRIVATE.txt')).write_text(traceback.format_exc());event(stage=u,status='qualification_failed',exception_type=type(e).__name__)
 (r/'cohort_terminal_PRIVATE.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'status':'FINITE_ALL_SIX_NATIVE_ACQUISITION_ENDED_REVIEW_INDIVIDUAL_QUALIFICATIONS'},indent=2)+'\n')
except Exception as e:
 (r/'cohort_failure_PRIVATE.txt').write_text(traceback.format_exc());(r/'cohort_terminal_PRIVATE.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'status':'FAILED_IN_ACQUISITION_OR_QUALIFICATION','exception_type':type(e).__name__},indent=2)+'\n');raise
