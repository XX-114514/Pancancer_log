"""Reuse all18 accepted mappings and first ongoing native retry; robust remaining5."""
from pathlib import Path
import os,json,time,hashlib,traceback
from datetime import datetime,timezone
from selective_public_zip_v3 import download_selected
b=Path(os.environ['BREAST_PARENT_WORK']);r=b/'local_native_review_PRIVATE';deadline=time.monotonic()+18*3600;eventfile=r/'cohort_events_v2_PRIVATE.jsonl'
def event(**z):
 with eventfile.open('a') as f:f.write(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),**z})+'\n')
x=sorted([v for v in json.loads((b/'sources_PRIVATE/anndata_objects.zip_inventory_PRIVATE.json').read_text())['members'] if '/xenium/' in v['member']],key=lambda v:v['member']);m=json.loads((r/'sources/mapping_all18/complete_PRIVATE.json').read_text());assert len(m['members'])==18
for v in m['members']:assert hashlib.sha256((r/'sources/mapping_all18'/v['local_name']).read_bytes()).hexdigest()==v['sha256']
event(stage='remaining5',status='ROBUST_BODY_RETRY_QUEUE_STARTED_EXISTING18_MAPPING_AND_FIRST_NATIVE_REUSED')
while not (b/'first_native_local_PRIVATE/complete_PRIVATE.json').exists():
 if time.monotonic()>deadline:raise TimeoutError('First native complete-source dependency deadline')
 time.sleep(30)
for i,z in enumerate(x[1:],2):
 u=f'S{i:02d}';out=r/'sources'/u
 if (out/'complete_PRIVATE.json').exists():event(stage=u,status='REUSE_COMPLETE_CRC_SHA_SOURCE');continue
 if time.monotonic()>deadline:raise TimeoutError('Finite acquisition deadline')
 try:
  event(stage=u,status='START_NATIVE_DOWNLOAD');download_selected('https://zenodo.org/records/18917100/files/anndata_objects.zip?download=1',[z['member']],out,int(z['compressed_bytes']+40*1024**2));event(stage=u,status='NATIVE_CRC_SHA_ACCEPTED')
 except Exception as e:(r/(u+'_download_v2_failure_PRIVATE.txt')).write_text(traceback.format_exc());event(stage=u,status='DOWNLOAD_FAILED_PARTIAL_RETAINED',exception_type=type(e).__name__)
(r/'cohort_terminal_v2_PRIVATE.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'status':'FINITE_SOURCE_STAGES_ENDED_CHECK_INDIVIDUAL_RECEIPTS','requested_sources':6,'all18_mapping_reused':True},indent=2)+'\n')
