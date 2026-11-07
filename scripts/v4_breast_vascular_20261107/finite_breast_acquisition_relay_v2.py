"""Finite public-data acquisition and descriptive schema qualification controller."""
from pathlib import Path
import sys,json,traceback,shutil
from datetime import datetime,timezone
from selective_public_zip import download_selected
from qualify_breast_xenium_v2 import run
r=Path(sys.argv[1]);config=json.loads((r/'config/archive_members_relay_v2_PRIVATE.json').read_text());assert shutil.disk_usage(r).free>500e9
log=r/'launch/events_relay_v2_PRIVATE.jsonl'
def event(**z):
 with log.open('a') as f:f.write(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),**z})+'\n')
try:
 m=config['mapping'];download_selected(m['url'],m['members'],r/'sources/mapping_v2',350*1024**2);event(stage='mapping',status='CRC_SHA_complete')
 for j in config['xenium']:
  u=j['unit_alias'];event(stage=u,status='start_download');out=r/'sources'/u
  download_selected(j['url'],[j['member']],out,int(j['compressed_bytes']+40*1024**2));event(stage=u,status='CRC_SHA_complete')
  try:run(r,u,out,'mapping_v2');event(stage=u,status='descriptive_complete')
  except Exception as e:event(stage=u,status='descriptive_failed',exception_type=type(e).__name__);(r/'launch'/(u+'_failure_PRIVATE.txt')).write_text(traceback.format_exc())
 (r/'launch/terminal_relay_v2.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'status':'FINITE_ALL_SIX_DOWNLOAD_STAGES_ENDED_REVIEW_UNIT_RECEIPTS'},indent=2))
except Exception as e:
 event(stage='controller',status='failed',exception_type=type(e).__name__);(r/'launch/controller_relay_v2_failure_PRIVATE.txt').write_text(traceback.format_exc());raise
