"""Finite dependency stage: qualify a fixed native object after CRC/SHA acceptance.
This writes scientific receipts only and does not publish or schedule notifications.
"""
from pathlib import Path
import os,time,json,hashlib,traceback,sys
from datetime import datetime,timezone
from qualify_breast_xenium_v2 import run
w=Path(os.environ['BREAST_NATIVE_WORK']);root=w/'local_native_review_PRIVATE';native=w/'first_native_local_PRIVATE';mapping=root/'sources/mapping';deadline=time.monotonic()+6*3600
state={'UTC':datetime.now(timezone.utc).isoformat(),'status':'WAITING_FOR_COMPLETE_NATIVE_AND_MAPPING_MEMBER_RECEIPTS','finite_max_seconds':6*3600,'no_publication_or_notification':True}
(root/'finite_stage_started_PRIVATE.json').write_text(json.dumps(state,indent=2)+'\n')
try:
 while time.monotonic()<deadline:
  if (native/'complete_PRIVATE.json').exists() and (mapping/'complete_PRIVATE.json').exists():break
  time.sleep(20)
 else:raise TimeoutError('Fixed acquisition dependency deadline reached; partial files preserved')
 for folder in [native,mapping]:
  z=json.loads((folder/'complete_PRIVATE.json').read_text())
  for v in z['members']:
   p=folder/v['local_name'];assert p.stat().st_size==v['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==v['sha256'] and v['member_CRC_validated']
 run(root,'S01',native,'mapping')
 (root/'finite_stage_terminal_PRIVATE.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'status':'DESCRIPTIVE_NATIVE_REVIEW_COMPLETED_IDENTITY_INFERENCE_NOT_FINAL','Python':sys.version},indent=2)+'\n')
except Exception as e:
 (root/'finite_stage_failure_PRIVATE.txt').write_text(traceback.format_exc());(root/'finite_stage_terminal_PRIVATE.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'status':'FAILED_DEPENDENCY_OR_QUALIFICATION','exception_type':type(e).__name__},indent=2)+'\n');raise
