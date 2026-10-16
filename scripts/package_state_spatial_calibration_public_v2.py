"""Publish reviewed aggregate evidence only, splitting TSVs without byte changes."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time,shutil
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);D=V/'state_spatial_partial_panel_calibration_20261016_v2';A=V/'state_spatial_partial_panel_independent_audit_20261016_v2'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for tick in range(480):
 if (A/'inventory.json').exists():break
 e=V/'state_spatial_partial_panel_independent_audit_20261016_v2.execution.json'
 if e.exists() and json.loads(e.read_text()).get('exit_code') is not None:raise RuntimeError('Independent calibration audit terminated without inventory')
 time.sleep(30)
else:raise TimeoutError('Four-hour independent-audit dependency wait')
ai=json.loads((A/'inventory.json').read_text());assert ai['status'].startswith('PASS_ALL91447') and ai['calibration_inventory_sha256']==sha(D/'inventory.json');di=json.loads((D/'inventory.json').read_text());O=V/'state_spatial_calibration_publication_20261016_v2';assert not O.exists();O.mkdir(mode=0o700);meta=[]
for origin,files in [(D,list(di['public_outputs'])),(A,['public_all_source_axis_independent_checks.tsv'])]:
 for name in files:
  source=origin/name;expected=di['public_outputs'][name] if origin==D else ai['public_output_sha256'];assert sha(source)==expected;data=source.read_bytes();lines=data.splitlines(keepends=True);header=lines[0];parts=[];buffer=[];size=len(header)
  for line in lines[1:]:
   assert len(header)+len(line)<3*2**20
   if buffer and size+len(line)>3*2**20:parts.append(header+b''.join(buffer));buffer=[];size=len(header)
   buffer.append(line);size+=len(line)
  if buffer or not parts:parts.append(header+b''.join(buffer))
  assert header+b''.join(part[len(header):] for part in parts)==data;folder=O/source.stem;folder.mkdir();out=[]
  for j,content in enumerate(parts):
   p=folder/f'part_{j:03d}.tsv';p.write_bytes(content);assert p.stat().st_size<=3*2**20;out.append(dict(relative_file=str(p.relative_to(O)),sha256=sha(p),bytes=p.stat().st_size,records=len(content.splitlines())-1))
  meta.append(dict(source_run=origin.name,source_file=name,original_sha256=expected,original_bytes=len(data),records=len(lines)-1,parts=out,byte_reconstruction='Use first part in full, then append each later part excluding its first header line. Exact original source SHA required.'))
for d,label in [(D,'calibration_inventory.json'),(A,'independent_audit_inventory.json')]:shutil.copyfile(d/'inventory.json',O/label)
out=dict(status='ALL_PUBLIC_AGGREGATE_TABLES_SHA_VERIFIED_SPLIT_UNDER3MIB_WITH_EXACT_BYTE_RECONSTRUCTION',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),calibration_inventory_sha256=sha(D/'inventory.json'),independent_audit_inventory_sha256=sha(A/'inventory.json'),sources=52,cells=91447,programmes=814,spatial_datasets=14,unique_symbol_axes=4,tables=meta,no_raw_counts_cell_rows_patient_rows_or_private_identity_tables=True,formal_states_or_spatial_biological_validation=False);(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='tables'}))
