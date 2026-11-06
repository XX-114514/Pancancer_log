"""Canonical SHA256 manifest of the immutable native Zarr source; read-only source scan."""
from pathlib import Path
import hashlib,json,time,datetime
r=Path(__file__).resolve().parents[1];src=r.parent/'20260818_pancancer_5M_merge_prep/02_merge_new_cohort/objects/5M_merge_concat.zarr';out=r/'private/native_zarr_file_SHA.tsv';start=time.time();total=0;count=0
with out.open('w') as f:
 f.write('sha256\tbytes\tpath\n')
 for p in sorted(src.rglob('*')):
  if not p.is_file():continue
  before=p.stat();h=hashlib.sha256()
  with p.open('rb') as handle:
   for chunk in iter(lambda:handle.read(8*1024*1024),b''):h.update(chunk)
  after=p.stat();assert before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns
  f.write(f'{h.hexdigest()}\t{before.st_size}\t{p.relative_to(src).as_posix()}\n');total+=before.st_size;count+=1
  if count%1000==0:print(count,total,round(time.time()-start),flush=True)
receipt={'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':'${PC_ROOT}/Project_v3/runs/20260818_pancancer_5M_merge_prep/02_merge_new_cohort/objects/5M_merge_concat.zarr','manifest_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'format':'sorted relative paths; SHA256<TAB>bytes<TAB>path; header; LF; SHA256 of complete manifest','files':count,'bytes':total,'seconds':time.time()-start,'source_unchanged_during_each_file_read':True};(r/'audit/native_expression_tree_SHA.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
