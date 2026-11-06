from pathlib import Path
import hashlib,json,datetime
r=Path(__file__).resolve().parents[1];files=[]
for d in ['results/global','results/cytokine','results/summary','results/figures','config','scripts','audit']:
 files.extend(p for p in (r/d).rglob('*') if p.is_file() and p.suffix not in ['.pyc'] and 'PRIVATE' not in p.name and p.name!='completed_outputs_SHA.tsv')
files.extend([r/'private/CMT_identity_stats.npz',r/'private/CMT_identity_stats_index.tsv',r/'private/full_count_block_SHA.tsv'])
out=r/'audit/completed_outputs_SHA.tsv'
with out.open('w') as f:
 f.write('path\tbytes\tsha256\n')
 for p in sorted(set(files)):
  h=hashlib.sha256()
  with p.open('rb') as b:
   for chunk in iter(lambda:b.read(4*1024*1024),b''):h.update(chunk)
  f.write(f'{p.relative_to(r).as_posix()}\t{p.stat().st_size}\t{h.hexdigest()}\n')
print(json.dumps({'UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':len(files),'manifest_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'B':'running_not_in_completed_output_manifest'}))
