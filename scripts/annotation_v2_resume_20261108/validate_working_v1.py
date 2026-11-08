import pathlib,sys,json,hashlib,collections,datetime
import pandas as pd
h,r=map(pathlib.Path,sys.argv[1:3]);out=r/'working_annotation_v1';assert json.loads((r/'consensus_replay_acceptance.json').read_text())['status']=='passed'
orig=pd.read_csv(h/'p0/cell_labels_r3.tsv.gz',sep='\t',usecols=['cell_id'],dtype=str).cell_id
expected=set(orig);assert len(expected)==len(orig)
main=set()
for pool in ['T_NK','B_Plasma','Stromal','Myeloid','Endothelial']:
 c=pd.read_csv(r/f'consensus_v1/{pool}.cells.tsv.gz',sep='\t',usecols=['cell_id','include_main_subtype']);main.update(c.loc[c.include_main_subtype,'cell_id'])
seen=set();counts=collections.Counter();nmain=0
for d in pd.read_csv(out/'cell_labels_working_v1.tsv.gz',sep='\t',dtype=str,keep_default_na=False,chunksize=200000):
 ids=set(d.cell_id);assert len(ids)==len(d) and not (seen&ids);seen.update(ids)
 actual=set(d.loc[d.include_main_subtype=='True','cell_id']);assert actual==(main&ids)
 assert not ((d.include_main_subtype=='True')&(d.annotation_status!='observed_stable')).any()
 assert not ((d.annotation_status=='transferred')&(pd.to_numeric(d.transfer_conf,errors='coerce')<.6)).any()
 assert (d.loc[d.annotation_status.isin(['set_aside_review','excluded_qc']),'include_lineage_total']=='False').all()
 nmain+=len(actual);counts.update(d.annotation_status)
assert seen==expected and nmain==len(main)
summary=json.loads((out/'summary.json').read_text());assert summary['status_counts']==dict(counts)
def sha(p):
 a=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):a.update(b)
 return a.hexdigest()
files=[]
for p in sorted(out.iterdir()):
 if p.is_file():files.append({'path':p.name,'bytes':p.stat().st_size,'sha256':sha(p)})
a={'status':'accepted_working_sidecar','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'rows':len(seen),'main_subtype_candidates':nmain,'cell_id_set_preserved':True,'consensus_main_membership_exact':True,'status_counts':dict(counts),'scientific_final':False,'files':files}
(out/'acceptance.json').write_text(json.dumps(a,indent=2));print(json.dumps(a),flush=True)
