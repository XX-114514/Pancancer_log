from pathlib import Path
import hashlib,json,pandas as pd
r=Path(__file__).resolve().parents[1];pc=r.parents[2];f=pc/'Project_v5/frozen/covarnet';orig=pc/'Project_v2/CoVarNet/outputs/Covarnet_5M_v5_malignancy_v1_20260825/tables';rows=[]
for name in ['nmf_W_subcluster_by_cm.tsv','nmf_H_cm_by_sample.tsv','covarnet_each_nodes.tsv','covarnet_each_edges.tsv']:
 a=f/name;b=orig/name;rows.append({'file':name,'frozen_sha256':hashlib.sha256(a.read_bytes()).hexdigest(),'original_sha256':hashlib.sha256(b.read_bytes()).hexdigest() if b.exists() else None})
for a in rows:a['match']=a['frozen_sha256']==a['original_sha256']
(r/'audit/original_vs_frozen.json').write_text(json.dumps(rows,indent=2));print(rows)
# Full original-LIANA file can support readouts only within old top10; not a full-LR replacement.
v=r.parent/'20260818_pancancer_5M_merge_prep/10_final_v5_malignancy_covarnet_liana/20260825_v5_final_v1'
p=v/'liana_cm_v3/tables/all_within_cm_interactions.tsv.gz';print('LIANA_FULL',p.exists(),p.stat().st_size if p.exists() else 0)
if p.exists():
 d=pd.read_csv(p,sep='\t',nrows=0);(r/'audit/liana_reuse_schema.json').write_text(json.dumps({'exists':True,'bytes':p.stat().st_size,'columns':list(d),'full_database_family':False,'reason':'old_top10_and_return_all_lrs_false; neither full CPDB universe nor patient inference'},indent=2))
