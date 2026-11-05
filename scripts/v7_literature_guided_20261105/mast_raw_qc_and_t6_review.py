import os,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))))
r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1'
o=r/'stage_d_raw_marker_review';o.mkdir(exist_ok=False)
a=ad.read_h5ad(r/'stage_b_mast_fullgene/mast_v7_fullgene_raw.h5ad');X=a.X.tocsr();genes=a.var_names
programs={'Mast':['TPSAB1','TPSB2','CPA3','KIT','MS4A2'],'T':['CD3D','CD3E','TRAC'],'Epithelial':['EPCAM','KRT8','KRT18','KRT19'],'Myeloid_shared':['LYZ','LST1','TYROBP']}
vals={};rows=[]
for threshold in [1,2,3]:
 for name,gs in programs.items():
  vals[(name,threshold)]=np.asarray((X[:,genes.get_indexer(gs)]>=threshold).sum(axis=1)).ravel()>=2
 for label,s in a.obs.groupby('final_annotation_v7',observed=True):
  ix=a.obs.index.get_indexer(s.index);z={'label_v7':label,'n_cells':len(ix),'min_raw_count_per_gene':threshold}
  for name in programs:z[name+'_ge2_genes_n']=int(vals[(name,threshold)][ix].sum())
  for name in ['T','Epithelial','Myeloid_shared']:z['Mast_AND_'+name+'_n']=int((vals[('Mast',threshold)]&vals[(name,threshold)])[ix].sum())
  rows.append(z)
q=pd.DataFrame(rows);q.to_csv(o/'mast_raw_same_cell_marker_sensitivity.tsv',sep='\t',index=False)
# Old candidate-panel review: all passing anchors are retained, regardless of winner.
u=ROOT/'runs/20260818_pancancer_5M_merge_prep/07_fine_annotation/20260905_uncertainty_resolution_v1'
d=pd.read_csv(u/'tables/expanded_panel_gene_evidence_v1.tsv.gz',sep='\t');d=d[d.cluster_uid.isin(['Tcell_4','Tcell_6'])]
d.to_csv(o/'T4_T6_legacy_all_panel_gene_evidence.tsv',sep='\t',index=False)
res=[]
for (cl,sub),s in d.groupby(['cluster_uid','subtype']):
 p=s[s.passes_threshold];res.append({'cluster_uid':cl,'subtype':sub,'n_panel_genes':len(s),'n_pass':len(p),'n_canonical_pass':int(p.is_canonical.sum()),'passed_genes':','.join(p.gene),'old_score_mean_passing_only':float(p.gene_evidence_score.mean()) if len(p) else None,'warning':'Diagnostic candidates only; sample support for non-winner not previously evaluated'})
pd.DataFrame(res).to_csv(o/'T4_T6_all_panel_candidates.tsv',sep='\t',index=False)
manifest={'stage':'D_completed_raw_marker_and_candidate_diagnostic','UTC':datetime.now(timezone.utc).isoformat(),'n_mast_cells':len(a),'n_annotation_changes':0,'negative_marker_policy':'Core raw source measurement mask not yet joined for extra genes; this stage uses positive detections only.','outputs':[{'name':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in o.iterdir()]}
(o/'complete.json').write_text(json.dumps(manifest,indent=2));print(json.dumps({'manifest':manifest,'raw_Mast':q.to_json(orient='records'),'T_candidates':res}))
