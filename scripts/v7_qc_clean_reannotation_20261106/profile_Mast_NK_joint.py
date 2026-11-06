import os
from pathlib import Path
PROJECT_ROOT=Path(os.environ["PROJECT_ROOT"])
import json
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad
run=PROJECT_ROOT/'runs/20261106_v7_QC_clean_reannotation_v1';a=ad.read_h5ad(run/'stage_c_clean_views/Mast/Mast_clean_raw.h5ad');q=pd.read_pickle(run/'stage_f_Mast_expanded_context_v1/additional_context_cellwise_PRIVATE.pkl.gz');assert a.obs_names.equals(q.index)
genes=['KLRD1','NCR1','GNLY','PRF1','NKG7','CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'];ix=a.var_names.get_indexer(genes);assert (ix>=0).all();x=a.X[:,ix].toarray()
q['NK_receptor_component_raw2']=(x[:,:2]>=2).any(axis=1);q['NK_effector_component_raw2']=(x[:,2:4]>=2).any(axis=1);q['NK_NKG7_raw2']=x[:,4]>=2;q['NK_three_component_positive_raw2']=q.NK_receptor_component_raw2&q.NK_effector_component_raw2&q.NK_NKG7_raw2
pct=pd.to_numeric(q.doublet_within_sample_method_percentile,errors='coerce');q['Mast_NK_joint_review']=q.NK_three_component_positive_raw2&(pct>=.95)
q.to_pickle(run/'stage_f_Mast_expanded_context_v1/NK_joint_evidence_PRIVATE.pkl.gz');rows=q.groupby('clean_cluster',observed=True).agg(n_cells=('NK_three_component_positive_raw2','size'),n_NK_three_component_positive=('NK_three_component_positive_raw2','sum'),n_joint_top95=('Mast_NK_joint_review','sum')).reset_index();rows.to_csv(run/'stage_f_Mast_expanded_context_v1/NK_joint_cluster_counts.tsv',sep='\t',index=False)
out={'status':'NK_joint_evidence_profile_not_yet_excluded','UTC':datetime.now(timezone.utc).isoformat(),'n_cells':len(q),'n_three_component_positive':int(q.NK_three_component_positive_raw2.sum()),'n_joint_top95':int(q.Mast_NK_joint_review.sum()),'rule':'At least one KLRD1/NCR1 receptor + at least one GNLY/PRF1 effector + NKG7 each >=2 raw UMI, plus original same-unit doublet percentile >=.95. Parent-Mast model input already marker-supported. Project positive diagnostic, not validated exclusive classifier.','counts':rows.to_dict('records')};(run/'stage_f_Mast_expanded_context_v1/NK_joint_profile.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
