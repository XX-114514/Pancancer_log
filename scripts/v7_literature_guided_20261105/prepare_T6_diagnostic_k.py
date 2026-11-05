import os,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad,h5py
from anndata.io import read_elem
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';o=r/'stage_k_T6_source_robust_diagnostic_v1';o.mkdir(exist_ok=False)
assert (r/'stage_j_T6_fullgene/export_complete.json').exists()
a=ad.read_h5ad(r/'stage_j_T6_fullgene/T6_v7_fullgene_raw.h5ad');assert a.shape==(156631,95038)
v=pd.read_pickle(ROOT/'runs/20260818_pancancer_5M_merge_prep/11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz');expected=v.index[v.cluster_uid.astype(str).eq('Tcell_6')&v.major_celltype_v7.eq('Tcell')];assert a.obs_names.equals(expected)
av=[];gses=sorted(a.obs.gse_id.astype(str).unique())
for gs in gses:
 p=ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad')
 with h5py.File(p,'r') as f:var=pd.Index(read_elem(f['var']).index.astype(str))
 av.append(a.var_names.isin(var))
av=np.asarray(av).T;coverage=av.mean(axis=1);allowed=(coverage>=.8)&~a.var_names.str.match('^HB[AB]');assert allowed.sum()>10000
pd.DataFrame({'gene_symbol':a.var_names[allowed]}).to_csv(o/'source_robust_model_genes.tsv',sep='\t',index=False)
manifest={'status':'T6_input_scope_verified_and_model_features_prepared','UTC':datetime.now(timezone.utc).isoformat(),'n_cells':a.n_obs,'n_genes':a.n_vars,'n_sources':len(gses),'n_source_robust_model_genes':int(allowed.sum()),'exact_frozen_V7_cluster_uid':'Tcell_6','exact_frozen_V7_major':'Tcell','frozen_scope_index_exact':True,'legacy_template_input_scope_description':str(a.uns['provenance'].get('input_scope','')),'scope_correction':'The prior generic exporter wording used major==Tcell_6. The actual executed mask was cluster_uid==Tcell_6 AND major==Tcell; exact full-rowset equality is now independently verified. Raw matrix is unchanged.','full_input_policy':'All 156631 frozen T6 cells retained in the initial diagnostic; per-cluster QC/competitor review follows. No Treg reference prefilter.','model_feature_policy':'Measured in >=80% actual T6 sources; HB omitted from HVG candidates only. Full RNA assay retained.','new_final_label_changes':0}
(o/'input_complete.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(manifest))
