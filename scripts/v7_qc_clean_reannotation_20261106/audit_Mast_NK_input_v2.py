import os
from pathlib import Path
PROJECT_ROOT=Path(os.environ["PROJECT_ROOT"])
import json
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad,h5py
run=PROJECT_ROOT/'runs/20261106_v7_QC_clean_reannotation_v1';g=run/'stage_g_Mast_NK_cleanup_v2';assert (g/'complete.json').exists();a=ad.read_h5ad(run/'stage_c_clean_views/Mast/Mast_clean_raw.h5ad');q=pd.read_pickle(run/'stage_c_clean_views/Mast_cleaning_decisions_PRIVATE.pkl.gz').loc[a.obs_names]
gs=['KLRD1','NCR1','GNLY','PRF1','NKG7'];x=a.X[:,a.var_names.get_indexer(gs)].toarray();take=(x[:,:2]>=2).any(axis=1)&(x[:,2:4]>=2).any(axis=1)&(x[:,4]>=2)&q.doublet_calibration_evaluable.to_numpy()&(pd.to_numeric(q.doublet_within_sample_method_percentile,errors='coerce').to_numpy()>=.95)
old=np.load(run/'stage_c_clean_views/global_QC_view_masks_PRIVATE.npz');new=np.load(g/'global_QC_view_masks_PRIVATE.npz');expected=q.global_position.to_numpy().astype(int)[take];actual=np.flatnonzero(new['primary_exclude']&~old['primary_exclude']);assert np.array_equal(np.sort(expected),actual) and np.array_equal(old['QC_holdout'],new['QC_holdout']) and not (new['primary_exclude']&new['QC_holdout']).any();assert np.array_equal(new['primary_keep'],~new['primary_exclude']&~new['QC_holdout'])
b=ad.read_h5ad(g/'Mast_clean_v2_raw.h5ad');assert b.obs_names.equals(a.obs_names[~take]) and b.var_names.equals(a.var_names);X=a.X[~take].tocsr();assert np.array_equal(X.data,b.X.data) and np.array_equal(X.indices,b.X.indices) and np.array_equal(X.indptr,b.X.indptr)
assert (int(new['primary_exclude'].sum()),int(new['QC_holdout'].sum()),int(new['primary_keep'].sum()))==(2003,236,4674548)
out={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'n_fresh_raw_recomputed_additional_exclusions':int(take.sum()),'n_Mast_model':b.n_obs,'n_global_exclude':int(new['primary_exclude'].sum()),'n_global_holdout':int(new['QC_holdout'].sum()),'n_global_keep':int(new['primary_keep'].sum()),'exact_original_count_subset':True,'calibrated_sample_method_background_required':True};(g/'independent_input_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
