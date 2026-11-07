"""Independent membership, CMT and native-normalization readback before CPDB."""
import json, hashlib
from datetime import datetime, timezone
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, zarr
from anndata.io import read_elem

r=PROJECT_ROOT/'runs/20261107_V4_CoVarNet_communication_v1'
base=PROJECT_ROOT/'runs/20261106_v7_selective_CoVarNet_paper_v2'
q=pd.read_pickle(r/'private/cell_bridge_PRIVATE.pkl.gz')
original=pd.read_pickle(base/'inputs/selected_cell_annotation_and_sample_PRIVATE.pkl.gz')
assert q.index.is_unique and len(q)==len(original)==560871 and set(q.index)==set(original.index)
for c in original.columns:assert np.array_equal(q[c].to_numpy(),original.loc[q.index,c].to_numpy()),c
H=pd.read_csv(base/'CM_results/H_PRIVATE.tsv',sep='\t',index_col=0)
assert q.CMT.equals(q.sampleID.map(H.idxmax()))
a=ad.read_h5ad(r/'inputs/full_common_PRIVATE.h5ad',backed='r')
b=ad.read_h5ad(r/'inputs/global_1000_PRIVATE.h5ad')
assert a.obs_names.equals(q.index) and b.obs_names.isin(a.obs_names).all()
assert set(a.obs.subCluster)==set(b.obs.subCluster)==set(original.subCluster)
assert b.obs.subCluster.value_counts().max()<=1000
g=zarr.open_group(str(PROJECT_ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r')
genes=pd.Index(read_elem(g['var/_index']));cols=genes.get_indexer(a.var_names);assert (cols>=0).all()
rng=np.random.default_rng(2026110702);checked=sorted(set(rng.choice(len(q),200,replace=False))|set(q.groupby('subCluster',sort=True).head(1).raw_row.map(lambda v:np.searchsorted(q.raw_row,v))))
maxerr=0
for ix in checked:
    row=int(q.raw_row.iloc[ix]);lo,hi=map(int,g['X/indptr'][row:row+2]);raw=sp.csr_matrix((g['X/data'][lo:hi],g['X/indices'][lo:hi],np.array([0,hi-lo])),shape=(1,len(genes)))
    expected=raw[:,cols].astype(float);expected.data=np.log1p(expected.data*10000/float(raw.sum()))
    diff=expected-a[ix].X;err=float(np.max(np.abs(diff.data))) if diff.nnz else 0;maxerr=max(maxerr,err)
assert maxerr<2e-6,maxerr
assert np.isfinite(b.X.data).all() and (b.X.data>=0).all()
eligible=pd.read_csv(r/'audit/CM_comparison_eligibility.tsv',sep='\t')
assert len(eligible)==14 and eligible.n_high_samples.sum()==164
proof={'status':'PASS_independent_sameV4_cells_labels_CMT_and_native_count_normalization',
    'UTC':datetime.now(timezone.utc).isoformat(),'n_cells':len(q),'n_labels':45,'n_samples':164,'n_CM':14,
    'n_independent_original_count_rows_checked':len(checked),'maximum_expression_error':maxerr,
    'n_common_measured_genes':a.n_vars,'n_global_cells':b.n_obs,'source_sample_as_patient_claim':False}
(r/'audit/input_independent_audit.json').open('x').write(json.dumps(proof,indent=2));print(json.dumps(proof))
