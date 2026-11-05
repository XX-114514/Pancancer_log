"""Independent arithmetic/source-mask readback of the new Mast review outputs."""
import os,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad,h5py
from anndata.io import read_elem
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';o=r/'stage_f_mast_review'
assert not (o/'independent_arithmetic_audit_v2.json').exists()
a=ad.read_h5ad(r/'stage_b_mast_fullgene/mast_v7_fullgene_raw.h5ad');X=a.X.tocsr();overlay=pd.read_pickle(o/'mast_cell_review_overlay_PRIVATE.pkl.gz')
assert overlay.index.equals(a.obs_names)
for c in ['major_celltype_v7','final_annotation_v7','annotation_confidence_v7']:assert overlay[c].astype(str).equals(a.obs[c].astype(str))
cl=overlay.mast_review_cluster_res0_4 if 'mast_review_cluster_res0_4' in overlay else overlay['mast_review_cluster_res0.4']
# Fresh source headers; full selected raw nonzero support must be source-measured.
header={};n_nonzeros_outside=0
for gs in a.obs.gse_id.astype(str).unique():
 p=ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad')
 with h5py.File(p,'r') as f:header[gs]=set(read_elem(f['var']).index.astype(str))
 ii=np.flatnonzero(a.obs.gse_id.astype(str).eq(gs).to_numpy());sub=X[ii];allowed=a.var_names.isin(header[gs]);n_nonzeros_outside+=int((~allowed[sub.indices]).sum())
assert n_nonzeros_outside==0
stats=pd.read_csv(o/'mast27_measured_marker_expression.tsv',sep='\t');n_checked=0
# Independent dense single-column + Pandas groupby calculation (not producer selector algebra).
for gene,s in stats.groupby('gene',sort=False):
 j=a.var_names.get_loc(gene);v=X[:,j].toarray().ravel();measured=a.obs.gse_id.astype(str).map(lambda gs:gene in header[gs]).to_numpy()
 for rec in s.to_dict('records'):
  take=cl.to_numpy()==rec['cluster'];den=int((measured&take).sum());pos=int(((v>0)&take).sum());assert den==rec['n_measured']
  if den:
   assert np.isclose(pos/den,rec['pct_expressing_measured'],rtol=1e-10,atol=1e-12)
   totals=np.asarray(X[take].sum(axis=1)).ravel();mean=float((10000*v[take]/np.maximum(totals,1)).sum()/den);assert np.isclose(mean,rec['mean_cp10k_measured'],rtol=2e-6,atol=1e-7)
  n_checked+=1
review=pd.read_csv(o/'mast27_cluster_review.tsv',sep='\t');assert len(review)==27 and review.n_cells.sum()==len(a)
for rec in review.to_dict('records'):
 take=cl.to_numpy()==rec['cluster'];g=a.obs.loc[take,'gse_id'].astype(str).value_counts();assert int(take.sum())==rec['n_cells'];assert len(g)==rec['n_gse'];assert np.isclose(float(g.max()/g.sum()),rec['max_gse_fraction'])
frozen=ROOT/'runs/20260818_pancancer_5M_merge_prep/11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz'
h=hashlib.sha256()
with frozen.open('rb') as f:
 for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
assert h.hexdigest()=='e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e'
result={'status':'PASS_complete_Mast_review_numeric_source_and_preservation_checks','UTC':datetime.now(timezone.utc).isoformat(),'n_cells_checked':len(a),'n_raw_nonzero_entries_source_mask_checked':int(X.nnz),'n_nonzeros_outside_measured_source_features':n_nonzeros_outside,'n_marker_cluster_records_independently_recomputed':n_checked,'n_cluster_count_source_records':len(review),'original_label_fields_exactly_preserved':True,'frozen_V7_SHA_unchanged':True,'scope':'Numeric/measurement/preservation acceptance; biological review routes remain provisional, no final subtype acceptance.'}
(o/'independent_arithmetic_audit_v2.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
