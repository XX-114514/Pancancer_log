import os
from pathlib import Path
PROJECT_ROOT=Path(os.environ["PROJECT_ROOT"])
import json,hashlib
import numpy as np,pandas as pd,anndata as ad,h5py
r=PROJECT_ROOT;run=r/'runs/20261106_v7_QC_clean_reannotation_v1';c=run/'stage_c_clean_views'
assert (c/'complete.json').exists(), 'C not complete; observe, do not rerun'
v=pd.read_pickle(r/'runs/20260818_pancancer_5M_merge_prep/11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz')
mask=np.load(c/'global_QC_view_masks_PRIVATE.npz');ex=np.zeros(len(v),bool);hold=ex.copy();p=json.loads((run/'QC_POLICY_v1.json').read_text());approved=set(p['accepted_supplement_positive_units'])
for name in ['T6','Mast']:
 q=pd.read_pickle(run/'stage_a_evidence'/(name+'_adjudication_evidence_PRIVATE.pkl.gz'));pos=v.index.get_indexer(q.index);j=q.joint_review_original.to_numpy();a=q.distinct_competing_anchor_supported.to_numpy();ex[pos[j&a]]=True;hold[pos[j&~a]]=True
b=pd.read_pickle(run/'stage_b_unscored_supplement/unscored_supplement_PRIVATE.pkl.gz')
for (gs,sample,method),z in b.groupby(['gse_id','sampleID','doublet_method'],observed=True,dropna=False):
 h=hashlib.sha256((str(gs)+'|'+str((sample,method))).encode()).hexdigest();ix=v.index.get_indexer(z.index[z.supplement_operational_exclude]);assert (ix>=0).all()
 if h in approved:ex[ix]=True
 else:hold[ix]=True
assert np.array_equal(ex,mask['primary_exclude']) and np.array_equal(hold,mask['QC_holdout'])
assert not (ex&hold).any() and np.array_equal(~ex&~hold,mask['primary_keep'])
assert (int(ex.sum()),int(hold.sum()),int(mask['primary_keep'].sum()))==(1967,236,4674584)
ih=hashlib.sha256(pd.util.hash_pandas_object(v.index,index=False).to_numpy().tobytes()).hexdigest();assert ih=='3f802eb8cc52f04beb27fa2987d2ffe3ddf75f4b4f51dec8e6448b26df8d1f16'
rows=[];lineages=[];old=r/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1'
for name in ['T6','Mast']:
 q=pd.read_pickle(c/(name+'_cleaning_decisions_PRIVATE.pkl.gz'));pos=v.index.get_indexer(q.index);assert np.array_equal(q.primary_QC_keep.to_numpy(),(~ex&~hold)[pos]);eligible=q.primary_QC_keep.to_numpy().copy()
 if name=='Mast':eligible&=q.review_route_v4.eq('marker_supported_Mast_coarse').to_numpy()
 assert np.array_equal(eligible,q.clean_lineage_model_eligible.to_numpy())
 raw=ad.read_h5ad(old/('stage_j_T6_fullgene/T6_v7_fullgene_raw.h5ad' if name=='T6' else 'stage_b_mast_fullgene/mast_v7_fullgene_raw.h5ad'));clean=ad.read_h5ad(c/name/(name+'_clean_raw.h5ad'));assert raw.obs_names.equals(q.index) and clean.obs_names.equals(q.index[eligible]) and raw.var_names.equals(clean.var_names)
 X=raw.X[eligible].tocsr();Y=clean.X.tocsr();assert X.shape==Y.shape and np.array_equal(X.indptr,Y.indptr) and np.array_equal(X.indices,Y.indices) and np.array_equal(X.data,Y.data)
 with h5py.File(c/name/(name+'_clean_raw_10x.h5'),'r') as f:
  assert tuple(f['matrix/shape'][:])==(Y.shape[1],Y.shape[0]);assert np.array_equal(f['matrix/data'][:],Y.data) and np.array_equal(f['matrix/indices'][:],Y.indices) and np.array_equal(f['matrix/indptr'][:],Y.indptr);assert np.array_equal(f['matrix/barcodes'].asstr()[:],clean.obs_names.to_numpy())
 cfg=json.loads((c/name/'refit_config.json').read_text());assert cfg['n_cells']==len(clean.obs) and cfg['n_genes']==95038
 for cl,t in q.groupby('review_cluster_before_cleanup',observed=True):
  pp=v.index.get_indexer(t.index);rows.append(dict(lineage=name,old_review_cluster=int(cl),n_before=len(t),n_exclude_primary=int(ex[pp].sum()),n_holdout_QC=int(hold[pp].sum()),n_primary_QC_keep=int((~ex&~hold)[pp].sum()),n_lineage_fit=int(t.clean_lineage_model_eligible.sum())))
 lineages.append(dict(lineage=name,n_before=len(q),n_exclude=int(ex[pos].sum()),n_holdout=int(hold[pos].sum()),n_fit=len(clean.obs),n_sources=cfg['n_sources'],model_gene_candidates=cfg['n_model_candidates'],DE_common_measured_genes=cfg['n_global_DE_candidates'],exact_counts=True))
pd.DataFrame(rows).to_csv(c/'old_cluster_cleanup_counts_audited.tsv',sep='\t',index=False)
out={'status':'PASS','frozen_index_hash':ih,'n_frozen':len(v),'n_excluded':int(ex.sum()),'n_holdout':int(hold.sum()),'n_keep':int((~ex&~hold).sum()),'lineages':lineages,'count_table_correction':'Original old_cluster_cleanup_counts.tsv holdout count used a mismatched label. Audited table uses exact independent masks. No mask or input error.'}
(c/'independent_input_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
