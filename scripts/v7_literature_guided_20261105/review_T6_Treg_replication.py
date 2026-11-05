import os,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad,zarr
from anndata.io import read_elem
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));B=ROOT/'runs/20260818_pancancer_5M_merge_prep'
r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';o=r/'stage_e_T6_Treg_candidate_review';o.mkdir(exist_ok=False)
df=pd.read_pickle(B/'11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz')
a=ad.read_h5ad(B/'07_fine_annotation/20260823_corrected_evidence_v1/objects/wide_marker_panel_full_library_masked_v2.h5ad');assert a.obs_names.equals(df.index)
assert list(a.varm['availability_by_gse'].columns.astype(str))==list(np.asarray(a.uns['availability_gse_order']).astype(str))
mask=df.major_celltype_v7.eq('Tcell');rows=np.flatnonzero(mask);obs=df.iloc[rows]
genes=['FOXP3','IKZF2','IL2RA','TNFRSF4','CTLA4','TIGIT'];X=a.X[rows][:,a.var_names.get_indexer(genes)].tocsr()
ix=obs.cluster_uid.astype(str).eq('Tcell_6').to_numpy();codes,units=pd.factorize(obs.gse_id.astype(str)+'::'+obs.sampleID.astype(str),sort=True);n=len(units)
H=sp.csr_matrix((np.ones(len(rows),np.int32),(codes,np.arange(len(rows)))),shape=(n,len(rows)))
inside=(H[:,ix]@X[ix]).toarray();outside=(H[:,~ix]@X[~ix]).toarray();ni=np.bincount(codes[ix],minlength=n);no=np.bincount(codes[~ix],minlength=n)
means_in=np.divide(inside,ni[:,None],out=np.zeros_like(inside),where=ni[:,None]>0);means_out=np.divide(outside,no[:,None],out=np.zeros_like(outside),where=no[:,None]>0)
# Match old sample support criterion for this formerly unexamined candidate.
avail=a.varm['availability_by_gse'].loc[genes];gse_for_unit=obs.groupby(codes,observed=True).gse_id.first().reindex(range(n)).astype(str)
measured=avail.loc[:,gse_for_unit.to_numpy()].to_numpy().T.astype(bool)
evaluable=(ni>=10)&(no>=10)&(measured.sum(axis=1)>=2)&(measured[:,:3].sum(axis=1)>=1)
delta=np.divide(((means_in-means_out)*measured).sum(axis=1),measured.sum(axis=1),out=np.zeros(n),where=measured.sum(axis=1)>0)
support=evaluable&(delta>.10)
co=np.asarray((X[ix,:3]>0).sum(axis=1)).ravel();t=sp.csr_matrix(a.X[rows[ix]][:,a.var_names.get_indexer(['CD3D','CD3E','TRAC'])]);tc=np.asarray((t>0).sum(axis=1)).ravel()
# Full-gene model labels are evidence only and never overwrite V7.
g=zarr.open_group(str(ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');assert pd.Index(read_elem(g['obs/_index'])).equals(df.index)
model=[]
for field in ['ev_celltypist_label','scatomic_full_label']:
 labels=np.asarray(read_elem(g['obs/'+field]))[rows[ix]]
 for label,count in pd.Series(labels).value_counts(dropna=False).head(15).items():model.append({'cluster_uid':'Tcell_6','field':field,'label':str(label),'n_cells':int(count)})
pd.DataFrame(model).to_csv(o/'T6_reference_model_aggregate.tsv',sep='\t',index=False)
result={'UTC':datetime.now(timezone.utc).isoformat(),'cluster_uid':'Tcell_6','n_cells':int(ix.sum()),'candidate':'CD4_Treg','legacy_full_panel_pass':6,'legacy_canonical_pass':3,'candidate_sample_n_evaluable':int(evaluable.sum()),'candidate_sample_n_support':int(support.sum()),'candidate_sample_support_fraction':float(support.sum()/evaluable.sum()),'n_ge2_Treg_canonical_genes':int((co>=2).sum()),'n_ge2_Treg_canonical_AND_ge2_T_identity':int(((co>=2)&(tc>=2)).sum()),'annotation_changes':0,'decision':'Treg-enriched candidate with within-cluster composition still to resolve; not an all-cell Treg assignment','sample_support_rule':'Within-sample T6 and other V7 T cells >=10; >=2 measured panel genes including >=1 canonical; mean log1p CP10K difference >0.10. Descriptive replication, not a hypothesis test.'}
(o/'complete.json').write_text(json.dumps(result,indent=2));print(json.dumps({'summary':result,'model_aggregate':model}))
