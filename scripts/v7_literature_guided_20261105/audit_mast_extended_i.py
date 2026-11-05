import os,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad,zarr,h5py
from anndata.io import read_elem
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';o=r/'stage_i_mast_extended_review_v1';assert not (o/'independent_audit.json').exists()
a=ad.read_h5ad(r/'stage_b_mast_fullgene/mast_v7_fullgene_raw.h5ad');q=pd.read_pickle(o/'Mast_annotation_review_v1_PRIVATE.pkl.gz');assert q.index.equals(a.obs_names)
for col in ['major_celltype_v7','final_annotation_v7','annotation_confidence_v7']:assert q[col].astype(str).equals(a.obs[col].astype(str))
# Independent complete-sample empirical CDF, not Pandas group rank from producer.
core=zarr.open_group(str(ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');pos=q.global_position.to_numpy().astype(int);assert pd.Index(read_elem(core['obs/_index']))[pos].equals(q.index)
c=pd.DataFrame({k:read_elem(core['obs/'+k]) for k in ['gse_id','sampleID','doublet_method','doublet_score']});c.doublet_score=pd.to_numeric(c.doublet_score,errors='coerce')
keys=c[['gse_id','sampleID','doublet_method']].astype(str);codes,_=pd.factorize(pd.MultiIndex.from_frame(keys),sort=False);group_indices=pd.Series(np.arange(len(c))).groupby(codes).indices;tc=codes[pos];scores=c.doublet_score.to_numpy();recalc=np.full(len(q),np.nan);sizes=np.zeros(len(q),int)
for code in np.unique(tc):
 members=group_indices[int(code)];values=np.sort(scores[members][np.isfinite(scores[members])]);ii=np.flatnonzero(tc==code);v=scores[pos[ii]];valid=np.isfinite(v);sizes[ii]=len(values)
 if len(values):recalc[ii[valid]]=(np.searchsorted(values,v[valid],side='left')+np.searchsorted(values,v[valid],side='right')+1)/(2*len(values))
assert np.allclose(recalc,q.doublet_within_sample_method_percentile.to_numpy(),rtol=0,atol=1e-12,equal_nan=True);assert np.array_equal(sizes,q.doublet_calibration_n.to_numpy())
unknown=q.sampleID.isna()|q.sampleID.astype(str).str.lower().isin(['','nan','none','unknown','na','not_available'])
# Independent sparse marker logic + exact policy/readback counts.
def vals(gs):return a.X[:,a.var_names.get_indexer(gs)].toarray()
parent=(vals(['TPSAB1','TPSB2'])>=2).any(axis=1)&(vals(['CPA3','KIT','MS4A2'])>=2).any(axis=1);assert np.array_equal(parent,q.Mast_tryptase_plus_support_raw2)
panels={'T':['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2'],'B':['MS4A1','CD79A','CD79B','CD19','BANK1'],'Plasma':['MZB1','JCHAIN','SDC1','DERL3'],'Myeloid':['LST1','CSF1R','CD14','CD163','MS4A7','C1QA'],'Stromal':['COL1A1','COL1A2','DCN','LUM'],'Epithelial':['EPCAM','KRT8','KRT18','KRT19','KRT14','KRT17'],'Endothelial':['PECAM1','VWF','CLDN5','KDR'],'Platelet':['PF4','PPBP','GP9','ITGA2B']}
comp=np.zeros(len(q),bool)
for name,gs in panels.items():
 x=vals(gs);flag=((x>0).sum(axis=1)>=3)|((x>=2).sum(axis=1)>=2);assert np.array_equal(flag,q['extended_'+name+'_program_review']);comp|=flag
joint=parent&comp&(sizes>=100)&np.isfinite(recalc)&(recalc>=.95);clear=parent&~comp
assert int(joint.sum())==1323 and int(clear.sum())==38959;assert np.array_equal(clear,q.reviewed_major_identity_v1.eq('Mast'))
md=pd.read_csv(r/'stage_h_mast_source_robust_refit_v1/refit/metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');fit=q.Mast_source_robust_refit_candidate_v2.to_numpy();assert md.index.equals(q.index[fit]);cl=md['source_snn_res.0.4'].to_numpy().astype(int)
headers={}
for gs in q.gse_id.astype(str).unique():
 with h5py.File(ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad'),'r') as f:headers[gs]=set(read_elem(f['var']).index.astype(str))
e=pd.read_csv(o/'Mast_H15_expanded_measured_markers.tsv',sep='\t');tot=np.asarray(a.X[fit].sum(axis=1)).ravel();checked=0
for gene,s in e.groupby('gene'):
 v=a.X[fit,a.var_names.get_loc(gene)].toarray().ravel();measured=q.iloc[np.flatnonzero(fit)].gse_id.astype(str).map(lambda x:gene in headers[x]).to_numpy()
 for row in s.to_dict('records'):
  ii=cl==row['cluster'];den=int((measured&ii).sum());assert den==row['n_measured'];assert int(ii.sum())==row['n_cells']
  if den:
   assert np.isclose(((v>0)&ii).sum()/den,row['pct_expressing_measured'],atol=1e-11,rtol=1e-10);assert np.isclose((10000*v[ii]/np.maximum(tot[ii],1)).sum()/den,row['mean_cp10k_measured'],rtol=2e-6,atol=1e-7)
  checked+=1
hvg=pd.read_csv(r/'stage_h_mast_source_robust_refit_v1/refit/HVG_source_robust.tsv',sep='\t');coverage=pd.read_csv(r/'stage_g_mast_marker_QC_v2/source_feature_coverage_PRIVATE.tsv',sep='\t').set_index('gene_symbol');assert len(hvg)==2000 and hvg.gene_symbol.is_unique;assert (coverage.loc[hvg.gene_symbol,'fraction_Mast_sources_measured']>=.8).all()
trans=pd.read_csv(o/'old_V7_to_Mast_reviewed_transition.tsv',sep='\t');assert trans.n_cells.sum()==45036
result={'status':'PASS_extended_Mast_review_numeric_and_scope_audit','UTC':datetime.now(timezone.utc).isoformat(),'n_cells':len(q),'n_sample_method_percentiles_independently_recomputed':len(q),'n_Mast_cells_with_missing_or_placeholder_recorded_sampleID':int(unknown.sum()),'n_extended_parent_programs_checked':len(panels),'n_marker_cluster_records_recomputed':checked,'n_HVG_source_eligibility_checked':len(hvg),'n_marker_supported_Mast_coarse':int(clear.sum()),'n_joint_doublet_candidates':int(joint.sum()),'original_V7_label_fields_exact':True,'scope':'Numeric and marker-evidence review acceptance only; no gold-standard doublet, purity, protein subtype or publication accuracy claim. Recorded sample units are not asserted independent patients.'}
(o/'independent_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
