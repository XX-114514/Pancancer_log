import os
from pathlib import Path
PROJECT_ROOT=Path(os.environ["PROJECT_ROOT"])
import os
os.environ["REVIEW_LINEAGE"]='Mast'
import os,json
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad,h5py
from anndata.io import read_elem
r=PROJECT_ROOT;run=r/'runs/20261106_v7_QC_clean_reannotation_v1';name=os.environ['REVIEW_LINEAGE'];d=run/'stage_i_Mast_NK_clean_refit_v2';e=run/'stage_j_Mast_NK_clean_marker_review_v3';terminal=json.loads((run/'launch'/('Mast_v2_refit_terminal.json')).read_text());assert terminal['exit_code']==0
cfg=json.loads((run/'stage_g_Mast_NK_cleanup_v2/refit_config.json').read_text());a=ad.read_h5ad(run/'stage_g_Mast_NK_cleanup_v2/Mast_clean_v2_raw.h5ad');md=pd.read_csv(d/'metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');assert a.obs_names.equals(md.index) and len(md)==cfg['n_cells'];q=pd.read_pickle(run/'stage_g_Mast_NK_cleanup_v2/Mast_cleaning_decisions_PRIVATE.pkl.gz').loc[md.index];assert q.primary_QC_keep.all() and q.clean_lineage_model_eligible.all()
cl=md['source_snn_res.0.4'].astype(int);counts=cl.value_counts().sort_index();res=pd.read_csv(d/'cluster_resolution_counts.tsv',sep='\t')
for graph,t in res.groupby('graph'):
 actual=md[graph].value_counts().sort_index();saved=t.set_index('cluster').n_cells.sort_index();assert actual.to_dict()==saved.to_dict() and t.n_cells.sum()==len(md)
src=pd.read_csv(d/'cluster_source_concentration.tsv',sep='\t')
for graph,t in src.groupby('graph'):
 for z in t.to_dict('records'):
  ss=md.loc[md[graph].eq(z['cluster']),'gse_id'].value_counts();assert len(ss)==z['n_gse'] and int(ss.sum())==z['n_cells'] and abs(ss.max()/ss.sum()-z['max_gse_fraction'])<1e-10
common=set(pd.read_csv(run/'stage_c_clean_views'/name/'DE_measured_all_sources_genes.tsv',sep='\t').gene_symbol);de=pd.read_csv(d/'RNA_FindAllMarkers_common_measured.tsv',sep='\t');assert set(de.gene)<=common and set(de.cluster.astype(int))<=set(counts.index)
assert np.isfinite(de[['avg_log2FC','p_val','p_val_adj','pct.1','pct.2']].to_numpy()).all() and de.p_val.between(0,1).all() and de.p_val_adj.between(0,1).all();assert not de.duplicated(['cluster','gene']).any()
hvg=pd.read_csv(d/'HVG_source_robust.tsv',sep='\t');assert len(hvg)==2000 and hvg.gene_symbol.is_unique and set(hvg.gene_symbol)<=set(pd.read_csv(run/'stage_c_clean_views'/name/'model_genes.tsv',sep='\t').gene_symbol)
m=pd.read_csv(e/'measured_marker_expression.tsv',sep='\t');assert not m.duplicated(['cluster','gene']).any();total=np.asarray(a.X.sum(axis=1)).ravel();available={}
for gs in md.gse_id.unique():
 with h5py.File(r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(str(gs)+'_symbol_counts.h5ad'),'r') as f:available[gs]=set(read_elem(f['var']).index.astype(str))
maxerr=0.0
for gene,t in m.groupby('gene'):
 ix=a.var_names.get_indexer([gene])[0];x=np.asarray(a.X[:,ix].toarray()).ravel() if ix>=0 else np.zeros(a.n_obs);measured=md.gse_id.map(lambda g:gene in available[g]).to_numpy();assert not ((x>0)&~measured).any()
 for row in t.to_dict('records'):
  take=cl.eq(row['cluster']).to_numpy();den=int((take&measured).sum());assert int(take.sum())==row['n_cells'] and den==row['n_measured'];nz=int((take&(x>0)).sum());assert nz==row['n_expressing']
  if den:
   mean=float((10000*x[take]/np.maximum(total[take],1)).sum()/den);err=abs(mean-row['mean_cp10k_measured']);maxerr=max(maxerr,err);assert err<1e-7 and abs(nz/den-row['pct_expressing_measured'])<1e-12
no_de=sorted(set(counts.index)-set(de.loc[(de.avg_log2FC>0)&(de.p_val_adj<.05),'cluster'].astype(int)))
out={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'lineage':name,'n_cells':len(md),'n_clusters_res0_4':len(counts),'n_marker_records':len(m),'max_absolute_mean_CP10K_recalculation_error':maxerr,'all_model_cells_exactly_QC_eligible':True,'all_resolution_and_source_counts_recomputed':True,'DE_restricted_to_common_measured_genes':True,'HVG2000_subset_of_source_robust_candidates':True,'clusters_no_positive_DE_padj05_under_current_settings':no_de,'final_fine_annotation_accepted':False};(e/'independent_postfit_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
