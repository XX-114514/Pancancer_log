"""Independent arithmetic/scope checks; not biological validation."""
import os,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad,zarr,h5py
from anndata.io import read_elem
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';o=r/'stage_l2_T6_marker_review_v1';dest=o/'independent_audit.json';assert not dest.exists()
obs=pd.read_pickle(r/'stage_l_T6_marker_review_v1/T6_marker_review_PRIVATE.pkl.gz');a=ad.read_h5ad(r/'stage_j_T6_fullgene/T6_v7_fullgene_raw.h5ad');assert obs.index.equals(a.obs_names) and len(obs)==156631
for c in a.obs.columns:pd.testing.assert_series_equal(obs[c],a.obs[c])
e=pd.read_csv(o/'T6_measured_marker_expression.tsv',sep='\t');genes=list(e.gene.unique());gses=sorted(obs.gse_id.astype(str).unique());avail={};hvg=pd.read_csv(r/'stage_k_T6_source_robust_diagnostic_v1/fit/HVG_source_robust.tsv',sep='\t').gene_symbol
hv=[]
for gs in gses:
 with h5py.File(ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad'),'r') as f:var=pd.Index(read_elem(f['var']).index.astype(str))
 avail[gs]=set(g for g in genes if g in var);hv.append(hvg.isin(var).to_numpy())
assert (np.asarray(hv).mean(axis=0)>=.8).all() and len(hvg)==2000
raw=a.X[:,a.var_names.get_indexer(genes)].toarray();total=np.asarray(a.X.sum(axis=1)).ravel();cl=obs['T6_review_cluster_source_res0.4'].to_numpy();gse=obs.gse_id.astype(str).to_numpy();checked=0;measurement={gene:np.isin(gse,[gs for gs in gses if gene in avail[gs]]) for gene in genes}
for c in np.unique(cl):
 for j,gene in enumerate(genes):
  take=(cl==c)&measurement[gene];v=raw[take,j];row=e[(e.cluster==c)&(e.gene==gene)].iloc[0];den=int(take.sum());assert den==row.n_measured and int((v>0).sum())==row.n_expressing
  if den:
   np.testing.assert_allclose(float((v*10000/total[take]).mean()),row.mean_cp10k_measured,rtol=1e-9,atol=1e-9);np.testing.assert_allclose(float((v>0).mean()),row.pct_expressing_measured,rtol=1e-9,atol=1e-9)
  checked+=1
core=zarr.open_group(str(ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');pos=obs.global_position.to_numpy().astype(int);assert pd.Index(read_elem(core['obs/_index']))[pos].equals(obs.index)
cal=pd.DataFrame({c:read_elem(core['obs/'+c]) for c in ['gse_id','sampleID','doublet_method','doublet_score']});cal['doublet_score']=pd.to_numeric(cal.doublet_score,errors='coerce');inverse=np.full(len(cal),-1,np.int64);inverse[pos]=np.arange(len(pos));recal=np.full(len(obs),np.nan);sizes=np.zeros(len(obs),np.int64)
scores=cal.doublet_score.to_numpy()
for _,rows in cal.groupby(['gse_id','sampleID','doublet_method'],dropna=False,observed=True).indices.items():
 rows=np.asarray(rows);selected=inverse[rows];ok=(selected>=0)&np.isfinite(scores[rows]);slots=selected[ok]
 if not len(slots):continue
 vals=np.sort(scores[rows][np.isfinite(scores[rows])]);q=scores[rows[ok]];recal[slots]=(np.searchsorted(vals,q,side='left')+np.searchsorted(vals,q,side='right')+1)/(2*len(vals));sizes[slots]=len(vals)
np.testing.assert_allclose(recal,obs.doublet_within_sample_method_percentile.to_numpy(),rtol=1e-12,atol=1e-12,equal_nan=True);assert np.array_equal((sizes>=100)&np.isfinite(recal),obs.doublet_calibration_evaluable.to_numpy())
spec=json.loads((r/'config/T6_marker_panels_v1.json').read_text());ix={g:i for i,g in enumerate(genes)};parent=(raw[:,[ix[g] for g in ['CD3D','CD3E','CD3G']]]>=2).any(axis=1)&(raw[:,[ix[g] for g in ['TRAC','TRBC1','TRBC2']]]>=2).any(axis=1);assert np.array_equal(parent,obs.T_parent_cross_component_raw2.to_numpy());anycomp=np.zeros(len(obs),bool)
for p in spec['panels']:
 if p['level']!='competing_parent':continue
 x=raw[:,[ix[g] for g in p['genes']]];flag=((x>0).sum(axis=1)>=3)|((x>=2).sum(axis=1)>=2);assert np.array_equal(flag,obs[p['id']+'_program_review'].to_numpy());anycomp|=flag
joint=parent&anycomp&(sizes>=100)&np.isfinite(recal)&(recal>=.95);assert np.array_equal(joint,obs.joint_competitor_doublet_candidate.to_numpy())
summary=pd.read_csv(o/'T6_cluster_marker_QC_summary.tsv',sep='\t');assert summary.n_cells.sum()==156631 and summary.n_joint_doublet_candidates.sum()==joint.sum()
result={'status':'PASS_T6_marker_QC_numeric_scope_audit','UTC':datetime.now(timezone.utc).isoformat(),'n_cells':len(obs),'n_marker_records_recomputed':checked,'n_doublet_background_percentiles_recomputed':len(obs),'n_HVG_source_eligibility_checked':len(hvg),'n_competitor_programs_checked':8,'n_joint_doublet_candidates':int(joint.sum()),'old_export_fields_exact':True,'counts_conserved':True,'limits':'Arithmetic and evidence-scope acceptance only; not validated Treg identity, purity, doublet accuracy or source/resolution biological stability.'}
dest.write_text(json.dumps(result,indent=2));print(json.dumps(result))
