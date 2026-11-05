"""Independent multiset, flag and original pre-filter-ledger checks. Read only."""
import json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,h5py,zarr
from anndata.io import read_elem
r=PROJECT_ROOT;b=r/'runs/20260818_pancancer_5M_merge_prep';u=r/'runs/20260725_120105_unified_pancancer_covarnet';g=zarr.open_group(str(r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');d=pd.DataFrame({c:read_elem(g['obs/'+c]) for c in ['gse_id','source_obs_name','cohort','doublet_score','doublet_status']});old=pd.read_csv(u/'06_full_qc_integration/phase06_production_qc_validation_33datasets_20260726_051121.tsv',sep='\t');new=pd.read_csv(b/'01_new_cohort_prep/manifests/dataset_qc_summary.tsv',sep='\t')
def resolve(x):
 p=Path(str(x));return p if p.is_absolute() else r.parent/p if str(p).startswith('Project_v3/') else r/p
nchecked=0;nfinite=0;old_raw_rows=0;old_doublet_removed=0;old_doublet_only=0;old_selected=[]
for gs,t in d.groupby('gse_id',observed=True):
 isold=t.cohort.astype(str).eq('existing_2M').all();p=resolve(old.loc[old.gse_id.eq(gs),'h5ad_output'].iloc[0]) if isold else resolve(new.loc[new.gse_id.eq(gs),'output_path'].iloc[0]);names=t.source_obs_name.astype(str).to_numpy();assert len(np.unique(names))==len(t)
 with h5py.File(p,'r') as f:
  ob=f['obs'];source=np.asarray(read_elem(ob[ob.attrs['_index']])).astype(str);assert len(source)==len(t) and len(np.unique(source))==len(source);assert np.array_equal(np.sort(source),np.sort(names))
  for c in ['predicted_doublet','doublet_removed']:
   assert c in ob;v=pd.Series(read_elem(ob[c])).astype(str).str.lower();assert v.isin(['false','0']).all()
  score=pd.Series(read_elem(ob['doublet_score']),index=source,dtype='float64').reindex(names).to_numpy();np.testing.assert_allclose(score,pd.to_numeric(t.doublet_score,errors='coerce').to_numpy(),rtol=0,atol=0,equal_nan=True);nfinite+=int(np.isfinite(score).sum());nchecked+=len(t)
 if isold:
  p=u/'06_full_qc_integration/per_dataset'/str(gs);files=list(p.glob('*cell_qc.tsv'));assert len(files)==1;head=pd.read_csv(files[0],sep='\t',nrows=0);cols=[c for c in ['predicted_doublet','doublet_removed','filter_reason'] if c in head];assert len(cols)==3;raw=pd.read_csv(files[0],sep='\t',usecols=cols);removed=raw.doublet_removed.astype(str).str.lower().isin(['true','1']);reasons=raw.filter_reason.astype(str);assert np.array_equal(removed.to_numpy(),reasons.str.contains('doublet').to_numpy());v=old[old.gse_id.eq(gs)].iloc[0];assert len(raw)==int(v.cells_pre_qc) and int(reasons.eq('pass').sum())==len(t);old_raw_rows+=len(raw);old_doublet_removed+=int(removed.sum());old_doublet_only+=int(reasons.eq('doublet').sum());old_selected.append(str(gs))
thresholds=pd.read_csv(b/'01_new_cohort_prep/manifests/sample_qc_thresholds.tsv',sep='\t');assert int(thresholds.n_doublets_removed.sum())==13821 and int(thresholds.n_predicted_doublets.sum())==13821
assert nchecked==4676787 and nfinite==4647788 and old_raw_rows==2148018 and old_doublet_removed==5518 and old_doublet_only==5177
result={'status':'PASS_independent_V7_doublet_provenance_audit','observed_UTC':datetime.now(timezone.utc).isoformat(),'n_current_cells_multiset_exact_with_QC_objects':nchecked,'n_current_sources':54,'all_source_cell_keys_unique':True,'n_current_cells_score_exact_including_missing':nchecked,'n_current_finite_historical_scores':nfinite,'n_current_missing_scores':nchecked-nfinite,'n_current_historical_predicted_positive':0,'n_current_historical_removed_flags':0,'n_old_raw_cell_ledger_rows_checked':old_raw_rows,'n_old_removed_by_doublet_rule':old_doublet_removed,'n_old_doublet_only_exclusions':old_doublet_only,'n_new_doublet_removals_reconciled':13821,'total_doublet_rule_exclusions_including_old_QC_overlap':old_doublet_removed+13821,'limits':'Verifies historical execution, filtering and membership. Does not validate biological singlet purity, independently rerun Scrublet, or assert all residual mixed RNA is a doublet.'}
print(json.dumps(result))
