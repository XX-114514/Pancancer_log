"""Same-V4 sample pseudobulk and complete observed imaging input inventory."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib
import numpy as np,pandas as pd,scipy.sparse as sp,anndata as ad,h5py
from scipy.io import mmwrite
from anndata.io import read_elem
r=PROJECT_ROOT/'runs/20261107_v4_spatialCellChat_bulk_v1';out=r/'bulk/signature_inputs';out.mkdir(exist_ok=False)
a=ad.read_h5ad(PROJECT_ROOT/'runs/20261107_V4_CoVarNet_communication_v1/inputs/full_common_PRIVATE.h5ad');assert a.shape==(560871,14051)
meta=a.obs[['sampleID','gse_id','CMT']].drop_duplicates();assert len(meta)==164 and meta.sampleID.is_unique
meta=meta.set_index('sampleID').sort_index();ind=meta.index.get_indexer(a.obs.sampleID);G=sp.csr_matrix((np.ones(a.n_obs),(ind,np.arange(a.n_obs))),shape=(164,a.n_obs));linear=a.X.copy();linear.data=np.expm1(linear.data);mean=(G@linear).toarray()/np.bincount(ind)[:,None];assert np.isfinite(mean).all() and (mean>=0).all()
mmwrite(out/'mean_CP10K_PRIVATE.mtx',sp.csr_matrix(mean.T));pd.Series(a.var_names).to_csv(out/'genes.tsv',sep='\t',index=False,header=False);meta.to_csv(out/'samples_PRIVATE.tsv',sep='\t');pd.DataFrame({'CMT':meta.CMT.value_counts().index,'n_samples':meta.CMT.value_counts().values}).to_csv(out/'CMT_signature_sample_scope.tsv',sep='\t',index=False)
p=PROJECT_ROOT.parent/'Project_v5';jobs=[]
for path in sorted((p/'spatial_expansion_20260903/processed/XENIUM_BREAST_10X').glob('*.h5ad')):
 jobs.append({'dataset':'XENIUM_BREAST_10X','platform':'Xenium','format':'h5ad','path':str(path),'coordinate_unit':'micrometre'})
for path in sorted((p/'xenium_cosmx_20260904/raw/xenium_io').glob('*/*cell_feature_matrix.h5')):
 jobs.append({'dataset':path.parent.name,'platform':'Xenium','format':'10x_h5','path':str(path),'coordinates_path':str(next(path.parent.glob('*_cells.csv.gz'))),'coordinate_unit':'micrometre'})
for path in sorted((p/'xenium_cosmx_20260904/raw/cosmx_nsclc').glob('*_flat.tar.gz')):
 jobs.append({'dataset':'COSMX_NSCLC_He2022','platform':'CosMx','format':'flat_tar','path':str(path),'coordinate_unit':'global_pixel','pixel_um':0.12028})
assert len(jobs)==13
for i,z in enumerate(jobs):z['unit_code']='I%03d'%(i+1);z['input_bytes']=Path(z['path']).stat().st_size
(r/'config/imaging_jobs_PRIVATE.json').open('x').write(json.dumps(jobs,indent=2));receipt={'status':'SAME164_SAMPLE_CMT_PSEUDOBULK_AND_ALL13_IMAGING_INPUTS_READY','UTC':datetime.now(timezone.utc).isoformat(),'n_samples':164,'n_genes':14051,'n_CM':14,'imaging_units':13,'Xenium':5,'CosMx':8,'no_clinical_outcomes_used':True,'score_not_exact_NMF_activity':True}
(out/'complete.json').open('x').write(json.dumps(receipt,indent=2));print(json.dumps(receipt),flush=True)
