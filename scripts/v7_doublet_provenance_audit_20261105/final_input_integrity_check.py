import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd,zarr
from anndata.io import read_elem
r=PROJECT_ROOT;p=r/'runs/20260818_pancancer_5M_merge_prep/11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz';h=hashlib.sha256()
with p.open('rb') as f:
 for part in iter(lambda:f.read(8*1024*1024),b''):h.update(part)
assert h.hexdigest()=='e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e'
g=zarr.open_group(str(r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r');v=pd.to_numeric(pd.Series(read_elem(g['obs/doublet_score'])),errors='coerce').to_numpy();status=pd.Series(read_elem(g['obs/doublet_status'])).astype(str).to_numpy();assert np.array_equal(np.isnan(v),status=='not_run') and not np.isinf(v).any()
print(json.dumps({'status':'PASS_frozen_input_and_full_score_status_integrity','observed_UTC':datetime.now(timezone.utc).isoformat(),'frozen_checkpoint_sha256':h.hexdigest(),'n_cells':len(v),'n_finite_scores':int(np.isfinite(v).sum()),'n_missing_scores':int(np.isnan(v).sum()),'n_infinite_scores':int(np.isinf(v).sum()),'missing_score_exactly_equals_not_run_status':True,'n_status_failure':int(pd.Series(status).str.contains('fail',case=False).sum()),'finite_min':float(np.nanmin(v)),'finite_max':float(np.nanmax(v))}))
