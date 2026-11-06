from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,json
from scipy.stats import ttest_ind,ttest_ind_from_stats
r=Path(__file__).resolve().parents[1];a=ad.read_h5ad(r/'inputs/global_1000.h5ad');assert a.obs_names.is_unique and a.var_names.is_unique and (a.obs.identity.value_counts()<=1000).all();assert np.isfinite(a.X.data).all() and a.X.data.min()>=0
rng=np.random.default_rng(1337);x=rng.normal(size=(47,103));y=rng.normal(size=(63,103));p1=ttest_ind(x,y,equal_var=False).pvalue;p2=ttest_ind_from_stats(x.mean(0),x.std(0,ddof=1),len(x),y.mean(0),y.std(0,ddof=1),len(y),equal_var=False).pvalue;assert np.max(np.abs(p1-p2))<1e-12
# Explicit CMT scale and tie invariants.
h=np.array([[3,1,1],[1,4,1]],float);assert np.array_equal(h.argmax(0),(h/h.sum(0)).argmax(0));assert np.isclose(h,h.max(0)).sum(0).tolist()==[1,1,2]
print(json.dumps({'status':'PASS','global_shape':a.shape,'max_identity_cells':int(a.obs.identity.value_counts().max()),'Welch_stats_max_error':float(np.max(np.abs(p1-p2))),'CMT_tie_and_column_scaling':'PASS'}))
