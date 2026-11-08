import json,sys,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from scipy.spatial.distance import cdist
r=Path(sys.argv[1]);o=r/'spatial_pilot_v1';d=pd.read_csv(o/'identity_candidates_PRIVATE.tsv.gz',sep='\t');rec=json.loads((o/'complete.json').read_text());assert len(d)==690322==rec['cells_input'];assert d.QC_keep.sum()==rec['QC_pass']
q=d.loc[d.QC_keep];ec=q.loc[q.identity_candidate.eq('Blood_EC'),['x_centroid','y_centroid']].to_numpy();err=[]
from scipy.spatial import cKDTree
for typ in ['B','NK','T']:
 xy=q.loc[q.identity_candidate.eq(typ),['x_centroid','y_centroid']].to_numpy();pick=np.random.default_rng(52).choice(len(xy),min(128,len(xy)),replace=False);a=cdist(xy[pick],ec).min(1);b=cKDTree(ec).query(xy[pick])[0];err.append(float(np.max(np.abs(a-b))));assert np.allclose(a,b,rtol=0,atol=1e-8)
 for row in rec['distance']:
  if row['identity_candidate']==typ:assert len(xy)==row['n'];assert abs(float(np.median(cKDTree(ec).query(xy)[0]))-row['median_nearest_EC_centroid_um'])<1e-8
res={'audit_pass':True,'checks':['690322 input cells','QC counts','384 independently brute-force Euclidean nearest distances','full B NK T distance median recomputation','source member CRC and SHA checked by producer'],'max_abs_distance_error_um':max(err),'scope':'software and arithmetic audit; identity and biological validation remain pending'}
(o/'independent_audit.json').write_text(json.dumps(res,indent=2));print(json.dumps(res))
