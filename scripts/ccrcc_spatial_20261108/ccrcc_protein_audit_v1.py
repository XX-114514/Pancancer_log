from pathlib import Path
import json,sys
import h5py,numpy as np,pandas as pd
from scipy.sparse import csc_matrix
r=Path(sys.argv[1]);out=r/'spatial_pilot_v1';c=pd.read_csv(out/'identity_candidates_PRIVATE.tsv.gz',sep='\t',index_col=0)
with h5py.File(r/'inputs/member_000_PRIVATE.h5') as f:
 m=f['matrix'];x=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr();g=[z.decode() for z in m['features']['name'][:]];ft=[z.decode() for z in m['features']['feature_type'][:]];ids=[z.decode() for z in m['barcodes'][:]]
assert list(c.index)==ids
ix=[i for i,z in enumerate(ft) if z=='Protein Expression'];p=x[:,ix].toarray();names=[g[i] for i in ix];print(names)
rows=[]
for k in ['B','T','NK','Blood_EC','Lymphatic','Mural','Other']:
 sel=c.QC_keep&(c.identity_candidate==k)
 for j,name in enumerate(names):
  rows.append({'RNA_identity_candidate':k,'n':int(sel.sum()),'protein':name,'median':float(np.median(p[sel,j])),'mean_log1p':float(np.log1p(p[sel,j]).mean())})
pd.DataFrame(rows).to_csv(out/'independent_protein_review.tsv',sep='\t',index=False)
# Orthogonal evidence only. No data-driven thresholds or identity changes.
q=pd.DataFrame(rows);print(q[q.protein.str.contains('CD20|CD31|CD3|CD16|Granzyme',case=False)].to_json(orient='records'))
