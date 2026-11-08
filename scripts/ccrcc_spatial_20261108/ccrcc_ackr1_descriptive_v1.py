"""Pre-specified ACKR1 context, descriptive only; no vessel segment reconstruction."""
import sys,json
from pathlib import Path
import numpy as np,pandas as pd,h5py
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree
r=Path(sys.argv[1]);o=r/'spatial_pilot_v1';c=pd.read_csv(o/'identity_candidates_PRIVATE.tsv.gz',sep='\t');qc=c.QC_keep.to_numpy();label=c.identity_candidate.to_numpy();xy=c[['x_centroid','y_centroid']].to_numpy()
with h5py.File(r/'inputs/member_000_PRIVATE.h5') as f:
 m=f['matrix'];x=csc_matrix((m['data'][:],m['indices'][:],m['indptr'][:]),shape=tuple(m['shape'][:])).T.tocsr();genes=[v.decode() for v in m['features']['name'][:]];ack=x[:,genes.index('ACKR1')].toarray().ravel()
ec=np.flatnonzero(qc&(label=='Blood_EC'));ii=np.flatnonzero(qc&~np.isin(label,['Blood_EC','Lymphatic']));dd,nearest=cKDTree(xy[ec]).query(xy[ii]);a=ack[ec[nearest]]>=1;rows=[]
for rad in [25,50,100,200]:
 for typ in ['B','NK','T']:
  v=[]
  for context,select in [('ACKR1_detected',a),('ACKR1_undetected',~a)]:
   sel=(dd<=rad)&select;den=int(sel.sum());num=int((sel&(label[ii]==typ)).sum());v.append((num,den));rows.append({'radius_um':rad,'identity_candidate':typ,'nearest_EC_context':context,'target_n':num,'all_nonEC_n':den,'fraction_nonEC':num/den if den else None,'context_EC_n':int((ack[ec]>=1).sum()) if context=='ACKR1_detected' else int((ack[ec]<1).sum())})
pd.DataFrame(rows).to_csv(o/'ACKR1_context_descriptive.tsv',sep='\t',index=False);print(pd.DataFrame(rows).query('radius_um==50').to_json(orient='records'))
