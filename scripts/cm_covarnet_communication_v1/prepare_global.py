"""Read native CSR in bounded blocks; fixed per-identity sample, all original samples; H eligibility is only required for comparisons."""
from pathlib import Path
import pandas as pd,numpy as np,scipy.sparse as sp,zarr,anndata as ad,json,hashlib,time
from anndata.io import read_elem
r=Path(__file__).resolve().parents[1]; p=r.parent/'20260818_pancancer_5M_merge_prep/02_merge_new_cohort/objects/5M_merge_concat.zarr'
q=pd.read_pickle(r/'private/cell_metadata.pkl.gz'); q=q.copy();q['identity']=q.final_annotation_with_malignancy_v6.astype(str)
rng=np.random.default_rng(1337);selected=[]
for ident,g in q.groupby('identity',sort=True):selected.extend(rng.choice(g.global_position.to_numpy(),min(1000,len(g)),replace=False))
selected=np.sort(selected);z=zarr.open(str(p),mode='r');m=z['X'];genes=read_elem(z['var']).index.astype(str);assert genes.is_unique
ip=np.asarray(m['indptr'][:]);blocks=[];checks=[];start=time.time()
for lo in range(0,int(m.attrs['shape'][0]),20000):
 hi=min(lo+20000,int(m.attrs['shape'][0]));sel=selected[(selected>=lo)&(selected<hi)]
 if not len(sel):continue
 a,b=int(ip[lo]),int(ip[hi]);data=np.asarray(m['data'][a:b]);ind=np.asarray(m['indices'][a:b],dtype='int32')
 assert np.isfinite(data).all() and (data>=0).all() and np.allclose(data,np.round(data),rtol=0,atol=1e-5),'non-count X'
 x=sp.csr_matrix((data,ind,ip[lo:hi+1]-a),shape=(hi-lo,len(genes)))[sel-lo].astype('float32');tot=np.asarray(x.sum(axis=1)).ravel();assert (tot>0).all()
 x.data*=np.repeat(1e4/tot,np.diff(x.indptr));np.log1p(x.data,out=x.data);blocks.append(x)
 checks.append({'start':lo,'end':hi,'nnz':len(data),'selected':len(sel),'data_sha256':hashlib.sha256(data.tobytes()).hexdigest()})
 if lo%200000==0:print('ROWS',hi,'elapsed',round(time.time()-start),flush=True)
X=sp.vstack(blocks,format='csr');obs=q.set_index('global_position').loc[selected,['identity','sampleID','CMT']].copy();obs.index=pd.Index(['c'+str(i) for i in selected]);a=ad.AnnData(X,obs=obs,var=pd.DataFrame(index=genes));a.uns['log1p']={'base':None};a.write_h5ad(r/'inputs/global_1000.h5ad',compression='gzip')
obs[['identity']].reset_index(names='Cell').rename(columns={'identity':'cell_type'}).to_csv(r/'inputs/global_meta.tsv',sep='\t',index=False)
pd.DataFrame(checks).to_csv(r/'private/global_counts_block_hashes.tsv',sep='\t',index=False)
manifest={'n_cells':len(obs),'n_genes':len(genes),'n_identities':obs.identity.nunique(),'seconds':time.time()-start,'counts_integer_verified_on_blocks':True,'normalization':'total_1e4_then_log1p','selection_sha256':hashlib.sha256(selected.tobytes()).hexdigest(),'h5ad_sha256':hashlib.sha256((r/'inputs/global_1000.h5ad').read_bytes()).hexdigest()}
(r/'audit/global_input_complete.json').write_text(json.dumps(manifest,indent=2));print(manifest,flush=True)
