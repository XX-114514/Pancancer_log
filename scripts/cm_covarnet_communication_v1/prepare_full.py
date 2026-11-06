"""One sequential native-count pass for full CMT LR input and sufficient statistics."""
from pathlib import Path
import numpy as np,pandas as pd,scipy.sparse as sp,zarr,anndata as ad,json,zipfile,hashlib,time
from anndata.io import read_elem
r=Path(__file__).resolve().parents[1];q=pd.read_pickle(r/'private/cell_metadata.pkl.gz');q['identity']=q.final_annotation_with_malignancy_v6.astype(str)
z=zarr.open(str(r.parent/'20260818_pancancer_5M_merge_prep/02_merge_new_cohort/objects/5M_merge_concat.zarr'),mode='r');m=z['X'];genes=read_elem(z['var']).index.astype(str);ip=m['indptr'][:]
with zipfile.ZipFile(r/'sources/cellphonedb_v4.1.0.zip') as db:gf=pd.read_csv(db.open('gene_table.csv'))
keep=genes.isin(gf.hgnc_symbol.dropna());pd.DataFrame({'gene':gf.hgnc_symbol.dropna().unique()}).assign(measured=lambda d:d.gene.isin(genes)).to_csv(r/'audit/LR_gene_measurability.tsv',sep='\t',index=False)
valid=q.CMT.notna().to_numpy();obs=q.loc[valid,['identity','sampleID','CMT']].copy();obs.index=pd.Index(['c'+str(i) for i in np.where(valid)[0]])
ckeys=pd.MultiIndex.from_frame(q[['CMT','identity']].fillna('NO_H'));cids,cu=pd.factorize(ckeys,sort=True);nc=len(cu);ng=len(genes)
sums=np.zeros((nc,ng));sq=np.zeros((nc,ng));det=np.zeros((nc,ng),dtype='int64');num=np.bincount(cids,minlength=nc)
blocks=[];rawcheck=[];t=time.time()
for lo in range(0,len(q),20000):
 hi=min(lo+20000,len(q));a,b=int(ip[lo]),int(ip[hi]);d=np.asarray(m['data'][a:b]);ii=np.asarray(m['indices'][a:b],dtype='int32');assert np.isfinite(d).all() and (d>=0).all() and np.allclose(d,np.round(d),rtol=0,atol=1e-5)
 x=sp.csr_matrix((d,ii,ip[lo:hi+1]-a),shape=(hi-lo,ng));total=np.asarray(x.sum(axis=1)).ravel();fac=np.divide(1e4,total,out=np.zeros_like(total),where=total>0);x=x.astype('float32');x.data*=np.repeat(fac,np.diff(x.indptr));np.log1p(x.data,out=x.data)
 agg=sp.csr_matrix((np.ones(hi-lo),(cids[lo:hi],np.arange(hi-lo))),shape=(nc,hi-lo));sums+=(agg@x).toarray();xx=x.copy();xx.data**=2;sq+=(agg@xx).toarray();xx.data[:]=1;det+=(agg@xx).toarray().astype('int64')
 blocks.append(x[valid[lo:hi]][:,keep]);rawcheck.append({'lo':lo,'hi':hi,'sha256':hashlib.sha256(d.tobytes()).hexdigest()})
 if lo%200000==0:print('ROWS',hi,'seconds',round(time.time()-t),flush=True)
X=sp.vstack(blocks,format='csr');adata=ad.AnnData(X,obs=obs,var=pd.DataFrame(index=genes[keep]));adata.uns['log1p']={'base':None};adata.write_h5ad(r/'inputs/full_LR.h5ad',compression='gzip');np.savez_compressed(r/'private/CMT_identity_stats.npz',sums=sums,sq=sq,det=det,n=num,genes=genes.to_numpy(dtype=str));pd.DataFrame(cu.tolist(),columns=['CMT','identity']).to_csv(r/'private/CMT_identity_stats_index.tsv',sep='\t',index=False);pd.DataFrame(rawcheck).to_csv(r/'private/full_count_block_SHA.tsv',sep='\t',index=False)
res={'cells':X.shape[0],'LR_genes':X.shape[1],'original_genes':ng,'seconds':time.time()-t,'raw_counts_all_blocks_integer_finite_nonnegative':True,'full_lr_sha256':hashlib.sha256((r/'inputs/full_LR.h5ad').read_bytes()).hexdigest()};(r/'audit/full_input_complete.json').write_text(json.dumps(res,indent=2));print(res,flush=True)
