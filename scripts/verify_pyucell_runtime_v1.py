"""Official CPU-average-tie UCell check; synthetic only, no biological scores."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,importlib.metadata,inspect
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse
from scipy.stats import rankdata
import pyucell as uc

def main():
    assert importlib.metadata.version('pyucell')=='0.8.0'
    raw=np.array([[5,5,2,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,1,2,3,4,5,6],[7,1,3,3,3,1,0,0]],dtype=float)
    genes=np.array(list('ABCDEFGH'));obs=pd.DataFrame(index=['fixture'+str(i) for i in range(4)]);sig={'AB':['A','B'],'ABC':['A','B','C'],'partial':['A','B','UNKNOWN']};maxrank=6
    def score(x,ids):
        a=ad.AnnData(X=sparse.csr_matrix(x),obs=pd.DataFrame(index=ids),var=pd.DataFrame(index=genes))
        uc.compute_ucell_scores(a,signatures=sig,max_rank=maxrank,ties_method='average',missing_genes='skip',chunk_size=2,n_jobs=1,device='cpu');return a
    a=score(raw,obs.index);expected=[]
    for cell in raw:
        ranks=np.full(8,maxrank,dtype=float);nz=np.where(cell>0)[0];r=rankdata(-cell[nz],method='average');ranks[nz]=np.where(r<=maxrank,r,maxrank)
        values=[]
        for ids in [[0,1],[0,1,2],[0,1]]:
            n=len(ids);u=sum(ranks[j] for j in ids)-n*(n+1)/2;values.append(1-u/(n*maxrank-n*(n+1)/2))
        expected.append(values)
    got=a.obs[['AB_UCell','ABC_UCell','partial_UCell']].to_numpy();assert np.allclose(got,expected,atol=1e-7);assert np.all((got>=0)&(got<=1));assert np.all(got[1]==0)
    b=score(np.vstack([raw,np.array([[9,2,0,1,1,0,0,3]])]),list(obs.index)+['extra']);assert np.array_equal(got,b.obs.iloc[:4][['AB_UCell','ABC_UCell','partial_UCell']].to_numpy())
    dense=ad.AnnData(X=raw.copy(),obs=obs.copy(),var=pd.DataFrame(index=genes));uc.compute_ucell_scores(dense,signatures=sig,max_rank=maxrank,ties_method='average',missing_genes='skip',chunk_size=3,n_jobs=1,device='cpu');assert np.array_equal(got,dense.obs[['AB_UCell','ABC_UCell','partial_UCell']].to_numpy())
    capped=ad.AnnData(X=sparse.csr_matrix(raw),obs=obs.copy(),var=pd.DataFrame(index=genes));uc.compute_ucell_scores(capped,signatures={'AB':['A','B']},max_rank=3,ties_method='average',missing_genes='skip',n_jobs=1,device='cpu');assert np.isfinite(capped.obs.AB_UCell).all()
    O=Path(os.environ['PANCANCER_PYUCELL_RUNTIME']);p=O/'synthetic_runtime_audit.json';assert not p.exists()
    out=dict(status='PASS_OFFICIAL_CPU_UCELL_SCALAR_TIES_ZERO_MISSING_BATCH_DENSE_SPARSE',observed_UTC=datetime.now(timezone.utc).isoformat(),version='0.8.0',verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),api_signature=str(inspect.signature(uc.compute_ucell_scores)),scalar_values_checked=12,max_scalar_difference=float(np.max(np.abs(got-np.array(expected)))),zero_measured_cell_score_zero=True,missing_gene_skip_confirmed=True,missing_is_not_measured_zero=True,batch_composition_invariant=True,dense_sparse_equal=True,primary_device='cpu',primary_ties_method='average',biological_scores_computed=False,frozen_state_dictionary=False,dependency_versions={x:importlib.metadata.version(x) for x in ['pyucell','anndata','numpy','scipy','scanpy','joblib','session-info2']})
    p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
