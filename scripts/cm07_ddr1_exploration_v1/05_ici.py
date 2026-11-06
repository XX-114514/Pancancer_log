import sys; from pathlib import Path; import numpy as np, pandas as pd
D=Path(sys.argv[1]); rng=np.random.default_rng(20261106); I=D/'imvigor210'
e=pd.read_csv(I/'data_mrna_seq_tpm.txt',sep='\t',low_memory=False); print(e.iloc[:2,:4].to_string()); s=e.columns[0]
e=e.dropna(subset=[s]).drop_duplicates(s).set_index(s); e=e.select_dtypes('number'); print(e.shape)
if 'DDR1' not in e.index: print('no DDR1'); sys.exit()
x=np.log2(e+1); x=x.T; z=(x-x.mean())/x.std()
COL=['COL1A1','COL1A2','COL3A1']; cyt=z[['GZMA','PRF1']].mean(1); L=z[COL].mean(1); R=z['DDR1']
pu=z[[g for g in ['EPCAM','KRT19','KRT8'] if g in z]].mean(1)
def ols(y,X):
    X=np.column_stack([np.ones(len(y))]+X); b,*_=np.linalg.lstsq(X,y,rcond=None); r=y-X@b; V=np.linalg.inv(X.T@X)*((r@r)/(len(y)-X.shape[1])); return b[1:],(b/np.sqrt(np.diag(V)))[1:]
b,t=ols(cyt.values,[L.values,R.values,(L*R).values]); print('IMvigor CYT: beta L,R,LxR',b.round(3),'t',t.round(2),'n',len(cyt),'r DDR1-epithelial',round(np.corrcoef(R,pu)[0,1],3))
m=x.mean(); sd=x.std(); ok=(sd>0)&~x.columns.isin(COL+['DDR1']); fv=np.column_stack([(m[ok]-m[ok].mean())/m[ok].std(),(sd[ok]-sd[ok].mean())/sd[ok].std()]); tg=np.array([(m['DDR1']-m[ok].mean())/m[ok].std(),(sd['DDR1']-sd[ok].mean())/sd[ok].std()])
pool=x.columns[ok][np.argsort(((fv-tg)**2).sum(1))[:500]]; pick=rng.choice(pool,300,replace=False)
nu=np.array([ols(cyt.values,[L.values,z[g].values,(L*z[g]).values])[1] for g in pick]); print('null t R median',np.median(nu[:,1]).round(2),'q05',np.quantile(nu[:,1],.05).round(2),'p_R_neg',(1+(nu[:,1]<=t[1]).sum())/301,'; LxR median',np.median(nu[:,2]).round(2),'q95',np.quantile(nu[:,2],.95).round(2),'p_LxR_pos',(1+(nu[:,2]>=t[2]).sum())/301)
