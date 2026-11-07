import json,sys,subprocess
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,scipy.sparse as sp
from scipy.stats import rankdata
r=PROJECT_ROOT/'runs/20261107_v4_Visium_author_validation_v2';base=PROJECT_ROOT/'runs/20261106_v7_selective_CoVarNet_paper_v2';rb=Path(sys.executable).parents[3]/'envs/scrna_r/bin/Rscript';out=r/'audit/math_fixture_PRIVATE';p=subprocess.run([str(rb),'--vanilla',str(r/'scripts/spatial_math_fixture.R'),str(r),str(base)],capture_output=True,text=True,timeout=300);(r/'audit/math_fixture_R_PRIVATE.log').open('x').write(p.stdout+p.stderr);assert p.returncode==0,p.stderr[-3000:];a=pd.read_csv(out/'abundance.tsv',sep='\t',index_col=0);w=pd.read_csv(out/'weights.tsv',sep='\t');G=sp.csr_matrix((w.weight,(w.i.to_numpy()-1,w.j.to_numpy()-1)),shape=(100,100));X=a.to_numpy();Z=(X-X.mean(axis=0))/X.std(axis=0,ddof=1);M=Z.T@(G@Z)/99;C=np.corrcoef(rankdata(X,axis=0).T);mem=pd.read_csv(base/'CM_results/each_node.tsv',sep='\t');errors=[]
for name,mat,directed in [('coloc',C,False),('moran',M,True)]:
 table=pd.read_csv(out/(name+'.tsv'),sep='\t');ii,jj=np.where(~np.eye(len(a.columns),dtype=bool)) if directed else np.triu_indices(len(a.columns),1)
 for z in table.itertuples():
  take=a.columns.isin(mem.loc[mem.cm.eq(z.cm),'subCluster']);inside=take[ii]&take[jj];bg=take[ii]|take[jj];score=np.mean(mat[ii[bg],jj[bg]]<=np.median(mat[ii[inside],jj[inside]]));errors.append(abs(score-z.score));assert abs(score-z.score)<1e-12
proof={'status':'PASS_literal_author_vs_independent_Spearman_and_directed_60NN_bivariate_Moran_math','UTC':datetime.now(timezone.utc).isoformat(),'fixture_only_not_real_spatial_validation':True,'n_synthetic_spots':100,'n_real_label_names':len(a.columns),'n_nonempty_CM_scores_checked':len(errors),'maximum_score_error':max(errors),'neighbors':60,'weights':'inverse square row normalized directed','moran_definition':'Z_transpose W Z /(n-1); n-1 sample standardization exact installed spdep T0'};(r/'audit/spatial_math_fixture_audit.json').open('x').write(json.dumps(proof,indent=2));print(json.dumps(proof))
