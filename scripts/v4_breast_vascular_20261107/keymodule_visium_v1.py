"""Composition-level extension, not a replacement of original CoVarNet CM scores."""
from pathlib import Path
import sys,json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy import sparse
from scipy.stats import rankdata
r=Path(sys.argv[1]);out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=False)
plan={'groups':{'B':['Bcell_MS4A1_CD79A','Bcell_Naive_IGHD','Bcell_Naive_TCL1A'],'EC':['Endothelial_Blood_PECAM1','Endothelial_VWF_PECAM1','Endothelial_Venous_ACKR1_ACKR1'],'NK':['NK_FCGR3A_cytotoxic_FCGR3A','NK_KLRD1_GNLY','NK_XCL1_CD56bright_XCL1']},'primary_edges':[['B','EC'],['EC','NK']],'exploratory_edge':['B','NK'],'posterior':['means','q05'],'aggregation':'sum disjoint fitted subtype abundances; q05 sum is sensitivity, NOT broad posterior 5th quantile','statistics':'Spearman and directed bivariate Moran observed statistic with previously audited author 60NN inverse squared row-standardized weights','no_pvalues':True,'scope':'all currently independently audited sections, no effect-based selection; incomplete cohort descriptive only','no_exact_subtype_claim':True}
(out/'configuration_before_effects.json').write_text(json.dumps(plan,indent=2)+'\n')
rows=[];receipts=[]
for receipt in sorted((r/'sections').glob('S*/independent_audit.json')):
 d=receipt.parent;z=json.loads(receipt.read_text());assert z['status'].startswith('PASS_');coords=pd.read_csv(d/'coordinates_PRIVATE.tsv.gz',sep='\t',index_col=0);w=pd.read_csv(d/'author_weights60_PRIVATE.tsv',sep='\t');ii=w.i.to_numpy()-1;jj=w.j.to_numpy()-1
 xy=coords.to_numpy(float);inv=1/((xy[ii]-xy[jj])**2).sum(axis=1);expected=inv/np.bincount(ii,weights=inv,minlength=len(coords))[ii];assert np.allclose(expected,w.weight,atol=1e-10)
 G=sparse.csr_matrix((w.weight,(ii,jj)),shape=(len(coords),len(coords)));assert np.allclose(np.asarray(G.sum(axis=1)).ravel(),1)
 for posterior in plan['posterior']:
  a=pd.read_csv(d/(posterior+'_abundance_PRIVATE.tsv.gz'),sep='\t',index_col=0);assert a.index.equals(coords.index);v=np.column_stack([a[nodes].sum(axis=1).to_numpy() for nodes in plan['groups'].values()]);assert np.isfinite(v).all() and (v>=0).all()
  for i,j in [(0,1),(1,2),(0,2)]:
   x,y=v[:,i],v[:,j];valid=np.std(x)>0 and np.std(y)>0;rho=mf=mr=np.nan
   if valid:
    rho=float(np.corrcoef(rankdata(x),rankdata(y))[0,1]);zx=(x-x.mean())/x.std(ddof=1);zy=(y-y.mean())/y.std(ddof=1);mf=float(zx@(G@zy)/(len(x)-1));mr=float(zy@(G@zx)/(len(x)-1))
    # Independent Pandas rank correlation and explicit neighbour-edge sum.
    assert abs(rho-pd.Series(x).corr(pd.Series(y),method='spearman'))<1e-10
    assert abs(mf-np.sum(zx[ii]*w.weight.to_numpy()*zy[jj])/(len(x)-1))<1e-10
   rows.append({'section_PRIVATE':d.name,'dataset':z['dataset'],'posterior':posterior,'edge':['B','EC','NK'][i]+'--'+['B','EC','NK'][j],'role':'exploratory' if (i,j)==(0,2) else 'original_edge_broadened','spearman':rho,'Moran_forward':mf,'Moran_reverse':mr})
 receipts.append({'audit_sha256':hashlib.sha256(receipt.read_bytes()).hexdigest()})
t=pd.DataFrame(rows);t.to_csv(out/'section_results_PRIVATE.tsv',sep='\t',index=False)
summary=t.groupby(['dataset','posterior','edge','role']).agg(n_sections=('spearman','size'),n_estimable=('spearman','count'),median_spearman=('spearman','median'),min_spearman=('spearman','min'),max_spearman=('spearman','max'),median_Moran_forward=('Moran_forward','median'),median_Moran_reverse=('Moran_reverse','median')).reset_index();summary.to_csv(out/'source_aggregate.tsv',sep='\t',index=False)
result={'UTC':datetime.now(timezone.utc).isoformat(),'n_audited_sections':len(receipts),'requested_sections':53,'status':'DESCRIPTIVE_COMPOSITION_EXTENSION_NUMERIC_CHECKS_PASS_NOT_PATIENT_LEVEL_REPLICATION','input_audit_sha256':receipts,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(out/'complete.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'receipt':result,'aggregate':summary.to_dict('records')}))
