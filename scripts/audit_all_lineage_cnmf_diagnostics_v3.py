"""Independent cdist/Hungarian, dummy regressions, scalar Spearman and incidence."""
from pathlib import Path
from datetime import datetime,timezone
from itertools import combinations
import os,json,hashlib
import numpy as np,pandas as pd,anndata as ad
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr
from scipy.sparse import csr_matrix
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);major=os.environ['PANCANCER_STATE_AUDIT_MAJOR'];V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);C=V/'all_lineage_cNMF_state_contrast_20261014_v2';record=next(r for r in json.loads((C/'progress.json').read_text())['summary'] if r['major']==major);root=C/record['private_relative_root'];D=V/'all_lineage_cNMF_diagnostics_20261014_v2'/record['private_relative_root'];iv=json.loads((D/'inventory.json').read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for name,h in iv['output_sha256'].items():assert sha(D/name)==h
for name,h in iv['input_output_sha256'].items():assert sha(root/name)==h
def frame(p):
 with np.load(p,allow_pickle=True) as z:return pd.DataFrame(z['data'],index=z['index'],columns=z['columns'])
def rho(x,y):
 return float(spearmanr(x,y).statistic) if np.ptp(x)>1e-12 and np.ptp(y)>1e-12 else np.nan
def eq(x,y):return (np.isnan(x) and np.isnan(y)) or np.isclose(x,y,atol=1e-7)
def dummy(v):return pd.get_dummies(pd.Series(v,dtype=str),dtype=float).to_numpy()
def resid(x,d):
 out=x-d@np.linalg.lstsq(d,x,rcond=None)[0]
 # Mathematically constant one-hot groups have exact zero residuals.
 for j in range(d.shape[1]):
  ii=np.flatnonzero(d[:,j]==1)
  if len(ii):
   if x.ndim==1:
    if np.ptp(x[ii])==0:out[ii]=0
   else:
    jj=np.flatnonzero(np.ptp(x[ii],axis=0)==0)
    if len(jj):out[np.ix_(ii,jj)]=0
 return out
A=R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json';item=next(r for r in json.loads(A.read_text())['files'] if r['major']==major);assert sha(item['path'])==item['sha256'];raw=ad.read_h5ad(item['path']);meta=raw.obs;ds=dummy(meta.gse_id.to_numpy());dn=dummy((meta.gse_id.astype(str)+'::'+meta.node_id.astype(str)).to_numpy());depth=np.column_stack([np.log1p(meta.full_RNA_total),np.log1p(meta.full_detected_genes)]);dr=resid(depth,dn)
mods=pd.read_csv(D/'all45_candidates.tsv',sep='\t');repeats=pd.read_csv(D/'all4500_repetition_matches.tsv',sep='\t');cross=pd.read_csv(D/'all650_crossK_comparisons.tsv',sep='\t');models={};maxdiff=0.;nrep=0;support_checks=0
for k in [10,15,20]:
 base=root/'cnmf/cnmf_tmp';spec=frame(base/f'cnmf.spectra.k_{k}.dt_0_5.consensus.df.npz');merged=frame(base/f'cnmf.spectra.k_{k}.merged.df.npz');density=frame(base/f'cnmf.local_density_cache.k_{k}.merged.df.npz');use=frame(base/f'cnmf.usages.k_{k}.dt_0_5.consensus.df.npz');coef=frame(base/f'cnmf.gene_spectra_score.k_{k}.dt_0_5.df.npz');positive=np.clip(coef.values,0,None);usage=use.values/use.values.sum(axis=1)[:,None];genes=coef.columns.to_numpy(str)
 for repeat in range(100):
  keys=[f'iter{repeat}_topic{j+1}' for j in range(k)];cos=1-cdist(merged.loc[keys].values,spec.values,metric='cosine');ri,ci=linear_sum_assignment(np.where(np.isfinite(cos),-cos,1000))
  for i,j in zip(ri,ci):
   q=repeats[(repeats.K==k)&(repeats['repeat']==repeat)&(repeats.original_topic==i+1)];assert len(q)==1;r=q.iloc[0];assert r.consensus_component==spec.index[j] and r.cosine_defined==bool(np.isfinite(cos[i,j])) and eq(cos[i,j],r.weight_cosine) and r.density_retained==bool(density.loc[keys[i]].iloc[0]<.5)
   if np.isfinite(cos[i,j]):maxdiff=max(maxdiff,abs(cos[i,j]-r.weight_cosine))
   nrep+=1
 source_res=resid(usage,ds);node_res=resid(usage,dn);den=((usage-usage.mean(axis=0))**2).sum(axis=0);sr=np.divide(den-(source_res**2).sum(axis=0),den,out=np.full(den.shape,np.nan),where=den>0);nr=np.divide(den-(node_res**2).sum(axis=0),den,out=np.full(den.shape,np.nan),where=den>0);tops=[];rows=[];cols=[]
 for j,component in enumerate(use.columns):
  rr=repeats[(repeats.K==k)&(repeats.consensus_component==component)];r=mods[(mods.K==k)&(mods.component==component)].iloc[0];finite=float(rr.cosine_defined.mean());med=float(np.nanmedian(rr.weight_cosine));p10=float(np.nanquantile(rr.weight_cosine,.1));ret=float(rr.density_retained.mean());assert eq(med,r.repeat_median_cosine) and eq(p10,r.repeat_p10_cosine) and eq(ret,r.repeat_density_retained_fraction) and eq(finite,r.finite_repeat_fraction);assert r.stable_repetition_rule==bool(finite>=.9 and med>=.8 and p10>=.7 and ret>=.5)
  assert eq(sr[j],r.source_r2) and eq(nr[j],r.source_node_r2) and eq(rho(node_res[:,j],dr[:,0]),r.RNA_partial_rank_after_source_node) and eq(rho(node_res[:,j],dr[:,1]),r.detected_partial_rank_after_source_node) and eq(rho(usage[:,j],depth[:,0]),r.normalized_usage_RNA_rho)
  order=np.argsort(-positive[j],kind='stable')[:50];order=order[positive[j,order]>0];names=genes[order];tops.append(set(names));assert list(names)==str(r.primary50_positive_coefficient_genes).split(';') and len(names)==r.signature_size;idx=raw.var_names.get_indexer(names);assert (idx>=0).all();rows.extend(idx.tolist());cols.extend([j]*len(idx))
 incidence=csr_matrix((np.ones(len(rows),np.int32),(rows,cols)),shape=(raw.n_vars,k));det=((raw.X>0).astype(np.int32)@incidence).toarray();frac=(raw.X@incidence).toarray()/meta.full_RNA_total.to_numpy()[:,None]
 for j,component in enumerate(use.columns):
  r=mods[(mods.K==k)&(mods.component==component)].iloc[0];n=int((det[:,j]>=5).sum());assert r.informative_gene_support_cells==n and eq(r.informative_gene_support_fraction,n/raw.n_obs) and eq(r.signature_RNA_fraction_median,float(np.median(frac[:,j])));support_checks+=1
 models[k]=dict(positive=positive,components=list(use.columns),tops=tops,source_res=source_res,node_res=node_res)
nxc=0
for ka,kb in combinations([10,15,20],2):
 a=models[ka];b=models[kb];cos=1-cdist(a['positive'],b['positive'],metric='cosine')
 for i,u in enumerate(a['components']):
  for j,v in enumerate(b['components']):
   r=cross[(cross.K_a==ka)&(cross.K_b==kb)&(cross.component_a==u)&(cross.component_b==v)].iloc[0];union=a['tops'][i]|b['tops'][j];jac=len(a['tops'][i]&b['tops'][j])/len(union) if union else np.nan;value=rho(a['source_res'][:,i],b['source_res'][:,j]);nv=rho(a['node_res'][:,i],b['node_res'][:,j]);assert eq(cos[i,j],r.positive_coef_cosine) and eq(jac,r.top50_jaccard) and eq(value,r.source_residual_usage_rho) and eq(nv,r.source_node_residual_usage_rho);nxc+=1
stable={(r.K,r.component):r.stable_repetition_rule for r in mods.itertuples()}
for r in mods.itertuples():
 count=0
 for x in cross.itertuples():
  if (x.K_a,x.component_a)==(r.K,r.component):other=(x.K_b,x.component_b)
  elif (x.K_b,x.component_b)==(r.K,r.component):other=(x.K_a,x.component_a)
  else:continue
  if x.positive_coef_cosine>=.8 and x.top50_jaccard>=.25 and x.source_residual_usage_rho>=.5 and stable[other]:count+=1
 assert count==r.crossK_stable_support_pairs and r.repeat_and_crossK_candidate==bool(r.stable_repetition_rule and count>0)
assert nrep==4500 and nxc==650 and support_checks==45
receipt=dict(status='PASS_ALL_MAJOR4500_REPEAT650_CROSSK45_SOURCE_NODE_DEPTH_NATIVE_SUPPORT_AND_RULES',major=major,observed_UTC=datetime.now(timezone.utc).isoformat(),audit_code_sha256=sha(__file__),producer_inventory_sha256=sha(D/'inventory.json'),repetition_matches=nrep,crossK_comparisons=nxc,module_native_support_and_regression_checks=support_checks,max_repeat_cosine_difference=maxdiff,K_statistics_not_independently_refit=True,source_distribution_quantiles_not_independently_recomputed_here=True,formal_states=False,mathematically_constant_group_residuals_exactzero=True,original_v1_failed_audit_retained=True)
O=V/'all_lineage_cNMF_diagnostic_audit_20261014_v3'/record['private_relative_root'];assert not O.exists();O.mkdir(parents=True);(O/'independent_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
