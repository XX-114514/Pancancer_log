"""All prespecified K/repeats for one completed original major, no state release."""
from pathlib import Path
from datetime import datetime,timezone
from itertools import combinations
import os,argparse,json,hashlib,gc
import numpy as np,pandas as pd,anndata as ad
from scipy.optimize import linear_sum_assignment
from scipy.stats import rankdata
from cnmf import cNMF
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def frame(p):
 with np.load(p,allow_pickle=True) as z:return pd.DataFrame(z['data'],index=z['index'],columns=z['columns'])
def unit(x):
 x=np.asarray(x,float);d=np.linalg.norm(x,axis=1,keepdims=True)
 return np.divide(x,d,out=np.full(x.shape,np.nan),where=d>0)
def ranks(x):
 y=rankdata(x,axis=0);y-=y.mean(axis=0);d=np.linalg.norm(y,axis=0)
 return np.divide(y,d,out=np.full(y.shape,np.nan),where=d>0)
def residual(x,g):
 y=np.asarray(x,float).copy()
 for v in np.unique(g):
  m=g==v;y[m]-=y[m].mean(axis=0)
 return y
def r2group(x,g):
 den=((x-x.mean(axis=0))**2).sum(axis=0);e=(residual(x,g)**2).sum(axis=0)
 return np.divide(den-e,den,out=np.full(den.shape,np.nan),where=den>0)
def top(weights,genes):return [set(genes[np.argsort(-x,kind='stable')[:50]][x[np.argsort(-x,kind='stable')[:50]]>0]) for x in weights]
def jac(x,y):return len(x&y)/len(x|y) if x|y else np.nan
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--major',required=True);ap.add_argument('--plan',required=True);args=ap.parse_args();major=args.major;P=Path(args.plan);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);C=V/'all_lineage_cNMF_state_contrast_20261014_v2'
 progress=json.loads((C/'progress.json').read_text());completed=[r for r in progress['summary'] if r['major']==major];assert len(completed)==1 and completed[0]['all_consensus_completed']
 root=C/completed[0]['private_relative_root'];before={str(p.relative_to(root)):sha(p) for p in root.rglob('*') if p.is_file()};assert len(list(root.rglob('*spectra.k_*.iter_*.df.npz')))==len(plan['K'])*plan['repeats']
 A=R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json';assert sha(A)==plan['accepted_input_inventory_sha256'];iv=json.loads(A.read_text());item=next(i for i in iv['files'] if i['major']==major);assert sha(item['path'])==item['sha256'];raw=ad.read_h5ad(item['path']);meta=raw.obs;source=meta.gse_id.to_numpy().astype(str);gn=(meta.gse_id.astype(str)+'::'+meta.node_id.astype(str)).to_numpy();depth=np.column_stack([np.log1p(meta.full_RNA_total),np.log1p(meta.full_detected_genes)])
 model=cNMF(output_dir=str(root),name='cnmf');stats=[];repeats=[];mods=[];within=[];models={}
 for k in plan['K']:
  base=root/'cnmf/cnmf_tmp';pre='cnmf.';spec=frame(base/(pre+f'spectra.k_{k}.dt_0_5.consensus.df.npz'));use=frame(base/(pre+f'usages.k_{k}.dt_0_5.consensus.df.npz'));coef=frame(base/(pre+f'gene_spectra_score.k_{k}.dt_0_5.df.npz'));tpm=frame(base/(pre+f'gene_spectra_tpm.k_{k}.dt_0_5.df.npz'));merged=frame(base/(pre+f'spectra.k_{k}.merged.df.npz'));density=frame(base/(pre+f'local_density_cache.k_{k}.merged.df.npz'))
  assert use.shape==(raw.n_obs,k) and coef.shape==(k,raw.n_vars) and tpm.shape==coef.shape
  assert use.index.equals(meta.index) and coef.index.equals(use.columns) and spec.index.equals(use.columns) and tpm.index.equals(coef.index);assert list(coef.columns)==list(raw.var_names)
  assert np.isfinite(coef.values).all() and np.isfinite(use.values).all() and (use.values>=0).all();assert merged.shape==(plan['repeats']*k,spec.shape[1]) and merged.columns.equals(spec.columns) and density.index.equals(merged.index)
  ks=model.consensus(k=k,density_threshold=.5,local_neighborhood_size=.3,show_clustering=False,skip_density_and_return_after_stats=True)
  stats.append(dict(major=major,K=k,repeats=plan['repeats'],HVGs=spec.shape[1],unfiltered_silhouette=float(ks.loc['silhouette','stats']),unfiltered_prediction_error=float(ks.loc['prediction_error','stats']),density_retained_fraction=float(density.iloc[:,0].lt(.5).mean())))
  assert (use.values.sum(axis=1)>0).all();usage=use.values/use.values.sum(axis=1,keepdims=True);positive=np.maximum(coef.values,0);genes=coef.columns.to_numpy(str);tops=top(positive,genes)
  for repeat in range(plan['repeats']):
   keys=[f'iter{repeat}_topic{j+1}' for j in range(k)];assert all(v in merged.index for v in keys);cos=unit(merged.loc[keys].values)@unit(spec.values).T
   # Undefined zero profiles remain NA, while a penalized assignment permits a full audit table.
   ri,ci=linear_sum_assignment(np.where(np.isfinite(cos),-cos,1000))
   for i,j in zip(ri,ci):repeats.append(dict(major=major,K=k,repeat=repeat,original_topic=int(i)+1,consensus_component=int(spec.index[j]),weight_cosine=float(cos[i,j]),cosine_defined=bool(np.isfinite(cos[i,j])),density_retained=bool(density.loc[keys[i]].iloc[0]<.5)))
  rr=pd.DataFrame(repeats);rr=rr[rr.K.eq(k)];sr=r2group(usage,source);nr=r2group(usage,gn);partial=ranks(residual(usage,gn)).T@ranks(residual(depth,gn));normalized_depth=ranks(usage).T@ranks(depth)
  for j,component in enumerate(use.columns):
   q=rr[rr.consensus_component.eq(component)];assert len(q)==plan['repeats'];finite=float(q.cosine_defined.mean());med=float(q.weight_cosine.median());p10=float(q.weight_cosine.quantile(.1));ret=float(q.density_retained.mean());stable=finite>=plan['min_finite_repeat_fraction'] and med>=.8 and p10>=.7 and ret>=.5;signature=list(genes[np.argsort(-positive[j],kind='stable')[:50]]);signature=[g for g in signature if positive[j,np.flatnonzero(genes==g)[0]]>0];cols=raw.var_names.get_indexer(signature);det=np.asarray((raw.X[:,cols]>0).sum(axis=1)).ravel() if len(cols) else np.zeros(raw.n_obs,int);fraction=np.asarray(raw.X[:,cols].sum(axis=1)).ravel()/meta.full_RNA_total.to_numpy() if len(cols) else np.zeros(raw.n_obs)
   mods.append(dict(candidate_id='ACN_'+major.replace(' ','_').replace('/','_')+f'_K{k}_{int(component):02}',major=major,K=k,component=int(component),finite_repeat_fraction=finite,repeat_median_cosine=med,repeat_p10_cosine=p10,repeat_density_retained_fraction=ret,stable_repetition_rule=bool(stable),source_r2=float(sr[j]),source_node_r2=float(nr[j]),RNA_partial_rank_after_source_node=float(partial[j,0]),detected_partial_rank_after_source_node=float(partial[j,1]),normalized_usage_RNA_rho=float(normalized_depth[j,0]),primary50_positive_coefficient_genes=';'.join(signature),signature_size=len(signature),informative_gene_support_cells=int((det>=5).sum()),informative_gene_support_fraction=float((det>=5).mean()),signature_RNA_fraction_median=float(np.median(fraction)),source_limited=meta.gse_id.nunique()<2,formal_state=False))
   for g in sorted(set(source)):
    mask=source==g;u=usage[mask,j];within.append(dict(major=major,K=k,component=int(component),gse_id=g,cells=int(mask.sum()),normalized_usage_median=float(np.median(u)),normalized_usage_q10=float(np.quantile(u,.1)),normalized_usage_q90=float(np.quantile(u,.9)),informative_gene_support_cells=int((det[mask]>=5).sum()),signature_RNA_fraction_median=float(np.median(fraction[mask]))))
  models[k]=dict(coef=positive,genes=genes,normalized=usage,components=use.columns.to_numpy(),top=tops,TPM_top=top(tpm.values,genes));gc.collect()
 cross=[]
 for ka,kb in combinations(plan['K'],2):
  a=models[ka];b=models[kb];assert np.array_equal(a['genes'],b['genes']);cos=unit(a['coef'])@unit(b['coef']).T;rho=ranks(residual(a['normalized'],source)).T@ranks(residual(b['normalized'],source));node_rho=ranks(residual(a['normalized'],gn)).T@ranks(residual(b['normalized'],gn))
  for i,u in enumerate(a['components']):
   for j,v in enumerate(b['components']):cross.append(dict(major=major,K_a=ka,component_a=int(u),K_b=kb,component_b=int(v),positive_coef_cosine=float(cos[i,j]),top50_jaccard=jac(a['top'][i],b['top'][j]),source_residual_usage_rho=float(rho[i,j]),source_node_residual_usage_rho=float(node_rho[i,j]),comparison_defined=bool(np.isfinite(cos[i,j]) and np.isfinite(rho[i,j]))))
 modules=pd.DataFrame(mods);cx=pd.DataFrame(cross);stablemap={(r.K,r.component):r.stable_repetition_rule for r in modules.itertuples()}
 for i,row in modules.iterrows():
  q=cx[((cx.K_a==row.K)&(cx.component_a==row.component))|((cx.K_b==row.K)&(cx.component_b==row.component))];supported=q[(q.positive_coef_cosine>=.8)&(q.top50_jaccard>=.25)&(q.source_residual_usage_rho>=.5)];count=0
  for r in supported.itertuples():count+=bool(stablemap[(r.K_b,r.component_b) if r.K_a==row.K else (r.K_a,r.component_a)])
  modules.loc[i,'crossK_stable_support_pairs']=count;modules.loc[i,'repeat_and_crossK_candidate']=bool(row.stable_repetition_rule and count>0)
 O=V/'all_lineage_cNMF_diagnostics_20261014_v2'/completed[0]['private_relative_root'];assert not O.exists();O.mkdir(parents=True)
 for name,x in [('all45_candidates.tsv',modules),('all4500_repetition_matches.tsv',pd.DataFrame(repeats)),('all650_crossK_comparisons.tsv',cx),('all3_K_diagnostics.tsv',pd.DataFrame(stats)),('all_source_distributions.tsv',pd.DataFrame(within))]:x.to_csv(O/name,sep='\t',index=False)
 after={str(p.relative_to(root)):sha(p) for p in root.rglob('*') if p.is_file()};assert before==after
 out=dict(status='COMPLETE_ORIGINAL_MAJOR_ALL_K_CNMFDIAGNOSTICS_NOT_STATES',major=major,observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),original_major_completed_record=completed[0],all_fit_outputs_before_after_unchanged=True,input_output_sha256=before,candidates=len(modules),repetition_matches=len(repeats),crossK_comparisons=len(cx),repeat_crossK_candidates=int(modules.repeat_and_crossK_candidate.sum()),clinical_LR_p_used=False,formal_states=False,output_sha256={p.name:sha(p) for p in O.iterdir() if p.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='input_output_sha256'}),flush=True)
if __name__=='__main__':main()
