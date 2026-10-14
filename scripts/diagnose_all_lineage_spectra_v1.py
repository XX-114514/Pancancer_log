"""Every block/factor/seed, identity-within-node diagnostics; no clinical/LR tests."""
from pathlib import Path
from datetime import datetime,timezone
from itertools import combinations
import argparse,os,json,hashlib,gc
import numpy as np,pandas as pd,anndata as ad
from scipy.optimize import linear_sum_assignment
from scipy.stats import rankdata

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def corr(x,y):
 x=rankdata(x);y=rankdata(y);x=x-x.mean();y=y-y.mean();den=np.linalg.norm(x)*np.linalg.norm(y)
 return float(x@y/den) if den>0 else np.nan
def residual(x,group):
 out=np.asarray(x,float).copy()
 for g in np.unique(group):
  keep=group==g;out[keep]-=out[keep].mean(axis=0)
 return out
def r2group(x,group):
 x=np.asarray(x,float);den=((x-x.mean(axis=0))**2).sum(axis=0);err=(residual(x,group)**2).sum(axis=0)
 return np.divide(den-err,den,out=np.full(den.shape,np.nan),where=den>0)
def normalized(x):
 x=np.asarray(x,float);den=np.linalg.norm(x,axis=1,keepdims=True)
 return np.divide(x,den,out=np.full(x.shape,np.nan),where=den>0)
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--block',required=True);ap.add_argument('--plan',required=True);args=ap.parse_args()
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=Path(args.plan);plan=json.loads(P.read_text());assert sha(__file__)==plan['producer_sha256'];block=args.block;assert block in plan['blocks']
 S=R/'all_lineage_Spectra_state_discovery_20261014_v1';root=S/block;D=json.loads((root/'fixed_feature_dictionary.json').read_text());seeds=plan['seeds'];models=[];obs=[];inputsha={}
 for seed in seeds:
  p=root/('seed'+str(seed));path=p/'factors_and_discovery_scores.npz'
  # A seed is complete only after producer publishes its progress record.
  progress=json.loads((S/'progress.json').read_text());assert any(r['block']==block and r['seed']==seed for r in progress['summary'])
  with np.load(path,allow_pickle=False) as z:models.append({k:z[k].copy() for k in z.files})
  obs.append(pd.read_csv(p/'private_discovery_metadata.tsv.gz',sep='\t',index_col=0,keep_default_na=False));inputsha[str(path.relative_to(S))]=sha(path);inputsha[str((p/'private_discovery_metadata.tsv.gz').relative_to(S))]=sha(p/'private_discovery_metadata.tsv.gz')
 ref=models[0];meta=obs[0];scope=ref['factor_scope'];assert set(scope)==({'block_global'}|set(D['majors']))
 for model,md in zip(models,obs):
  for key in ['cell_ids','genes','factor_scope']:assert np.array_equal(model[key],ref[key])
  assert np.array_equal(md.index.to_numpy(str),ref['cell_ids']) and md.equals(meta)
  assert np.isfinite(model['factors']).all() and np.isfinite(model['cell_scores']).all() and (model['factors']>=0).all() and (model['cell_scores']>=0).all()
 maps={};pairs=[];scores=[];tops=[np.argsort(-m['factors'],axis=1,kind='stable')[:,:50] for m in models];units=[normalized(m['factors']) for m in models]
 for a,b in combinations(range(3),2):
  mapping={}
  for s in ['block_global',*D['majors']]:
   ids=np.flatnonzero(scope==s);cos=units[a][ids]@units[b][ids].T;assert np.isfinite(cos).all();ri,ci=linear_sum_assignment(-cos)
   for u,v in zip(ids[ri],ids[ci]):mapping[int(u)]=int(v)
  maps[a,b]=mapping
  for u,v in mapping.items():
   ga=set(tops[a][u]);gb=set(tops[b][v]);pairs.append(dict(block=block,seed_a=seeds[a],seed_b=seeds[b],factor_a=u,factor_b=v,scope=scope[u],cosine=float(units[a][u]@units[b][v]),top50_jaccard=len(ga&gb)/len(ga|gb)))
   for major in D['majors'] if scope[u]=='block_global' else [scope[u]]:
    keep=meta.major.eq(major).to_numpy();x=models[a]['cell_scores'][keep,u];y=models[b]['cell_scores'][keep,v];src=meta.gse_id.to_numpy()[keep];gn=(meta.gse_id.astype(str)+'::'+meta.node_id.astype(str)).to_numpy()[keep]
    scores.append(dict(block=block,seed_a=seeds[a],seed_b=seeds[b],factor_a=u,factor_b=v,major=major,cells=int(keep.sum()),pooled_score_rho=corr(x,y),source_residual_score_rho=corr(residual(x,src),residual(y,src)),source_node_residual_score_rho=corr(residual(x,gn),residual(y,gn))))
 matched=[np.arange(len(scope)),np.array([maps[0,1][j] for j in range(len(scope))]),np.array([maps[0,2][j] for j in range(len(scope))])]
 pooled_weights=np.median(np.stack([units[i][matched[i]] for i in range(3)]),axis=0);pooled_top=np.argsort(-pooled_weights,axis=1,kind='stable')[:,:50]
 A=R/'all_lineage_state_native_input_audit_20261014_v2/accepted_input_inventory.json';assert sha(A)==plan['accepted_input_inventory_sha256'];accepted=json.loads(A.read_text());lookup={item['major']:item for item in accepted['files']}
 diagnostics=[];within=[];support={}
 for major in D['majors']:
  item=lookup[major];assert sha(item['path'])==item['sha256'];raw=ad.read_h5ad(item['path']);keep=meta.major.eq(major).to_numpy();md=meta.loc[keep];assert np.array_equal(raw.obs_names,md.index)
  ids=np.flatnonzero((scope=='block_global')|(scope==major));y=ref['cell_scores'][keep][:,ids];src=md.gse_id.to_numpy();gn=(md.gse_id.astype(str)+'::'+md.node_id.astype(str)).to_numpy();yr=residual(y,gn)
  depth=np.column_stack([np.log1p(md.full_RNA_total),np.log1p(md.full_detected_genes)]);dr=residual(depth,gn);fit=np.linalg.lstsq(dr,yr,rcond=None)[0];den=((y-y.mean(axis=0))**2).sum(axis=0);extra=np.divide((yr**2).sum(axis=0)-((yr-dr@fit)**2).sum(axis=0),den,out=np.full(den.shape,np.nan),where=den>0);sr=r2group(y,src);nr=r2group(y,gn)
  for jj,j in enumerate(ids):
   names=ref['genes'][pooled_top[j]];cols=raw.var_names.get_indexer(names);assert (cols>=0).all();det=np.asarray((raw.X[:,cols]>0).sum(axis=1)).ravel();rna=np.asarray(raw.X[:,cols].sum(axis=1)).ravel()/md.full_RNA_total.to_numpy();n_support=int((det>=plan['informative_genes']).sum());support[j,major]=n_support
   diagnostics.append(dict(block=block,anchor_factor=int(j),scope=scope[j],major=major,cells=len(md),sources=md.gse_id.nunique(),original_nodes=md.node_id.nunique(),source_r2=float(sr[jj]),source_node_r2=float(nr[jj]),node_increment_after_source=float(nr[jj]-sr[jj]),RNA_detected_increment_after_source_node=float(extra[jj]),RNA_partial_rank_after_source_node=corr(yr[:,jj],dr[:,0]),detected_partial_rank_after_source_node=corr(yr[:,jj],dr[:,1]),informative_gene_threshold=plan['informative_genes'],informative_gene_support_cells=n_support,informative_gene_support_fraction=n_support/len(md),signature_RNA_fraction_median=float(np.median(rna)),cell_scores_variable=bool(np.ptp(y[:,jj])>0),source_residual_scores_variable=bool(np.ptp(residual(y[:,jj],src))>0),source_limited=md.gse_id.nunique()<2))
   for sid in sorted(md.gse_id.unique()):
    which=md.gse_id.eq(sid).to_numpy();values=y[which,jj];within.append(dict(block=block,anchor_factor=int(j),scope=scope[j],major=major,gse_id=sid,cells=int(which.sum()),score_median=float(np.median(values)),score_q10=float(np.quantile(values,.1)),score_q90=float(np.quantile(values,.9)),RNA_rho=corr(values,depth[which,0]),detected_genes_rho=corr(values,depth[which,1]),informative_gene_support_cells=int((det[which]>=plan['informative_genes']).sum()),signature_RNA_fraction_median=float(np.median(rna[which]))))
  del raw;gc.collect()
 candidates=[];pframe=pd.DataFrame(pairs);sframe=pd.DataFrame(scores);dframe=pd.DataFrame(diagnostics)
 for j in range(len(scope)):
  js=[int(x[j]) for x in matched];cycle=maps[1,2][js[1]]==js[2];pm=[];sm=[]
  for a,b in combinations(range(3),2):
   pm.extend(pframe[(pframe.seed_a==seeds[a])&(pframe.seed_b==seeds[b])&(pframe.factor_a==js[a])].to_dict('records'));sm.extend(sframe[(sframe.seed_a==seeds[a])&(sframe.seed_b==seeds[b])&(sframe.factor_a==js[a])].to_dict('records'))
  mincos=min(r['cosine'] for r in pm);minjac=min(r['top50_jaccard'] for r in pm);weight_ok=cycle and mincos>=plan['min_cosine'] and minjac>=plan['min_top50_jaccard'];major_rows=[]
  majors=D['majors'] if scope[j]=='block_global' else [scope[j]]
  for major in majors:
   r=[row['source_residual_score_rho'] for row in sm if row['major']==major];assert len(r)==3;value=min(r) if np.isfinite(r).all() else np.nan;informative=support[j,major]>=plan['min_informative_cells'];major_rows.append(dict(major=major,min_rho=value,informative=informative,supported=bool(informative and np.isfinite(value) and value>=plan['min_source_residual_score_rho'])))
  supported=[v['major'] for v in major_rows if v['supported']];required=2 if scope[j]=='block_global' else 1;shared=bool(weight_ok and scope[j]=='block_global' and len(supported)>=2);context=bool(weight_ok and len(supported)>=1);main=shared if scope[j]=='block_global' else context
  dg=dframe[dframe.anchor_factor.eq(j)];complexity=bool((dg.RNA_partial_rank_after_source_node.abs()>=.5).any() or (dg.detected_partial_rank_after_source_node.abs()>=.5).any());names=list(ref['genes'][pooled_top[j]]);gs=set(names);prior_scope='global' if scope[j]=='block_global' else scope[j];prior=[]
  for name,genes in D['reference'][prior_scope].items():
   ss=set(genes);overlap=len(gs&ss)
   if overlap:prior.append((overlap/len(gs|ss),overlap,name))
  best=sorted(prior,key=lambda r:(-r[0],-r[1],r[2]))[0] if prior else (0,0,'none')
  candidates.append(dict(candidate_id='AS_'+block+'_'+str(j+1).zfill(3),block=block,anchor_factor=j,scope=scope[j],seed2_factor=js[1],seed3_factor=js[2],cycle_consistent=cycle,min_weight_cosine=mincos,min_top50_jaccard=minjac,weight_stability_pass=bool(weight_ok),required_major_support=required,supported_majors=';'.join(supported),supported_major_count=len(supported),stable_shared_or_lineage_candidate=bool(main),stable_context_candidate=context,stable_shared_candidate=shared,complexity_review=complexity,source_limited_block=len(set(meta.gse_id))<2,best_prior_descriptor=best[2],prior_top50_overlap=best[1],prior_descriptor_is_confirmed_state=False,primary50_median_L2_gene_weights=';'.join(names),all_major_score_support_diagnostics=json.dumps(major_rows),formal_state=False))
 O=R/'all_lineage_Spectra_diagnostics_20261014_v1'/block;assert not O.exists();O.mkdir(parents=True)
 for name,rows in [('all_factor_matches.tsv',pairs),('all_major_score_reproducibility.tsv',scores),('all_source_node_complexity_diagnostics.tsv',diagnostics),('all_source_major_distributions.tsv',within),('all_candidate_diagnostics.tsv',candidates)]:pd.DataFrame(rows).to_csv(O/name,sep='\t',index=False)
 assert len(pairs)==3*len(scope) and len(candidates)==len(scope)
 out=dict(status='COMPLETE_BLOCK_MULTISEED_DIAGNOSTICS_NOT_FORMAL_STATES',block=block,observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(P),input_seed_sha256=inputsha,candidates=len(scope),shared_or_lineage_stable_candidates=sum(c['stable_shared_or_lineage_candidate'] for c in candidates),context_stable_candidates=sum(c['stable_context_candidate'] for c in candidates),model_scope_all_majors_reported=True,global_factor_activity_not_required_in_all_majors=True,original_labels_unchanged=True,no_cell_level_p_or_clinical_LR_selection=True,formal_states=False,output_sha256={p.name:sha(p) for p in O.iterdir() if p.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
