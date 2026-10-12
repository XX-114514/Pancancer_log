"""Recount all masks/keys and independently check nuisance projection via lstsq.
No producer import; no response coefficient or p-value calculated.
"""
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def standard(a):
 a=np.asarray(a,float);s=np.sqrt(np.mean((a-a.mean(0))**2,axis=0));return (a-a.mean(0))/np.where(s>1e-12,s,1)
def independent_status(Z,x):
 coef,_,rank,_=np.linalg.lstsq(Z,x,rcond=1e-10);res=x-Z@coef
 if res@res<1e-10*len(x):return 'exposure_unidentified'
 X=np.column_stack((Z,x));_,_,fullrank,sv=np.linalg.lstsq(X,np.ones(len(x)),rcond=1e-10)
 if fullrank!=rank+1:return 'exposure_unidentified'
 if len(x)<=fullrank+2:return 'residual_df_insufficient'
 if sv[0]/sv[fullrank-1]>1e8:return 'condition_insufficient'
 leverage=np.sum(X*np.linalg.pinv(X,rcond=1e-10).T,axis=1)
 return 'high_leverage' if leverage.max()>.95 else 'design_estimable'
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';N=P/'major_pair_native_assembly_v1';T=P/'major_receiver_TF_preflight_v3';O=P/'major_LR_TF_design_preflight_v1';C=Path(os.environ['PANCANCER_SIGNALING_PRIOR']);plan=json.loads(Path(os.environ['PANCANCER_LR_TF_DESIGN_PLAN']).read_text());assert sha(__file__)==plan['auditor_sha256'] and not (O/'independent_audit.json').exists();inv=json.loads((O/'inventory.json').read_text())
 for n,h in inv['outputs'].items():assert sha(O/n)==h
 counts_path=AR/'20261010_expression_annotation_downstream_v1/inputs/sample_major_subcluster_counts.tsv';assert sha(counts_path)==plan['counts_sha256'];counts=pd.read_csv(counts_path,sep='\t');totals=counts.groupby('sampleID').n_cells.sum();pri=pd.read_csv(C/'private_all_receptor_TF_paths.tsv',sep='\t');assert sha(C/'private_all_receptor_TF_paths.tsv')==plan['prior_paths_sha256'];byR={r:g[['TF','expected_sign','min_edges']] for r,g in pri.groupby('receptor')};resource=pd.read_csv(Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']),sep='\t');contexts=pd.read_csv(T/'all_corrected_context_readiness.tsv',sep='\t',keep_default_na=False);summary=pd.read_csv(O/'all_context_model_readiness.tsv',sep='\t',keep_default_na=False).set_index('context_id');multisource=defaultdict(set);designsource=defaultdict(set);allcases=0;alldesigns=0
 for row in contexts.to_dict('records'):
  cid=row['context_id'];a=np.load(T/'private_context_arrays'/(cid+'.npz'));out=np.load(O/'private_context_arrays'/(cid+'.npz'));native=np.load(N/'private_group_arrays'/(row['group_id']+'.npz'));mask=native['corrected_material_for_next_inference']==row['material'];samples=a['samples'];TFs=a['TFs'];assert np.array_equal(samples,out['samples']) and np.array_equal(samples,native['samples'][mask]);tc={t:i for i,t in enumerate(TFs)};ready_tf=set(TFs[a['TF_ready']]);records=[]
  for li in np.flatnonzero(a['LR_RNA_support_ready']):
   r=resource.iloc[li].receptor
   if r in byR:
    for tf,sg,depth in byR[r].itertuples(index=False,name=None):
     if tf in ready_tf:records.append((int(li),tc[tf],int(sg),int(depth)))
  expected=set(records);observed=set(zip(out['candidate_LR_index'].tolist(),out['candidate_TF_index'].tolist(),out['prior_expected_sign'].tolist(),out['prior_min_edges'].tolist()));assert expected==observed and len(observed)==len(records);assert len(observed)==int(summary.loc[cid].candidate_family_cases);active=sorted(set(x[0] for x in records));ns=native['source_n'][mask];nt=native['target_n'][mask];base=np.column_stack((ns/totals.loc[samples],nt/totals.loc[samples],np.log1p(ns),np.log1p(nt),native['source_mean_detected_nonLR_genes'][mask],native['target_mean_detected_nonLR_genes'][mask]));pc=out['source_target_mixture_PCs'];assert pc.shape==(len(samples),4)
  # PCA independently verified with the sample covariance eigenvectors, sign-invariant.
  for offset,(major,n) in enumerate([(row['source'],ns),(row['target'],nt)]):
   g=counts[counts.majorCluster==major];f=g.pivot_table(index='sampleID',columns='subCluster',values='n_cells',aggfunc='sum',fill_value=0).reindex(samples,fill_value=0).to_numpy();assert np.array_equal(f.sum(1),n);ff=f/n[:,None];A=ff-ff.mean(0);ev,u=np.linalg.eigh(A@A.T);order=np.argsort(ev)[::-1];V=pc[:,offset*2:offset*2+2];positive=ev[order[:2]]>1e-20;rank=int(positive.sum());assert np.allclose(V.T@V,np.diag(np.maximum(ev[order[:2]],0)) if len(ev)>=2 else np.diag([max(ev[order[0]],0),0]),rtol=2e-6,atol=1e-10)
   if rank:
    U=u[:,order[:rank]];assert np.allclose(V-U@(U.T@V),0,atol=2e-6)
  names=out['status_names'].tolist();primary={};mixture={}
  for li in active:
   x=standard(native['ligand_meanlog'][mask,li]);Z=np.column_stack((np.ones(len(samples)),standard(np.column_stack((base,native['receptor_meanlog'][mask,li])))));Zm=np.column_stack((Z,standard(pc)));s=independent_status(Z,x);sm=independent_status(Zm,x);primary[li]=s;mixture[li]=sm;assert names[int(out['primary_design_status'][li])]==s and names[int(out['mixture_design_status'][li])]==sm;alldesigns+=2
  assert int(summary.loc[cid].primary_design_estimable_cases)==sum(primary[x[0]]=='design_estimable' for x in records);assert int(summary.loc[cid].mixture_design_estimable_cases)==sum(mixture[x[0]]=='design_estimable' for x in records)
  for li,tfi,sg,depth in records:
   key=(row['source'],row['target'],resource.iloc[li].ligand,resource.iloc[li].receptor,TFs[tfi]);multisource[key].add(row['gse_id'])
   if primary[li]=='design_estimable':designsource[key].add(row['gse_id'])
  allcases+=len(records)
 reported=pd.read_csv(O/'private_all_interface_TF_source_coverage.tsv.gz',sep='\t');assert len(reported)==len(multisource)
 for x in reported.itertuples(index=False):
  key=(x.source,x.target,x.ligand,x.receptor,x.TF);assert len(multisource[key])==x.n_GSE_RNA_TF_eligible and len(designsource[key])==x.n_GSE_primary_design_estimable
 assert allcases==inv['full_frozen_candidate_family_cases'];j=dict(status='PASS_ALL_CANDIDATE_MASKS_DESIGN_STATUSES_MIXTURE_SUBSPACES_AND_SOURCE_COUNTS',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),all_family_cases_checked=allcases,all_LR_designs_checked=alldesigns,all_contexts_checked=len(contexts),all_source_interface_keys_checked=len(multisource),no_effects_p_or_functional_validation=True);(O/'independent_audit.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j))
if __name__=='__main__':main()
