"""Freeze automatic candidate families and nuisance spans, no association p.
All 318 corrected contexts; measured eligibility stays separate from estimability.
"""
from pathlib import Path
from collections import defaultdict
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def z(x):
 a=np.asarray(x,float);s=a.std(0);return (a-a.mean(0))/np.where(s>1e-12,s,1)
def pcs(fractions):
 a=fractions-fractions.mean(0);u,s,v=np.linalg.svd(a,full_matrices=False);o=np.zeros((len(a),2))
 for k in range(min(2,len(s))):
  if s[k]>1e-10:
   sign=1 if v[k,np.argmax(abs(v[k]))]>=0 else -1;o[:,k]=u[:,k]*s[k]*sign
 return o

def design(Z,x):
 Q,s,_=np.linalg.svd(Z,full_matrices=False);rank=int((s>s[0]*1e-10).sum());Q=Q[:,:rank];xr=x-Q@(Q.T@x);ss=float(xr@xr)
 if ss<1e-10*len(x):return 'exposure_unidentified',rank,None
 full=np.column_stack([Z,x]);U,sv,_=np.linalg.svd(full,full_matrices=False);rr=int((sv>sv[0]*1e-10).sum());hat=np.sum(U[:,:rr]**2,axis=1);mh=float(hat.max())
 if rr!=rank+1:return 'exposure_unidentified',rank,mh
 if len(x)<=rr+2:return 'residual_df_insufficient',rank,mh
 if sv[0]/sv[rr-1]>1e8:return 'condition_insufficient',rank,mh
 if mh>.95:return 'high_leverage',rank,mh
 return 'design_estimable',rank,mh

def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';N=P/'major_pair_native_assembly_v1';T=P/'major_receiver_TF_preflight_v3';O=P/'major_LR_TF_design_preflight_v1';C=Path(os.environ['PANCANCER_SIGNALING_PRIOR']);planp=Path(os.environ['PANCANCER_LR_TF_DESIGN_PLAN']);plan=json.loads(planp.read_text());assert sha(__file__)==plan['producer_sha256'] and not O.exists();start=datetime.now(timezone.utc).isoformat()
 for f,h in [(N/'inventory.json',plan['native_inventory_sha256']),(T/'inventory.json',plan['TF_inventory_sha256']),(T/'independent_audit.json',plan['TF_audit_sha256']),(C/'inventory.json',plan['prior_inventory_sha256']),(C/'independent_audit.json',plan['prior_audit_sha256']),(C/'private_all_receptor_TF_paths.tsv',plan['prior_paths_sha256'])]:assert sha(f)==h
 ti=json.loads((T/'inventory.json').read_text());ni=json.loads((N/'inventory.json').read_text());counts_path=AR/'20261010_expression_annotation_downstream_v1/inputs/sample_major_subcluster_counts.tsv';assert sha(counts_path)==plan['counts_sha256'];counts=pd.read_csv(counts_path,sep='\t');totals=counts.groupby('sampleID').n_cells.sum();fractions={}
 for major,g in counts.groupby('majorCluster'):fractions[major]=g.pivot_table(index='sampleID',columns='subCluster',values='n_cells',aggfunc='sum',fill_value=0)
 resource=pd.read_csv(Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']),sep='\t');assert sha(Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']))==plan['LR_sha256'];prior=pd.read_csv(C/'private_all_receptor_TF_paths.tsv',sep='\t');TFs=np.load(T/'private_source_TF_arrays/GSE116256.npz')['TFs'];tc={t:i for i,t in enumerate(TFs)};byR={r:g for r,g in prior.groupby('receptor')};global_pairs=[]
 for li,r in enumerate(resource.receptor):
  if r in byR:
   for p in byR[r].itertuples(index=False):global_pairs.append((li,tc[p.TF],p.expected_sign,p.min_edges))
 gp=np.array(global_pairs,dtype=np.int32);assert len(set((int(x),int(y)) for x,y in gp[:,:2]))==len(gp);contexts=pd.read_csv(T/'all_corrected_context_readiness.tsv',sep='\t',keep_default_na=False);O.mkdir();(O/'private_context_arrays').mkdir();summary=[];multi=defaultdict(set);designmulti=defaultdict(set);all_family_cases=0
 for row in contexts.to_dict('records'):
  cid=row['context_id'];f=T/'private_context_arrays'/(cid+'.npz');assert sha(f)==ti['outputs'][str(f.relative_to(T))];a=np.load(f);nf=N/'private_group_arrays'/(row['group_id']+'.npz');assert sha(nf)==ni['outputs'][str(nf.relative_to(N))];native=np.load(nf);mask=native['corrected_material_for_next_inference']==row['material'];assert np.array_equal(a['samples'],native['samples'][mask]);samples=a['samples'];chosen=gp[a['LR_RNA_support_ready'][gp[:,0]]&a['TF_ready'][gp[:,1]]];Nlib=len(samples);primary=np.full(4624,-1,dtype=np.int8);mixture=primary.copy();status_names=['design_estimable','exposure_unidentified','residual_df_insufficient','condition_insufficient','high_leverage'];active_lr=np.unique(chosen[:,0]) if len(chosen) else np.array([],int);source_n=native['source_n'][mask];target_n=native['target_n'][mask]
  base=np.column_stack([source_n/totals.loc[samples].to_numpy(),target_n/totals.loc[samples].to_numpy(),np.log1p(source_n),np.log1p(target_n),native['source_mean_detected_nonLR_genes'][mask],native['target_mean_detected_nonLR_genes'][mask]])
  mix=[]
  for major,n in [(row['source'],source_n),(row['target'],target_n)]:
   f=fractions[major].reindex(samples,fill_value=0).to_numpy();assert np.array_equal(f.sum(1),n);mix.append(pcs(f/n[:,None]))
  pc=np.column_stack(mix);ranks=[];hats=[]
  for li in active_lr:
   x=z(native['ligand_meanlog'][mask,li]);Z=np.column_stack([np.ones(Nlib),z(np.column_stack([base,native['receptor_meanlog'][mask,li]]))]);Zm=np.column_stack([Z,z(pc)]);s,rank,hat=design(Z,x);sm,rm,hm=design(Zm,x);primary[li]=status_names.index(s);mixture[li]=status_names.index(sm);ranks.append((int(li),rank,rm));hats.append((int(li),np.nan if hat is None else hat,np.nan if hm is None else hm))
  np.savez_compressed(O/'private_context_arrays'/(cid+'.npz'),candidate_LR_index=chosen[:,0],candidate_TF_index=chosen[:,1],prior_expected_sign=chosen[:,2].astype(np.int8),prior_min_edges=chosen[:,3].astype(np.int8),primary_design_status=primary,mixture_design_status=mixture,source_target_mixture_PCs=pc,design_ranks=np.array(ranks,dtype=int).reshape(-1,3),design_max_leverages=np.array(hats,float).reshape(-1,3),samples=samples,status_names=np.array(status_names,dtype=str))
  result=row|dict(candidate_family_cases=len(chosen),unique_LR_exposures=len(active_lr),prior_sign_ambiguous_cases=int((chosen[:,2]==0).sum()) if len(chosen) else 0,primary_design_estimable_cases=int((primary[chosen[:,0]]==0).sum()),mixture_design_estimable_cases=int((mixture[chosen[:,0]]==0).sum()))
  for k,name in enumerate(status_names):result['primary_'+name+'_LRs']=int((primary[active_lr]==k).sum());result['mixture_'+name+'_LRs']=int((mixture[active_lr]==k).sum())
  summary.append(result);all_family_cases+=len(chosen)
  for li,tfi,sg,dep in chosen:
   key=(row['source'],row['target'],int(li),int(tfi));multi[key].add(row['gse_id'])
   if primary[li]==0:designmulti[key].add(row['gse_id'])
  print(json.dumps(dict(context=cid,family_cases=len(chosen),primary_estimable_cases=result['primary_design_estimable_cases'],mixture_estimable_cases=result['mixture_design_estimable_cases'])),flush=True)
 pd.DataFrame(summary).to_csv(O/'all_context_model_readiness.tsv',sep='\t',index=False);rows=[]
 for key,sources in sorted(multi.items()):
  source,target,li,tfi=key;rows.append(dict(source=source,target=target,ligand=resource.iloc[li].ligand,receptor=resource.iloc[li].receptor,TF=TFs[tfi],n_GSE_RNA_TF_eligible=len(sources),n_GSE_primary_design_estimable=len(designmulti[key])))
 pd.DataFrame(rows).to_csv(O/'private_all_interface_TF_source_coverage.tsv.gz',sep='\t',index=False);public=[]
 for (source,target),g in pd.DataFrame(rows).groupby(['source','target']):public.append(dict(source=source,target=target,unique_LR_TF_cases=len(g),RNA_TF_eligible_ge2GSE=int((g.n_GSE_RNA_TF_eligible>=2).sum()),RNA_TF_eligible_ge3GSE=int((g.n_GSE_RNA_TF_eligible>=3).sum()),primary_design_estimable_ge2GSE=int((g.n_GSE_primary_design_estimable>=2).sum()),primary_design_estimable_ge3GSE=int((g.n_GSE_primary_design_estimable>=3).sum())))
 pd.DataFrame(public).to_csv(O/'cell_pair_multisource_model_readiness.tsv',sep='\t',index=False);j=dict(status='COMPLETE_DESIGN_PREFLIGHT_PENDING_INDEPENDENT_RECOUNT',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),external_global_LR_TF_pairs=len(gp),context_roster=len(contexts),full_frozen_candidate_family_cases=all_family_cases,primary_design_estimable_cases=sum(x['primary_design_estimable_cases'] for x in summary),mixture_design_estimable_cases=sum(x['mixture_design_estimable_cases'] for x in summary),unique_cellpair_LR_TF_interfaces=len(rows),interfaces_RNA_TF_eligible_ge2GSE=sum(len(v)>=2 for v in multi.values()),interfaces_primary_design_estimable_ge2GSE=sum(len(v)>=2 for v in designmulti.values()),no_effects_p_or_validation=True,no_constant_covariate_deleted=True,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()});(O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps({k:v for k,v in j.items() if k!='outputs'}))
if __name__=='__main__':main()
