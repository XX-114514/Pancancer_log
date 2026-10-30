"""Frozen CM network edges within source/cancer/material/treatment contexts.

Cluster-bootstrap intervals describe conditional composition associations.
Candidate material blocks are not certified patients; discovery is reused.
"""
import argparse,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from TME_block_rank_bootstrap_v1 import block_multiplicities,weighted_spearman

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bound(b):
 p=Path(b['path']);assert sha(p)==b['SHA256'];return p
def binding(p):return dict(path=str(p),SHA256=sha(p))
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,indent=2,allow_nan=False)
def wait(root,start,seconds):
 p=Path(root)/'provenance/supervisor_terminal_v1.json'
 while not p.exists():
  assert time.monotonic()-start<seconds;time.sleep(5)
 assert load(p)['actual_child_exit_code']==0 and load(p).get('science_terminal_exists',True)
def meta_summary(rows):
 rho=np.array([r['Spearman'] for r in rows],float);z=np.arctanh(np.clip(rho,-1+1e-7,1-1e-7));mu=float(np.tanh(z.mean()));n=len(z)
 leave=[float(np.tanh(np.delete(z,i).mean())) for i in range(n)] if n>=2 else []
 variances=np.array([r['bootstrap_Fisher_z_variance'] if r['bootstrap_Fisher_z_variance'] is not None else np.nan for r in rows]);valid=np.isfinite(variances)&(variances>0);I2=None
 if valid.sum()>=2:
  v=variances[valid];zz=z[valid];w=1/v;center=np.sum(zz*w)/w.sum();Q=float(np.sum(w*(zz-center)**2));I2=max(0.,(Q-(len(zz)-1))/Q) if Q>0 else 0.
 return dict(n_sources=n,equal_source_Fisher_mean_Spearman=mu,n_positive_sources=int((rho>0).sum()),n_negative_sources=int((rho<0).sum()),leave_one_source_Spearman=leave,all_leave_one_direction_matches=bool(leave and all(np.sign(x)==np.sign(mu) for x in leave)),I_squared_descriptive=I2,bootstrap_source_errors_not_independent_patient_certification=True)

def main(R):
 cfg=load(R/'config/spec_v1.json');start=time.monotonic();wait(cfg['math_test_stage'],start,cfg['dependency_wait_seconds']);test=Path(cfg['math_test_stage'])/'outputs/accepted_block_rank_bootstrap_math_v1.json';assert load(test)['all_tests_passed'];assert sha(Path(cfg['math_test_stage'])/'scripts/TME_block_rank_bootstrap_v1.py')==sha(R/'scripts/TME_block_rank_bootstrap_v1.py')
 sys.path.insert(0,str(Path(cfg['joined'])/'scripts'));from extract_full_source_local_marker_evidence_v1 import storage_gate
 storage_gate(R,cfg);results=[];origins=[];metas=[];replays=[];base_context_columns=['source','cancer_type','material_class','treatment_phase']
 for scheme in cfg['schemes']:
  root=Path(scheme['root']);assert load(root/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
  ac=load(root/'config/spec_v1.json')['input_acceptance'];a=load(bound(ac));fq=pd.read_csv(bound(a['frequency']),sep='\t',index_col=0);meta=pd.read_csv(bound(a['metadata']),sep='\t');assert fq.columns[1:].tolist()==meta.sampleID.tolist();memb=load(bound(a['annotation_acceptance']));node={n['name']:n for n in memb['nodes']}
  for c in base_context_columns:meta[c]=meta[c].fillna('unknown').astype(str)
  edges={}
  for seed in cfg['seeds']:
   for flavor in ['author_compatibility','raw_Pearson_BH_sensitivity']:
    f=root/'outputs'/('seed_'+str(seed))/(flavor+'_each_edge_v1.tsv');tab=pd.read_csv(f,sep='\t')
    for r in tab.to_dict('records'):
     pair=tuple(sorted([r['subCluster1'],r['subCluster2']]));edges.setdefault(pair,[]).append(dict(seed=seed,flavor=flavor,module=r['cm']))
  for pair,origin in sorted(edges.items()):origins.append(dict(scheme=scheme['id'],node_a=pair[0],node_b=pair[1],origins=origin,duplicate_shared_module_records_count_once=True))
  groups=[]
  for context,positions in meta.groupby(base_context_columns,sort=True,dropna=False).indices.items():
   ix=np.asarray(positions,int);key=dict(zip(base_context_columns,context));nb=meta.iloc[ix].patient_candidate_block.nunique();basis=meta.iloc[ix].patient_basis_available.fillna(False).astype(str).str.lower().isin(['true','1','yes']).sum()
   seed=cfg['seed']+int(hashlib.sha256('|'.join(context).encode()).hexdigest()[:8],16);weights=block_multiplicities(meta.iloc[ix].patient_candidate_block.to_numpy(),cfg['n_bootstrap'],seed) if len(ix)>=cfg['minimum_units'] and nb>=cfg['minimum_blocks'] else None
   groups.append((key,ix,nb,int(basis),weights))
  for pair in sorted(edges):
   x=fq.loc[pair[0]].iloc[1:].to_numpy(float);y=fq.loc[pair[1]].iloc[1:].to_numpy(float);edge_rows=[]
   for key,ix,nb,basis,weights in groups:
    r=dict(scheme=scheme['id'],node_a=pair[0],node_b=pair[1],major_a=fq.loc[pair[0]].iloc[0],major_b=fq.loc[pair[1]].iloc[0],**key,n_work_units=len(ix),n_candidate_blocks=nb,n_units_with_patient_basis_field=basis,patient_basis_field_is_not_manual_patient_review=True,same_parent_compositional_closure_applies=fq.loc[pair[0]].iloc[0]==fq.loc[pair[1]].iloc[0])
    if weights is None:r.update(status='insufficient_units_or_candidate_blocks',Spearman=None)
    elif np.ptp(x[ix])==0 or np.ptp(y[ix])==0:r.update(status='no_within_context_variation',Spearman=None)
    else:
     rho=float(spearmanr(x[ix],y[ix]).statistic);rr=weighted_spearman(x[ix],y[ix],weights);valid=rr[np.isfinite(rr)];enough=len(valid)>=np.ceil(cfg['minimum_finite_bootstrap_fraction']*cfg['n_bootstrap']);ci=np.quantile(valid,[.025,.975]).tolist() if enough else [None,None]
     z=np.arctanh(np.clip(valid,-1+1e-7,1-1e-7));r.update(status='estimable' if enough else 'insufficient_variable_bootstrap_replicates',Spearman=rho,CI_low=ci[0],CI_high=ci[1],n_finite_bootstraps=len(valid),bootstrap_Fisher_z_variance=float(z.var(ddof=1)) if len(z)>1 else None,n_observed_zero_a=int((x[ix]==0).sum()),n_observed_zero_b=int((y[ix]==0).sum()))
     # Independently expand material rows for the first valid positive/negative
     # case, in addition to the synthetic unequal-block/tie counterexamples.
     sign='positive' if rho>0 else 'negative' if rho<0 else 'zero'
     if not any(v['sign']==sign for v in replays):
      for b in range(min(3,len(weights))):
       expanded=np.repeat(np.arange(len(ix)),weights[b]);want=spearmanr(x[ix][expanded],y[ix][expanded]).statistic;np.testing.assert_allclose(rr[b],want,equal_nan=True,atol=1e-12)
      replays.append(dict(sign=sign,expanded_material_bootstrap_replay=True))
    results.append(r);edge_rows.append(r)
   valid=[r for r in edge_rows if r['status']=='estimable']
   # Same cancer, material and phase first. Different cancer support is kept
   # as a distinct cross-context summary, never called same-cancer replication.
   frame=pd.DataFrame(valid)
   if len(frame):
    for context,group in frame.groupby(['cancer_type','material_class','treatment_phase'],sort=True,dropna=False):
     rows=group.to_dict('records')
     if len(rows)>=2:metas.append(dict(scheme=scheme['id'],node_a=pair[0],node_b=pair[1],scope='same_recorded_cancer_material_phase',context=dict(zip(['cancer_type','material_class','treatment_phase'],context)),unknown_context_not_certified_as_same_biology=any(str(v).lower() in {'unknown','unavailable','not_reported','nan','na',''} for v in context),**meta_summary(rows)))
    # At most one equal contribution per source for the broader summary.
    bysource=[]
    for source,group in frame.groupby('source',sort=True):
     rs=group.to_dict('records');zz=np.arctanh(np.clip([v['Spearman'] for v in rs],-1+1e-7,1-1e-7));bysource.append(dict(source=source,Spearman=float(np.tanh(zz.mean())),bootstrap_Fisher_z_variance=None))
    if len(bysource)>=2:metas.append(dict(scheme=scheme['id'],node_a=pair[0],node_b=pair[1],scope='multiple_sources_across_recorded_contexts',n_recorded_cancers=frame.cancer_type.nunique(),not_same_cancer_replication=True,**meta_summary(bysource)))
  save(R/'checkpoints'/('scheme_'+scheme['id']+'_complete_v1.json'),dict(n_unique_edges=len(edges),n_rows=len(results)));storage_gate(R,cfg)
 table=R/'outputs/private_CM_within_source_context_effects_v1.tsv';pd.DataFrame(results).to_csv(table,sep='\t',index=False);mp=R/'outputs/private_CM_cross_source_context_summary_v1.json';save(mp,metas)
 op=R/'outputs/accepted_CM_source_context_recurrence_v1.json';save(op,dict(schema='TME_CM_frozen_edge_source_context_recurrence_v1',n_schemes=len(cfg['schemes']),n_unique_edges=len(origins),n_context_effect_records=len(results),n_estimable_records=sum(r['status']=='estimable' for r in results),effects=binding(table),cross_source_summaries=binding(mp),edge_module_origins=origins,exact_bootstrap_math_acceptance=binding(test),direct_expanded_real_case_replays=replays,n_bootstrap=cfg['n_bootstrap'],conditional_on_discovery_selected_edges=True,no_independent_discovery_or_patient_validation_claim=True,unknown_source_design_and_patient_review_limits_not_removed_by_adjustment=True,clinical_outcomes_used=False))
 storage_gate(R,cfg);save(R/'provenance/terminal_v1.json',dict(status='FROZEN_CM_EDGE_SOURCE_CONTEXT_COMPOSITION_RECURRENCE_COMPLETE',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),acceptance_SHA256=sha(op)))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);R=Path(p.parse_args().output)
 try:main(R)
 except BaseException:
  import traceback
  save(R/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
