"""Automatic full-family covariance-aware recheck; exploratory same-data validation.
Every nominal main LR p<.05, all matched TFs and all original controls retained.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,importlib.util,gc
import numpy as np,pandas as pd
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';N=P/'major_pair_native_assembly_v1';T=P/'major_receiver_TF_preflight_v3';D=P/'major_LR_TF_design_preflight_v1';S=P/'major_LR_receiver_TF_screen_v1';O=P/'major_LR_receiver_TF_wild_v1';planp=Path(os.environ['PANCANCER_LR_TF_WILD_PLAN']);plan=json.loads(planp.read_text());assert sha(__file__)==plan['producer_sha256'] and not O.exists();cp=Path(os.environ['PANCANCER_CR2_CORE']);wp=Path(os.environ['PANCANCER_WILD_CORE']);assert sha(cp)==plan['CR2_core_sha256'] and sha(wp)==plan['wild_core_sha256'];core=load('cr2',cp);wild=load('wild',wp)
 for p,h in [(S/'inventory.json',plan['screen_inventory_sha256']),(S/'independent_numeric_audit.json',plan['screen_audit_sha256']),(S/'all177393_LR_union_models.tsv.gz',plan['screen_table_sha256'])]:assert sha(p)==h
 table=pd.read_csv(S/'all177393_LR_union_models.tsv.gz',sep='\t',na_values=['NA'],keep_default_na=False);selected=table[(table.variant=='abundance_complexity')&(table.joint_p_exploratory<.05)].copy();assert len(selected)==2023;assert not selected.duplicated(['context_id','LR_index']).any();contexts=pd.read_csv(D/'all_context_model_readiness.tsv',sep='\t',keep_default_na=False).set_index('context_id');counts_path=AR/'20261010_expression_annotation_downstream_v1/inputs/sample_major_subcluster_counts.tsv';assert sha(counts_path)==plan['counts_sha256'];totals=pd.read_csv(counts_path,sep='\t').groupby('sampleID').n_cells.sum();O.mkdir();(O/'private_weight_audit_cases').mkdir();out=table.copy();out['wild_status']='not_selected_exploratory_screen';out['wild_joint_p']=np.nan;out['wild_draws']=np.nan;out['wild_extreme']=np.nan;out['wild_enumerated_sign_pairs']=False;out['wild_nonfinite_statistics']=np.nan;started=datetime.now(timezone.utc).isoformat();rows=[];case_examples=0
 for cid,g in selected.groupby('context_id',sort=True):
  row=contexts.loc[cid];a=np.load(D/'private_context_arrays'/(cid+'.npz'));b=np.load(T/'private_context_arrays'/(cid+'.npz'));nat=np.load(N/'private_group_arrays'/(row.group_id+'.npz'));mask=nat['corrected_material_for_next_inference']==row.material;samples=a['samples'];assert np.array_equal(samples,b['samples']);blocks=b['blocks'];ns=nat['source_n'][mask];nt=nat['target_n'][mask];n=len(samples);base=np.column_stack((ns/totals.loc[samples],nt/totals.loc[samples],np.log1p(ns),np.log1p(nt),nat['source_mean_detected_nonLR_genes'][mask],nat['target_mean_detected_nonLR_genes'][mask]))
  for entry in g.itertuples(index=False):
   li=int(entry.LR_index);ix=np.flatnonzero(a['candidate_LR_index']==li);tfidx=a['candidate_TF_index'][ix];x=core.standardize(nat['ligand_meanlog'][mask,li]);Z=np.column_stack((np.ones(n),core.standardize(np.column_stack((base,nat['receptor_meanlog'][mask,li])))));Zm=np.column_stack((Z,core.standardize(a['source_target_mixture_PCs'])))
   for variant in ['abundance_complexity','subtype_mixture','full_gene_sensitivity']:
    index=out.index[(out.context_id==cid)&(out.LR_index==li)&(out.variant==variant)];assert len(index)==1;i=int(index[0]);state=out.loc[i,'status'];result=dict(context_id=cid,LR_index=li,variant=variant,matched_TFs=len(tfidx),state='screen_'+state)
    if state=='inference_estimable':
     ZZ=Zm if variant=='subtype_mixture' else Z;cr=core.prepare_cr2(ZZ,x,blocks,min_df=4);assert cr['status']=='inference_estimable';Y=core.standardize(b['sensitivity_TF_scores' if variant=='full_gene_sensitivity' else 'primary_TF_scores'][:,tfidx]);stats=core.estimate_cr2(cr,Y);prep=wild.prepare(ZZ,x,blocks,cr)
     if prep['status']=='prepared':
      seed=int(hashlib.sha256((cid+'|'+str(li)+'|'+variant+'|wildmax20261012').encode()).hexdigest()[:16],16);rng=np.random.default_rng(seed);res=wild.test(prep,Y,stats['stat'],rng,draws=plan['bootstrap_draws'],exact_G_threshold=plan['enumerate_G_le'],batch=256);result.update(state='wild_rechecked_exploratory',seed=str(seed),**res);out.loc[i,['wild_joint_p','wild_draws','wild_extreme','wild_enumerated_sign_pairs','wild_nonfinite_statistics']]=[res['p'],res['draws'],res['extreme'],res['enumerated_all_sign_pairs'],res['nonfinite_bootstrap_statistics']]
      # Deterministic first eligible case per context, capped only audit fixtures, not science.
      if variant=='abundance_complexity' and case_examples<plan['audit_fixture_cases']:
       weights=np.random.default_rng(seed^314159).choice(np.array([-1.,1.]),size=(32,len(np.unique(blocks))));fast=wild.bootstrap_max_statistics(prep,Y,weights);np.savez_compressed(O/'private_weight_audit_cases'/(cid+'_LR'+str(li)+'.npz'),Z=ZZ,x=x,blocks=blocks,Y=Y,weights=weights,fast_max_statistics=fast);case_examples+=1
     else:result['state']=prep['status']
    out.loc[i,'wild_status']=result['state'];rows.append(result);print(json.dumps({k:v for k,v in result.items() if k!='finite_sample_method'}),flush=True)
  gc.collect()
 assert len(rows)==3*2023
 for v,g in out.groupby('variant'):
  p=g.wild_joint_p.fillna(1).to_numpy(float);q=core.bh(p);valid=g.wild_status.eq('wild_rechecked_exploratory').to_numpy();out.loc[g.index,'wild_BH_exploratory']=np.where(valid,q,np.nan);out.loc[g.index,'wild_BY_exploratory']=np.where(valid,np.minimum(1,q*sum(1/i for i in range(1,len(g)+1))),np.nan)
 out.to_csv(O/'all177393_LR_wild_models.tsv.gz',sep='\t',index=False,na_rep='NA');pd.DataFrame(rows).to_csv(O/'all6069_selected_recheck_records.tsv.gz',sep='\t',index=False,na_rep='NA');summary=[]
 for (variant,sid),g in out.groupby(['variant','gse_id']):summary.append(dict(variant=variant,gse_id=sid,original_family_cases=len(g),wild_rechecked=int(g.wild_status.eq('wild_rechecked_exploratory').sum()),BH05_exploratory=int((g.wild_BH_exploratory<.05).sum()),BY05_exploratory=int((g.wild_BY_exploratory<.05).sum())))
 pd.DataFrame(summary).to_csv(O/'source_wild_summary.tsv',sep='\t',index=False);j=dict(status='COMPLETE_COVARIANCE_AWARE_RECHECK_PENDING_NUMERIC_AUDIT',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),selected_main_LR_cases=2023,recheck_rows=len(rows),original_family_each_variant=59131,summary=summary,independent_source_or_binding_validation=False,adaptive_screening_same_data=True,all_same_TFs_controls_endpoints_retained=True,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()});(O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps({k:v for k,v in j.items() if k!='outputs'}))
if __name__=='__main__':main()
