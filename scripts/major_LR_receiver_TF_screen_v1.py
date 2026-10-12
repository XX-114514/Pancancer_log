"""All frozen LR/context families: exploratory CR2 union over receiver TFs.
Not confirmed FDR or receptor activation. Follow-ups need block resampling/source checks.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,importlib.util,gc
import numpy as np,pandas as pd
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';N=P/'major_pair_native_assembly_v1';T=P/'major_receiver_TF_preflight_v3';D=P/'major_LR_TF_design_preflight_v1';O=P/'major_LR_receiver_TF_screen_v1';planp=Path(os.environ['PANCANCER_LR_TF_SCREEN_PLAN']);plan=json.loads(planp.read_text());assert sha(__file__)==plan['producer_sha256'] and not O.exists();corep=Path(os.environ['PANCANCER_CR2_CORE']);assert sha(corep)==plan['core_sha256'];spec=importlib.util.spec_from_file_location('core',corep);core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
 for f,h in [(D/'inventory.json',plan['design_inventory_sha256']),(D/'independent_audit.json',plan['design_audit_sha256']),(N/'inventory.json',plan['native_inventory_sha256']),(T/'inventory.json',plan['TF_inventory_sha256'])]:assert sha(f)==h
 di=json.loads((D/'inventory.json').read_text());ni=json.loads((N/'inventory.json').read_text());ti=json.loads((T/'inventory.json').read_text());counts_path=AR/'20261010_expression_annotation_downstream_v1/inputs/sample_major_subcluster_counts.tsv';assert sha(counts_path)==plan['counts_sha256'];counts=pd.read_csv(counts_path,sep='\t');totals=counts.groupby('sampleID').n_cells.sum();resourcep=Path(os.environ['PANCANCER_FROZEN_LR_RESOURCE']);assert sha(resourcep)==plan['LR_resource_sha256'];resource=pd.read_csv(resourcep,sep='\t');contexts=pd.read_csv(D/'all_context_model_readiness.tsv',sep='\t',keep_default_na=False);O.mkdir();(O/'private_component_effects').mkdir();rows=[];start=datetime.now(timezone.utc).isoformat();variants=['abundance_complexity','subtype_mixture','full_gene_sensitivity']
 for row in contexts.to_dict('records'):
  cid=row['context_id'];af=D/'private_context_arrays'/(cid+'.npz');bf=T/'private_context_arrays'/(cid+'.npz');nf=N/'private_group_arrays'/(row['group_id']+'.npz')
  for f,inv,root in [(af,di,D),(bf,ti,T),(nf,ni,N)]:assert sha(f)==inv['outputs'][str(f.relative_to(root))]
  a=np.load(af);b=np.load(bf);nat=np.load(nf);mask=nat['corrected_material_for_next_inference']==row['material'];assert np.array_equal(a['samples'],b['samples']) and np.array_equal(a['samples'],nat['samples'][mask]);blocks=b['blocks'];n=len(blocks);sn=nat['source_n'][mask];tn=nat['target_n'][mask];base=np.column_stack((sn/totals.loc[a['samples']],tn/totals.loc[a['samples']],np.log1p(sn),np.log1p(tn),nat['source_mean_detected_nonLR_genes'][mask],nat['target_mean_detected_nonLR_genes'][mask]));case_lr=a['candidate_LR_index'];case_tf=a['candidate_TF_index'];casecount=len(case_lr);effects={v:{k:np.full(casecount,np.nan) for k in ['beta','se','component_p_diagnostic','lower_marginal','upper_marginal','df']} for v in variants};states={v:[] for v in variants};TFs=b['TFs'];active=np.unique(case_lr);meta={k:row[k] for k in ['context_id','group_id','gse_id','stratum','cancerType','sample_class','source','target','material','n_libraries','n_analysis_blocks','n_reviewed_blocks']}
  for li in active:
   where=np.flatnonzero(case_lr==li);tfidx=case_tf[where];K=len(tfidx);x=core.standardize(nat['ligand_meanlog'][mask,li]);Z=np.column_stack((np.ones(n),core.standardize(np.column_stack((base,nat['receptor_meanlog'][mask,li])))));Zm=np.column_stack((Z,core.standardize(a['source_target_mixture_PCs'])));nativeY=b['primary_TF_scores'][:,tfidx];fullY=b['sensitivity_TF_scores'][:,tfidx];prepared={}
   for variant in variants:
    dstatus=int(a['mixture_design_status' if variant=='subtype_mixture' else 'primary_design_status'][li]);result=meta|dict(variant=variant,LR_index=int(li),ligand=resource.iloc[li].ligand,receptor=resource.iloc[li].receptor,n_compatible_TFs=K,status='design_'+a['status_names'][dstatus],effective_df=np.nan,max_cluster_hat=np.nan,joint_p_exploratory=np.nan,lead_TF='',lead_beta=np.nan,lead_se=np.nan,lead_marginal_lower=np.nan,lead_marginal_upper=np.nan,lead_prior_expected_sign=np.nan,joint_family_size=59131,confirmed_communication=False)
    if dstatus==0:
     key='mixture' if variant=='subtype_mixture' else 'base'
     if key not in prepared:prepared[key]=core.prepare_cr2(Zm if key=='mixture' else Z,x,blocks,min_df=plan['min_effective_df'])
     prep=prepared[key];result['status']=prep['status'];result['max_cluster_hat']=prep.get('max_cluster_hat',np.nan);result['effective_df']=prep.get('df',np.nan)
     if 'coef_vector' in prep:
      Y=core.standardize(fullY if variant=='full_gene_sensitivity' else nativeY);assert np.isfinite(Y).all();e=core.estimate_cr2(prep,Y)
      for k,c in [('beta','beta'),('se','se'),('component_p_diagnostic','p'),('lower_marginal','lower'),('upper_marginal','upper')]:effects[variant][k][where]=e[c]
      effects[variant]['df'][where]=prep['df']
      if prep['status']=='inference_estimable':
       p=np.where(np.isfinite(e['p']),e['p'],1.);result['joint_p_exploratory']=float(min(1,K*p.min()));lead=int(np.argmin(p));result.update(lead_TF=str(TFs[tfidx[lead]]),lead_beta=float(e['beta'][lead]),lead_se=float(e['se'][lead]),lead_marginal_lower=float(e['lower'][lead]),lead_marginal_upper=float(e['upper'][lead]),lead_prior_expected_sign=int(a['prior_expected_sign'][where[lead]]))
       if not np.isfinite(e['p']).any():result['status']='outcome_residual_variance_insufficient';result['joint_p_exploratory']=np.nan
    rows.append(result);states[variant].append(result['status'])
  payload=dict(candidate_LR_index=case_lr,candidate_TF_index=case_tf,TFs=TFs,prior_expected_sign=a['prior_expected_sign'])
  for v in variants:
   for k,values in effects[v].items():payload[v+'__'+k]=values
  np.savez_compressed(O/'private_component_effects'/(cid+'.npz'),**payload);print(json.dumps(dict(context=cid,LR_units=len(active),inference_estimated={v:states[v].count('inference_estimable') for v in variants})),flush=True);del effects,payload;gc.collect()
 out=pd.DataFrame(rows);assert len(out)==3*59131 and not out.duplicated(['context_id','LR_index','variant']).any()
 for v,g in out.groupby('variant',sort=False):
  assert len(g)==59131;p=g.joint_p_exploratory.fillna(1).to_numpy();q=core.bh(p);by=core.bh(p)*sum(1/i for i in range(1,len(g)+1));good=g.status.eq('inference_estimable').to_numpy();out.loc[g.index,'joint_BH_exploratory']=np.where(good,q,np.nan);out.loc[g.index,'joint_BY_exploratory']=np.where(good,np.minimum(1,by),np.nan)
 out.to_csv(O/'all177393_LR_union_models.tsv.gz',sep='\t',index=False,na_rep='NA');summary=[]
 for (v,sid),g in out.groupby(['variant','gse_id']):summary.append(dict(variant=v,gse_id=sid,LR_units=len(g),inference_estimated=int(g.status.eq('inference_estimable').sum()),BH05_exploratory=int((g.joint_BH_exploratory<.05).sum()),BY05_exploratory=int((g.joint_BY_exploratory<.05).sum())))
 pd.DataFrame(summary).to_csv(O/'source_screen_summary.tsv',sep='\t',index=False);out[(out.variant=='abundance_complexity')&(out.joint_BH_exploratory<.05)].to_csv(O/'all_primary_BH05_exploratory_records.tsv',sep='\t',index=False,na_rep='NA');j=dict(status='COMPLETE_EXPLORATORY_CR2_SCREEN_PENDING_NUMERIC_AUDIT_AND_BLOCK_VALIDATION',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),rows=len(out),joint_family_per_variant=59131,summary=summary,physical_binding_or_causal_response_validated=False,approximation_calibration_limit_retained=True,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()});(O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps({k:v for k,v in j.items() if k!='outputs'}))
if __name__=='__main__':main()
