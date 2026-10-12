"""All-effect audit using complete X pseudoinverse and direct block residual CR2.
No producer/core import. No individual rows exported.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
from scipy.stats import t as td
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def z(a):
 a=np.asarray(a,float);m=a.mean(0);s=np.sqrt(np.mean((a-m)**2,axis=0));return (a-m)/np.where(s>1e-12,s,1)
def adjust(p):
 p=np.asarray(p,float);i=np.argsort(p,kind='stable');q=p[i]*len(p)/np.arange(1,len(p)+1);a=np.empty(len(p));a[i]=np.minimum(1,np.minimum.accumulate(q[::-1])[::-1]);return a

def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);P=R/'immune_same_name_family_source_review_v1';N=P/'major_pair_native_assembly_v1';T=P/'major_receiver_TF_preflight_v3';D=P/'major_LR_TF_design_preflight_v1';O=P/'major_LR_receiver_TF_screen_v1';plan=json.loads(Path(os.environ['PANCANCER_LR_TF_SCREEN_PLAN']).read_text());assert sha(__file__)==plan['auditor_sha256'] and not (O/'independent_numeric_audit.json').exists();inv=json.loads((O/'inventory.json').read_text())
 for n,h in inv['outputs'].items():assert sha(O/n)==h
 table=pd.read_csv(O/'all177393_LR_union_models.tsv.gz',sep='\t',na_values=['NA'],keep_default_na=False);assert len(table)==177393 and not table.duplicated(['context_id','LR_index','variant']).any();counts=pd.read_csv(AR/'20261010_expression_annotation_downstream_v1/inputs/sample_major_subcluster_counts.tsv',sep='\t');total=counts.groupby('sampleID').n_cells.sum();contexts=pd.read_csv(D/'all_context_model_readiness.tsv',sep='\t',keep_default_na=False);all_values=0;all_LR=0;maxdiff=0.
 for row in contexts.to_dict('records'):
  cid=row['context_id'];a=np.load(D/'private_context_arrays'/(cid+'.npz'));b=np.load(T/'private_context_arrays'/(cid+'.npz'));e=np.load(O/'private_component_effects'/(cid+'.npz'));nat=np.load(N/'private_group_arrays'/(row['group_id']+'.npz'));mask=nat['corrected_material_for_next_inference']==row['material'];ns=nat['source_n'][mask];nt=nat['target_n'][mask];samples=a['samples'];n=len(samples);base=np.column_stack((ns/total.loc[samples],nt/total.loc[samples],np.log1p(ns),np.log1p(nt),nat['source_mean_detected_nonLR_genes'][mask],nat['target_mean_detected_nonLR_genes'][mask]));blocks=b['blocks'];assert np.array_equal(e['candidate_LR_index'],a['candidate_LR_index']) and np.array_equal(e['candidate_TF_index'],a['candidate_TF_index']);ct=table[table.context_id==cid].set_index(['LR_index','variant'])
  for li in np.unique(a['candidate_LR_index']):
   where=np.flatnonzero(a['candidate_LR_index']==li);tfidx=a['candidate_TF_index'][where];K=len(where);Z=np.column_stack((np.ones(n),z(np.column_stack((base,nat['receptor_meanlog'][mask,li])))));Zm=np.column_stack((Z,z(a['source_target_mixture_PCs'])));x=z(nat['ligand_meanlog'][mask,li])
   for variant in ['abundance_complexity','subtype_mixture','full_gene_sensitivity']:
    all_LR+=1;reported=ct.loc[(int(li),variant)];dstatus=int(a['mixture_design_status' if variant=='subtype_mixture' else 'primary_design_status'][li]);q=np.column_stack((Zm if variant=='subtype_mixture' else Z,x))
    if dstatus!=0:
     assert reported.status=='design_'+a['status_names'][dstatus] and not np.isfinite(float(reported.joint_p_exploratory));continue
    inverse=np.linalg.pinv(q,rcond=1e-10);c=inverse[-1];H=q@inverse;H=(H+H.T)/2;M=np.eye(n)-H;Y=z(b['sensitivity_TF_scores' if variant=='full_gene_sensitivity' else 'primary_TF_scores'][:,tfidx]);beta=(inverse@Y)[-1];res=Y-q@(inverse@Y);scores=[];Acols=[];singular=False;maxhat=0
    for block in np.unique(blocks):
     ix=np.flatnonzero(blocks==block);vals,V=np.linalg.eigh(np.eye(len(ix))-H[np.ix_(ix,ix)]);maxhat=max(maxhat,float(1-vals.min()))
     if vals.min()<=1e-8:singular=True;break
     Ag=(V/np.sqrt(vals))@V.T;scores.append(c[ix]@Ag@res[ix]);vec=np.zeros(n);vec[ix]=Ag@c[ix];Acols.append(M@vec)
    if singular:
     assert reported.status=='cluster_adjustment_singular';continue
    B=np.column_stack(Acols);E=B.T@B;df=float(np.trace(E)**2/np.sum(E*E));se=np.sqrt(np.sum(np.array(scores)**2,axis=0));p=2*td.sf(np.divide(abs(beta),se,out=np.full_like(beta,np.nan),where=se>1e-14),df);crit=td.ppf(.975,df);expect={'beta':beta,'se':se,'component_p_diagnostic':p,'lower_marginal':beta-crit*se,'upper_marginal':beta+crit*se,'df':np.full(K,df)}
    for k,value in expect.items():
     actual=e[variant+'__'+k][where];assert np.allclose(actual,value,rtol=3e-6,atol=3e-7,equal_nan=True),(cid,int(li),variant,k);finite=np.isfinite(actual)&np.isfinite(value)
     if finite.any():maxdiff=max(maxdiff,float(np.max(abs(actual[finite]-value[finite]))));all_values+=int(finite.sum())
    assert np.isclose(float(reported.effective_df),df,rtol=3e-6,atol=3e-7)
    if df<plan['min_effective_df']:assert reported.status=='effective_df_insufficient' and not np.isfinite(float(reported.joint_p_exploratory))
    elif not np.isfinite(p).any():assert reported.status=='outcome_residual_variance_insufficient'
    else:
     assert reported.status=='inference_estimable';pp=np.where(np.isfinite(p),p,1);joint=min(1,K*pp.min());assert np.isclose(float(reported.joint_p_exploratory),joint,rtol=3e-6,atol=3e-7);lead=int(np.argmin(pp));assert reported.lead_TF==str(b['TFs'][tfidx[lead]]) and np.isclose(float(reported.lead_beta),beta[lead],rtol=3e-6,atol=3e-7)
  print(json.dumps(dict(context=cid,LR_variants_checked=all_LR,component_values_checked=all_values)),flush=True)
 for v,g in table.groupby('variant'):
  assert len(g)==59131;p=pd.to_numeric(g.joint_p_exploratory,errors='coerce').fillna(1).to_numpy();q=adjust(p);good=g.status.eq('inference_estimable').to_numpy();actual=pd.to_numeric(g.joint_BH_exploratory,errors='coerce').to_numpy();assert np.allclose(actual[good],q[good],atol=1e-12) and np.isnan(actual[~good]).all();by=np.minimum(1,q*sum(1/i for i in range(1,len(g)+1)));actual=pd.to_numeric(g.joint_BY_exploratory,errors='coerce').to_numpy();assert np.allclose(actual[good],by[good],atol=1e-12)
 assert all_LR==177393;result=dict(status='PASS_ALL_FULLX_CR2_EFFECTS_DF_LR_UNION_AND_MULTIPLICITY_EXPLORATORY_ONLY',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),all_LR_variants_checked=all_LR,all_finite_component_values_checked=all_values,max_absolute_numeric_difference=maxdiff,approximation_and_confirmation_limits_retained=True,no_independent_validation_or_causal_claim=True);(O/'independent_numeric_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
