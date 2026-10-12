"""Exploratory ligand-associated expression of unsigned prior targets.

This score is not NicheNet ligand activity, receptor-specific activation,
new formal cell states, causal binding, or independent validation.
All620 fixed interface/node cases and all three variants retain missing cases.
"""
from pathlib import Path
from datetime import datetime, timezone
import os,json,hashlib,gc
import numpy as np
import pandas as pd
import pyreadr
from scipy.stats import t as tdist

PANELS={'IFN_response8':['ISG15','IFIT1','IFIT2','IFIT3','MX1','OAS1','STAT1','IRF7'],
'cytotoxicity7':['PRF1','GZMB','GZMH','GNLY','NKG7','CTSW','GZMA'],
'exhaustion_related6':['PDCD1','LAG3','HAVCR2','TOX','CXCL13','LAYN'],
'cycling8':['MKI67','TOP2A','STMN1','PCNA','TYMS','UBE2C','HMGB2','TUBA1B'],
'naive_memory_related6':['IL7R','CCR7','TCF7','LEF1','LTB','SELL'],
'antigen_presentation8':['HLA-DRA','HLA-DRB1','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQB1','CD74','CIITA']}


def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for v in iter(lambda:f.read(8388608),b''):h.update(v)
 return h.hexdigest()


def standardize(a):
 a=np.asarray(a,float);sd=a.std(axis=0,ddof=1)
 return np.divide(a-a.mean(axis=0),sd,out=np.zeros_like(a),where=sd>1e-12),sd


def bh(p):
 p=np.asarray(p,float);order=np.argsort(p,kind='stable');q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];out=np.empty_like(q);out[order]=np.minimum(q,1);return out


def evaluate(a,y,key,plan):
 n,k=a.shape[0],a.shape[1]+1
 if not np.isfinite(a).all():return dict(status='MISSING_REQUIRED_COVARIATE')
 if n<=k:return dict(status='INSUFFICIENT_RESIDUAL_DF',design_columns=k)
 z,sd=standardize(a);X=np.column_stack([np.ones(n),z]);sy=y.std(ddof=1)
 if not np.isfinite(sy) or sy<=1e-12:return dict(status='TARGET_EXPRESSION_SCORE_CONSTANT')
 Y=(y-y.mean())/sy;s=np.linalg.svd(X,compute_uv=False);rank=int((s>1e-10*s[0]).sum())
 result=dict(design_columns=k,design_rank=rank,constant_covariate_indices=';'.join(map(str,np.flatnonzero(sd<=1e-12))))
 if rank<k:return dict(result,status='REQUIRED_DESIGN_NOT_FULL_RANK')
 cond=float(s[0]/s[-1]);result['design_condition']=cond
 if cond>=1e8:return dict(result,status='REQUIRED_DESIGN_ILL_CONDITIONED')
 inverse=np.linalg.pinv(X,rcond=1e-10);h=np.sum(X*inverse.T,axis=1);result['max_leverage']=float(h.max())
 if h.max()>.95:return dict(result,status='HIGH_INFLUENCE_DESIGN')
 b=inverse@Y;res=Y-X@b;bread=np.linalg.inv(X.T@X);scores=X*res[:,None]
 covariance=bread@(scores.T@scores)@bread*(n/(n-1))*((n-1)/(n-k));se=float(np.sqrt(max(covariance[1,1],0)))
 beta=float(b[1]);crit=tdist.ppf(.975,n-1)
 result.update(status='ESTIMABLE_EXPLORATORY_ASSOCIATION',sender_beta=beta,CR1_se=se,CR1_ci_low=float(beta-crit*se),CR1_ci_high=float(beta+crit*se),analytic_p=float(2*tdist.sf(abs(beta/se),n-1)) if se>0 else None)
 B=plan['resampling']['repeats'];seed=(plan['resampling']['seed']+int(hashlib.sha256(key.encode()).hexdigest()[:8],16))%(2**32);rng=np.random.default_rng(seed)
 reduced=np.delete(X,1,axis=1);fitted=reduced@np.linalg.lstsq(reduced,Y,rcond=None)[0];nullres=Y-fitted
 permutations=np.stack([rng.permutation(n) for _ in range(B)])
 perm_y=fitted[None,:]+nullres[permutations];perm_beta=perm_y@inverse[1]
 result['residual_permutation_p']=float((1+(np.abs(perm_beta)>=abs(beta)-1e-12).sum())/(B+1))
 draws=rng.integers(0,n,(B,n));valid_b=[]
 for lo in range(0,B,100):
  pick=draws[lo:lo+100];xb=X[pick];yb=Y[pick];sing=np.linalg.svd(xb,compute_uv=False);valid=(sing[:,-1]>1e-10*sing[:,0]) & (sing[:,0]/np.maximum(sing[:,-1],1e-300)<1e8)
  if valid.any():
   xv=xb[valid];iv=np.linalg.pinv(xv,rcond=1e-10);lv=np.sum(xv*np.swapaxes(iv,1,2),axis=2).max(axis=1);bb=np.einsum('bij,bj->bi',iv,yb[valid])[:,1];valid_b.extend(bb[lv<=.95].tolist())
 result.update(resampling_seed=int(seed),bootstrap_valid=len(valid_b),bootstrap_requested=B,score_and_scaling_frozen_during_resampling=True)
 if len(valid_b)>=plan['resampling']['minimum_valid_bootstrap']:
  low,high=np.quantile(valid_b,[.025,.975]);result.update(bootstrap_ci_low=float(low),bootstrap_ci_high=float(high),bootstrap_interval_status='AVAILABLE')
 else:result.update(bootstrap_ci_low=None,bootstrap_ci_high=None,bootstrap_interval_status='INSUFFICIENT_VALID_RESAMPLED_DESIGNS')
 return result,X,Y


def main():
 started=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';O=P/'fullgene_receiver_target_association_v1';F=P/'fullgene_receiver_target_readiness_v1';N=P/'supported_LR_native_cell_panel_v1';C=P/'source_reviewed_family_LR_cache_v1';Q=P/'fullgene_receiver_response_preflight_v1';U=P/'patient_block_provenance_overlay_v1';planpath=R/'fullgene_receiver_target_code_v1/fullgene_receiver_target_association_plan_v1.json';plan=json.loads(planpath.read_text());assert not O.exists()
 assert sha(__file__)==plan['code_sha256'][Path(__file__).name]
 paths={'preflight_inventory':Q/'inventory.json','preflight_audit':Q/'independent_audit.json','family_cache_inventory':C/'inventory.json','family_cache_audit':C/'independent_audit.json','target_prior':F/'verified_prior_inputs_v2/ligand_target_matrix_nsga2r_final.rds'}
 for k,p in paths.items():assert sha(p)==plan['input_sha256'][k],k
 fi=json.loads((F/'inventory.json').read_text());assert sha(F/'independent_audit.json')==plan['fullgene_audit_sha256']
 for n,h in fi['outputs'].items():assert sha(F/n)==h
 ci=json.loads((C/'inventory.json').read_text())
 for p,h in ci['outputs'].items():assert sha(p)==h
 ui=json.loads((U/'inventory.json').read_text());scopepath=U/'private_all1290_patient_provenance_overlay.tsv';assert sha(scopepath)==ui['outputs'][scopepath.name];scope=pd.read_csv(scopepath,sep='\t',keep_default_na=False).set_index('sampleID')
 prior=next(iter(pyreadr.read_r(str(paths['target_prior'])).values()));prior.index=prior.index.astype(str);prior.columns=prior.columns.astype(str)
 ptable=pd.read_csv(F/'all336_prior_coverage.tsv',sep='\t',keep_default_na=False);cases=pd.read_csv(Q/'all_interface_original_node_preflight.tsv',sep='\t',keep_default_na=False);assert len(cases)==620
 with np.load(C/'source_reviewed_native_family_expression_cache.npz',allow_pickle=False) as c:
  cm=int(np.flatnonzero(c['modes']=='marker_supported_main')[0]);cf=list(c['families']);cs=pd.Index(c['samples']);cn=c['n_cells'][cm];frac=c['within_major_fractions'][cm]
 O.mkdir();(O/'private_designs').mkdir();rows=[];score_details=[];tested=0
 records=json.loads((N/'inventory.json').read_text())['source_records']
 for rec in records:
  sid=rec['gse_id'];b=dict(np.load(F/'private_source_unit_arrays'/(sid+'.npz'),allow_pickle=False));samples=b['samples'];ns=len(samples);nodes=list(b['nodes']);gi={g:i for i,g in enumerate(b['genes'])};count=b['n_cells'].reshape(len(nodes),ns);cp=b['mean_CP10K'].reshape(len(nodes),ns,len(gi));lg=b['mean_log1p_CP10K'].reshape(len(nodes),ns,len(gi));rna=b['full_RNA_total_sum'].reshape(len(nodes),ns);ng=b['full_detected_gene_sum'].reshape(len(nodes),ns);sc=scope.loc[samples];ix=cs.get_indexer(samples);assert (ix>=0).all()
  for fi2,fam in enumerate(cf):assert np.array_equal(count[np.asarray(b['node_family'])==fam].sum(axis=0),cn[fi2,ix])
  for idx,row in cases[cases.gse_id.eq(sid)].iterrows():
   test_id='NODE%04d'%(idx+1);base={k:row[k] for k in ['interface_id','gse_id','stratum','sender_family','receiver_family','ligand','receptor','receiver_original_node','LR_status']};base.update(test_id=test_id,scope_selected_from_current_data=True,response_is_unsigned_prior_target_expression_not_ligand_activity=True)
   if row.status!='COUNT_VARIATION_DIAGNOSTIC_DESIGN_SCREEN_PASSED_NOT_INFERENCE':
    for variant in ['M0_expression_depth','M1_abundance_identity','M2_sender_shared_programmes']:rows.append(dict(base,variant=variant,status='PREFLIGHT_'+row.status))
    continue
   pr=ptable[(ptable.gse_id.eq(sid)) & ptable.sender_family.eq(row.sender_family) & ptable.receiver_family.eq(row.receiver_family) & ptable.ligand.eq(row.ligand) & ptable.receptor.eq(row.receptor)].iloc[0];targets=pr.top200_after_exclusions.split(';');native_targets=[g for g in targets if g in gi];weights=prior.loc[native_targets,row.ligand].to_numpy(float);total=prior.loc[targets,row.ligand].sum();weight_coverage=float(weights.sum()/total);weights/=weights.sum()
   si=np.flatnonzero(b['node_family']==row.sender_family);ri=nodes.index(row.receiver_original_node);sn=count[si].sum(axis=0);rn=count[ri];sr=rna[si].sum(axis=0);sg=ng[si].sum(axis=0)
   sender_cp=np.divide(np.nansum(cp[si]*count[si,:,None],axis=0),sn[:,None],out=np.full((ns,len(gi)),np.nan),where=sn[:,None]>0);sender_lg=np.divide(np.nansum(lg[si]*count[si,:,None],axis=0),sn[:,None],out=np.full((ns,len(gi)),np.nan),where=sn[:,None]>0)
   context=sc.stratum.eq(row.stratum).to_numpy() & sc.analysis_eligible.to_numpy(dtype=bool);ok=context & (sn>=20) & (rn>=20);ids=np.flatnonzero(ok);meta=sc.iloc[ids];assert len(ids)==int(row.n_libraries_endpoint_cells_ge20)
   if meta.effective_resampling_block.nunique()!=len(meta):
    for variant in ['M0_expression_depth','M1_abundance_identity','M2_sender_shared_programmes']:rows.append(dict(base,variant=variant,status='MULTI_LIBRARY_BLOCK_PERMUTATION_NOT_PRESPECIFIED'))
    continue
   expression=lg[ri,ok][:,[gi[g] for g in native_targets]];z,sd=standardize(expression);Y=z@weights
   base.update(n_libraries=len(ids),n_patient_blocks=meta.effective_resampling_block.nunique(),native_target_genes=len(native_targets),variable_target_genes=int((sd>1e-12).sum()),native_prior_weight_coverage=weight_coverage)
   Lgenes=row.ligand.split('_');Rgenes=row.receptor.split('_');lv=np.log1p(sender_cp[ok][:,[gi[g] for g in Lgenes]]).mean(axis=1);rv=np.log1p(cp[ri,ok][:,[gi[g] for g in Rgenes]]).mean(axis=1)
   a0=np.column_stack([lv,rv,np.log1p(sn[ok]),np.log1p(rn[ok]),np.log1p(sr[ok]/sn[ok]),np.log1p(rna[ri,ok]/rn[ok]),np.log1p(sg[ok]/sn[ok]),np.log1p(ng[ri,ok]/rn[ok])]);names0=['L','R','log_sender_n','log_receiver_n','log_sender_RNA','log_receiver_RNA','log_sender_nGenes','log_receiver_nGenes']
   sf=cf.index(row.sender_family);rf=cf.index(row.receiver_family);rfn=cn[rf,ix];rfrac=frac[rf,ix][ok]*rn[ok]/rfn[ok];assert (rfrac<=1+1e-9).all()
   present=[j for j in si if count[j,context].sum()>0];reference=present[-1];other=[j for j in present if j!=reference];mix=np.column_stack([count[j,ok]/sn[ok] for j in other]) if other else np.empty((len(ids),0))
   a1=np.column_stack([a0,frac[sf,ix][ok],rfrac,mix]);names1=names0+['sender_within_major_fraction','receiver_node_within_major_fraction']+['sender_mix_'+nodes[j] for j in other]
   exclude=set(Lgenes+Rgenes+['CCL5','IL7R','MS4A1','FCER1G']);programs=[];program_fields={}
   applicable=list(PANELS) if row.sender_family.startswith('T_cell_') else ['IFN_response8','cycling8','antigen_presentation8']
   for label in applicable:
    panel=PANELS[label]
    remain=[g for g in panel if g not in exclude];measured=[g for g in remain if g in gi];valid=len(measured)>=4 and len(measured)>=.6*len(remain)
    programs.append(sender_lg[ok][:,[gi[g] for g in measured]].mean(axis=1) if valid else np.full(len(ids),np.nan));program_fields[label]=';'.join(measured)
   a2=np.column_stack([a1,np.column_stack(programs)]);names2=names1+['sender_program_'+p for p in applicable]
   detail=dict(test_id=test_id,interface_id=row.interface_id,receiver_original_node=row.receiver_original_node,target_genes=';'.join(native_targets),sender_mixture_reference=nodes[reference],sender_mixture_members=';'.join(nodes[j] for j in other),sender_programme_retained_genes=program_fields,weights_unsigned=True,weight_coverage=weight_coverage)
   score_details.append(detail);design=dict(raw_target_expression_score=Y,sample_keys=samples[ok],patient_blocks=meta.effective_resampling_block.to_numpy(dtype=str),target_genes=np.asarray(native_targets),target_prior_weights=weights,target_means=expression.mean(axis=0),target_sds=sd,raw_target_gene_expression=expression)
   for variant,a,names in [('M0_expression_depth',a0,names0),('M1_abundance_identity',a1,names1),('M2_sender_shared_programmes',a2,names2)]:
    design['raw_a_'+variant]=a;design['names_'+variant]=np.asarray(['intercept']+names)
    evaluated=evaluate(a,Y,test_id+'|'+variant,plan)
    if isinstance(evaluated,tuple):result,X,y=evaluated;design['X_'+variant]=X;design['y_'+variant]=y;tested+=1
    else:result=evaluated
    result.update(base,variant=variant,covariate_names=';'.join(['intercept']+names),sender_mixture_reference=nodes[reference]);rows.append(result)
   np.savez_compressed(O/'private_designs'/(test_id+'.npz'),**design)
  del b,cp,lg,sender_cp,sender_lg;gc.collect()
 result=pd.DataFrame(rows);assert len(result)==1860 and not result[['test_id','variant']].duplicated().any()
 result['diagnostic_BH_q']=np.nan
 for variant in result.variant.unique():
  mask=result.variant.eq(variant);p=pd.to_numeric(result.loc[mask,'residual_permutation_p'],errors='coerce');adjusted=bh(p.fillna(1).to_numpy());result.loc[mask,'diagnostic_BH_q']=np.where(p.notna(),adjusted,np.nan)
 result.to_csv(O/'all1860_target_association_models.tsv',sep='\t',index=False,na_rep='NA');(O/'score_provenance.json').write_text(json.dumps(score_details,indent=2)+'\n')
 summary=dict(frozen_node_cases=620,model_records=1860,estimated_models=tested,status_by_variant={v:g.status.value_counts().to_dict() for v,g in result.groupby('variant')},diagnostic_q05_by_variant={v:int(g.diagnostic_BH_q.lt(.05).sum()) for v,g in result.groupby('variant')},bootstrap_CI_available_by_variant={v:int(g.get('bootstrap_interval_status',pd.Series(dtype=str)).eq('AVAILABLE').sum()) for v,g in result.groupby('variant')},confirmatory_or_independent_validation=False)
 inventory=dict(status='COMPLETE_TARGET_ASSOCIATION_PENDING_INDEPENDENT_AUDIT',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planpath),summary=summary,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()},new_exploratory_models=True,formal_states=False,canonical_Zarr_changed=False,scientific_release=False,axis_validation_complete=False)
 (O/'inventory.json').write_text(json.dumps(inventory,indent=2)+'\n');print(json.dumps(inventory['summary']),flush=True)


if __name__=='__main__':main()
