"""Full nuisance-span CM coefficient, strict original design retained in output.

Exploratory association, not causal communication/CM-driver validation.
All required controls remain in their column space; no p-based variable choice.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os
import numpy as np,pandas as pd
from scipy.stats import t
from statsmodels.stats.multitest import multipletests

VARIANTS=['endpoint_excluded_library','original_H_library']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def seed(key):return int(hashlib.sha256(repr(key).encode()).hexdigest()[:8],16)

def setup(raw,blocks,key,draws=2000):
    n=len(raw);means=raw.mean(axis=0);sd=raw.std(axis=0,ddof=1);const=sd<=1e-12;z=np.divide(raw-means,sd,out=np.zeros_like(raw),where=~const);Z=np.column_stack([np.ones(n),z[:,1:]]);u,s,_=np.linalg.svd(Z,full_matrices=False);rz=int((s>np.finfo(float).eps*max(Z.shape)*s[0]).sum());Q=u[:,:rz];x=z[:,0];v=x-Q@(Q.T@x);ss=v@v
    result=dict(n_rows=n,n_analysis_blocks=len(set(blocks)),nuisance_rank=rz,constant_nuisance_indices=';'.join(str(i) for i in np.flatnonzero(const[1:])),resampling_seed=seed(key))
    if const[0] or ss<=1e-10*max(1.,x@x):return result|dict(status='CM_EXPOSURE_NOT_IDENTIFIABLE')
    B=np.column_stack([Q,x]);sv=np.linalg.svd(B,compute_uv=False);condition=float(sv[0]/sv[-1]);rank=rz+1;inv=np.linalg.inv(B.T@B);hat=np.einsum('ni,ij,nj->n',B,inv,B);lev=float(hat.max());result.update(design_rank=rank,effective_basis_condition=condition,max_leverage=lev)
    if n<20 or len(set(blocks))<10:return result|dict(status='INSUFFICIENT_ROWS_OR_BLOCKS')
    if condition>=1e8 or lev>.95:return result|dict(status='ILL_CONDITIONED_OR_HIGH_INFLUENCE_DESIGN')
    order=list(dict.fromkeys(blocks));indices=[np.flatnonzero(blocks==k) for k in order];G=len(order);rng=np.random.default_rng(result['resampling_seed']);bins={}
    for i,ids in enumerate(indices):bins.setdefault(len(ids),[]).append(i)
    permutations=np.tile(np.arange(n),(draws,1))
    for d in range(draws):
        for ids in bins.values():
            perm=rng.permutation(ids)
            for target,source in zip(ids,perm):permutations[d,indices[target]]=indices[source]
    permutation_linear=np.empty((draws,n));a=v/ss
    for d in range(draws):permutation_linear[d,permutations[d]]=a
    blockdraw=rng.integers(0,G,size=(draws,G));blockweights=np.zeros((draws,G),dtype=np.int16)
    for d in range(draws):blockweights[d]=np.bincount(blockdraw[d],minlength=G)
    rowblock=np.array([order.index(k) for k in blocks]);weights=blockweights[:,rowblock].astype(float);gram=np.einsum('dn,ni,nj->dij',weights,B,B);eig=np.linalg.eigvalsh(gram);valid=(eig[:,0]>1e-12*eig[:,-1])&(eig[:,-1]/np.maximum(eig[:,0],1e-300)<1e16);ids=np.flatnonzero(valid);bootlinear=[]
    if len(ids):
        gi=np.linalg.inv(gram[ids]);h=np.einsum('ni,dij,nj->dn',B,gi,B);h=np.where(weights[ids]>0,h,0);ok=h.max(axis=1)<=.95;ids=ids[ok];gi=gi[ok];bootlinear=np.einsum('di,ni,dn->dn',gi[:,-1,:],B,weights[ids])
    result.update(status='ESTIMABLE_FULL_CONTROL_SPAN_ASSOCIATION',boot_valid=len(ids),permutation_size_bins=';'.join(str(k)+':'+str(len(vv)) for k,vv in sorted(bins.items())),fixed_singleton_size_bins=sum(len(vv)==1 for vv in bins.values()))
    result['_internal']=(Q,B,inv,a,indices,permutation_linear,np.asarray(bootlinear),rank)
    return result

def evaluate(template,y):
    r={k:v for k,v in template.items() if k!='_internal'}
    if template['status']!='ESTIMABLE_FULL_CONTROL_SPAN_ASSOCIATION':return r
    if not np.isfinite(y).all() or y.std(ddof=1)<=1e-12:return r|dict(status='NO_RESPONSE_VARIATION')
    yy=(y-y.mean())/y.std(ddof=1);Q,B,inv,a,blocks,pl,bl,rank=template['_internal'];beta=float(a@yy);fit=B@(inv@(B.T@yy));err=yy-fit;n=len(yy);G=len(blocks);variance=sum(float(a[ids]@err[ids])**2 for ids in blocks)*G/(G-1)*(n-1)/(n-rank);se=float(np.sqrt(max(0,variance)))
    if se<=1e-12:return r|dict(status='DEGENERATE_CLUSTER_VARIANCE')
    residual=yy-Q@(Q.T@yy);permbeta=pl@residual;p=float((1+np.count_nonzero(abs(permbeta)>=abs(beta)-1e-12))/(len(permbeta)+1));ci=t.ppf(.975,G-1)*se;bootstrap=bl@yy if len(bl) else np.array([])
    r.update(beta=beta,CR1_se=se,CR1_ci_low=beta-ci,CR1_ci_high=beta+ci,analytic_p=float(2*t.sf(abs(beta/se),G-1)),residual_block_permutation_p=p,boot_requested=2000,bootstrap_ci_status='AVAILABLE' if len(bootstrap)>=1800 else 'INSUFFICIENT_VALID_RESAMPLED_DESIGNS',bootstrap_ci_low=float(np.quantile(bootstrap,.025)) if len(bootstrap)>=1800 else None,bootstrap_ci_high=float(np.quantile(bootstrap,.975)) if len(bootstrap)>=1800 else None)
    return r

def main():
    start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';L=P/'CM_edge_LR_extraction_v2';D=P/'CM_LR_model_design_preflight_v1';O=P/'CM_LR_full_control_association_v1';planpath=Path(os.environ['PANCANCER_CM_ASSOCIATION_PLAN']);plan=json.loads(planpath.read_text());assert not O.exists() and sha(__file__)==plan['code_sha256'][Path(__file__).name]
    for k,p in [('extraction_inventory',L/'inventory.json'),('extraction_audit',L/'independent_audit.json'),('design_inventory',D/'inventory.json'),('design_audit',D/'independent_audit.json')]:assert sha(p)==plan['input_sha256'][k]
    li=json.loads((L/'inventory.json').read_text());cases=pd.read_csv(L/'all314432_consensus_case_coverage.tsv.gz',sep='\t',keep_default_na=False);assert sha(L/'all314432_consensus_case_coverage.tsv.gz')==li['outputs']['all314432_consensus_case_coverage.tsv.gz'];cases=cases[cases.data_coverage_gate_only.eq(True)];assert len(cases)==5377;diagnostics=pd.read_csv(D/'all40491_model_design_diagnostics.tsv.gz',sep='\t').set_index(['case_id','cm','variant']);rows=[];O.mkdir();(O/'private_resampling_examples').mkdir();saved=0
    for gid,g in cases.groupby('group_id',sort=True):
        file=L/'private_group_arrays'/(gid+'.npz');assert sha(file)==li['outputs'][str(file.relative_to(L))];b=dict(np.load(file,allow_pickle=False));cache={};index={(l,r):i for i,(l,r) in enumerate(zip(b['ligands'],b['receptors']))};cmidx={c:i for i,c in enumerate(b['CMs'])}
        for case in g.itertuples():
            j=index[case.ligand,case.receptor];mask=np.isfinite(b['scores'][:,j,0]);Y=np.log1p(b['scores'][mask,j,0]);blocks=b['patient_blocks'][mask];controls=np.column_stack([b['source_within_major'][mask],b['target_within_major'][mask],np.log1p(b['source_n'][mask]),np.log1p(b['target_n'][mask])]);review=b['reviewed_patient'][mask]
            for cm in case.CMs.split(';'):
                for variant in VARIANTS:
                    key=(gid,np.packbits(mask).tobytes(),cm,variant);activity=b['endpoint_excluded_activity' if variant=='endpoint_excluded_library' else 'original_H'][mask,cmidx[cm]];raw=np.column_stack([activity,controls])
                    if key not in cache:cache[key]=setup(raw,blocks,key)
                    result=evaluate(cache[key],Y);base={k:getattr(case,k) for k in ['case_id','group_id','pair_id','stratum','gse_id','cancerType','source','target','ligand','receptor','binding_status']};base.update(cm=cm,variant=variant,original_strict_design_status=diagnostics.loc[(case.case_id,cm,variant),'status'],n_reviewed_source_patient_blocks=len(set(blocks[review])),scope_selected_from_current_data=True,full_required_nuisance_span_retained=True,independent_validation=False);rows.append(base|result)
                    if saved<6 and result['status']=='ESTIMABLE_FULL_CONTROL_SPAN_ASSOCIATION':
                        np.savez_compressed(O/'private_resampling_examples'/('example'+str(saved+1)+'.npz'),raw_predictors=raw,raw_response=Y,patient_blocks=blocks,samples=b['samples'][mask],variant=variant,case_id=case.case_id,cm=cm,resampling_seed=result['resampling_seed']);saved+=1
        print('finished_group',gid,'models',len(rows),flush=True)
    assert len(rows)==26994;d=pd.DataFrame(rows)
    for variant,g in d.groupby('variant'):
        p=g.residual_block_permutation_p.fillna(1).to_numpy();q=multipletests(p,method='fdr_bh')[1];d.loc[g.index,'full_family_BH_q']=np.where(g.residual_block_permutation_p.notna(),q,np.nan)
    d.to_csv(O/'all26994_full_control_CM_LR_associations.tsv.gz',sep='\t',index=False,na_rep='NA');summary={v:dict(models=len(g),status_counts=g.status.value_counts().to_dict(),diagnostic_q05=int(g.full_family_BH_q.lt(.05).sum()),bootstrap_CI_available=int(g.bootstrap_ci_status.eq('AVAILABLE').sum())) for v,g in d.groupby('variant')}
    inv=dict(status='COMPLETE_EXPLORATORY_ASSOCIATION_PENDING_INDEPENDENT_AUDIT',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planpath),summary=summary,new_coefficients=True,primary_strict_statuses_preserved=True,all_required_controls_in_span=True,scientific_or_causal_validation=False,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()})
    (O/'inventory.json').write_text(json.dumps(inv,indent=2)+'\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
