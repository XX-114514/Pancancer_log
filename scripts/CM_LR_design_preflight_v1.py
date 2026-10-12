"""Full eligible CM/LR design diagnostics, no coefficients or p-values.

Strict original design retained. Nuisance-span identifiability is diagnostic
only; no required covariate is removed or inference rescued in this run.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd

NAMES=['intercept','CM_activity','source_within_major','target_within_major','log_source_cells','log_target_cells']
VARIANTS=['endpoint_excluded_library','original_H_library','endpoint_excluded_block_mean_diagnostic']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def standardize(a):
    sd=np.std(a,axis=0,ddof=1);constant=sd<=1e-12
    x=np.divide(a-a.mean(axis=0),sd,out=np.zeros_like(a,dtype=float),where=~constant)
    return np.column_stack([np.ones(len(a)),x]),constant

def diagnostics(a,blocks):
    n=len(a);X,const=standardize(a);U,s,V=np.linalg.svd(X,full_matrices=False);tol=np.finfo(float).eps*max(X.shape)*s[0];rank=int((s>tol).sum());Z=X[:,[0,2,3,4,5]];Uz,sz,_=np.linalg.svd(Z,full_matrices=False);ztol=np.finfo(float).eps*max(Z.shape)*sz[0];zr=int((sz>ztol).sum());v=X[:,1]-Uz[:,:zr]@(Uz[:,:zr].T@X[:,1]);residual=float(v@v);condition=float(s[0]/s[-1]) if rank==6 else None;leverage=float(np.max(np.sum(U[:,:rank]**2,axis=1)));unique=not const[0] and residual>1e-10*max(1.,float(X[:,1]@X[:,1]));G=len(set(blocks))
    if n<20 or G<10:status='INSUFFICIENT20_ROWS_OR10_BLOCKS'
    elif const[0]:status='CONSTANT_CM_ACTIVITY'
    elif const[1:].any():status='CONSTANT_REQUIRED_NUISANCE_COLUMN'
    elif rank<6:status='COLLINEAR_REQUIRED_DESIGN'
    elif condition>=1e8:status='ILL_CONDITIONED_REQUIRED_DESIGN'
    elif leverage>.95:status='HIGH_INFLUENCE_DESIGN'
    else:status='STRICT_DESIGN_SCREEN_PASSED_NOT_INFERENCE'
    return dict(n_rows=n,n_analysis_blocks=G,n_repeated_library_rows=n-G,design_columns=6,design_rank=rank,nuisance_span_rank=zr,constant_columns=';'.join(np.array(NAMES[1:])[const]),design_condition=condition,max_leverage=leverage,CM_residual_sum_squares_given_all_nuisance_columns=residual,CM_coefficient_identifiable_in_full_nuisance_span=unique,status=status)

def main():
    start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';L=P/'CM_edge_LR_extraction_v2';O=P/'CM_LR_model_design_preflight_v1';planpath=Path(os.environ['PANCANCER_CM_DESIGN_PLAN']);plan=json.loads(planpath.read_text());assert not O.exists() and sha(__file__)==plan['code_sha256'][Path(__file__).name]
    inv=json.loads((L/'inventory.json').read_text());assert sha(L/'inventory.json')==plan['inventory_sha256'] and sha(L/'independent_audit.json')==plan['audit_sha256'] and sha(L/'all314432_consensus_case_coverage.tsv.gz')==plan['case_table_sha256'];assert json.loads((L/'independent_audit.json').read_text())['status'].startswith('PASS_ALL314432')
    for n,h in inv['outputs'].items():assert sha(L/n)==h
    cases=pd.read_csv(L/'all314432_consensus_case_coverage.tsv.gz',sep='\t',keep_default_na=False);cases=cases[cases.data_coverage_gate_only.eq(True)];assert len(cases)==5377;rows=[];cache={};mapped=0
    for gid,g in cases.groupby('group_id',sort=True):
        b=dict(np.load(L/'private_group_arrays'/(gid+'.npz'),allow_pickle=False));lrindex={(l,r):i for i,(l,r) in enumerate(zip(b['ligands'],b['receptors']))};cmindex={c:i for i,c in enumerate(b['CMs'])}
        for case in g.itertuples():
            j=lrindex[case.ligand,case.receptor];obs=np.isfinite(b['scores'][:,j,0]);Y=np.log1p(b['scores'][obs,j,0]);blocks=b['patient_blocks'][obs];review=b['reviewed_patient'][obs];assert obs.sum()==case.n_observed_libraries and len(set(blocks))==case.n_observed_blocks;controls=np.column_stack([b['source_within_major'][obs],b['target_within_major'][obs],np.log1p(b['source_n'][obs]),np.log1p(b['target_n'][obs])]);maskkey=np.packbits(obs).tobytes()
            for cm in case.CMs.split(';'):
                mapped+=1
                for variant in VARIANTS:
                    key=(gid,maskkey,cm,variant);activity=(b['original_H'] if variant=='original_H_library' else b['endpoint_excluded_activity'])[obs,cmindex[cm]];raw=np.column_stack([activity,controls]);yy=Y;bb=blocks;rr=review
                    if variant=='endpoint_excluded_block_mean_diagnostic':
                        order=list(dict.fromkeys(blocks));raw=np.stack([raw[blocks==k].mean(axis=0) for k in order]);yy=np.array([Y[blocks==k].mean() for k in order]);rr=np.array([review[blocks==k].all() for k in order]);bb=np.array(order)
                    if key not in cache:
                        cache[key]=diagnostics(raw,bb) if np.isfinite(raw).all() else dict(status='CM_ACTIVITY_UNAVAILABLE',n_rows=len(raw),n_analysis_blocks=len(set(bb)))
                    result=cache[key].copy();result['response_variable']=bool(np.std(yy,ddof=1)>1e-12)
                    if not result['response_variable']:result['status']='NO_RESPONSE_VARIATION'
                    rows.append({k:getattr(case,k) for k in ['case_id','group_id','pair_id','stratum','gse_id','cancerType','source','target','ligand','receptor']}|dict(cm=cm,variant=variant,n_reviewed_rows=int(rr.sum()),n_reviewed_source_patient_blocks=len(set(bb[rr])),coverage_only_source_scope_selected_from_current_data=True)|result)
    assert mapped==13497 and len(rows)==40491
    O.mkdir();d=pd.DataFrame(rows);d.to_csv(O/'all40491_model_design_diagnostics.tsv.gz',sep='\t',index=False,na_rep='NA');summary={v:dict(records=len(g),status_counts=g.status.value_counts().to_dict(),partial_CM_identifiable_diagnostic=int(g.CM_coefficient_identifiable_in_full_nuisance_span.eq(True).sum())) for v,g in d.groupby('variant')}
    out=dict(status='COMPLETE_DESIGN_PREFLIGHT_PENDING_INDEPENDENT_AUDIT',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planpath),extraction_inventory_sha256=sha(L/'inventory.json'),coverage_cases=5377,mapped_CM_cases=13497,diagnostic_records=40491,distinct_X_templates=len(cache),summary=summary,primary_coefficient_estimated=False,new_p_values=False,required_covariates_removed=False,partial_span_not_an_inference_acceptance=True,outputs={'all40491_model_design_diagnostics.tsv.gz':sha(O/'all40491_model_design_diagnostics.tsv.gz')})
    (O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
