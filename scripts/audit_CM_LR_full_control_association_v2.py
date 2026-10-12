"""Independent full-X Moore-Penrose coefficient/CR1/BH audit, no producer import."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os
import numpy as np,pandas as pd
from scipy.stats import t
from scipy.linalg import qr
from statsmodels.stats.multitest import multipletests

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def exact_int(value):
    v=float(value);assert np.isfinite(v) and v.is_integer();return int(v)

def main():
    start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';L=P/'CM_edge_LR_extraction_v2';O=P/'CM_LR_full_control_association_v1';planpath=Path(os.environ['PANCANCER_CM_ASSOCIATION_PLAN']);plan=json.loads(planpath.read_text());repair=json.loads(Path(os.environ['PANCANCER_CM_ASSOCIATION_AUDIT_REPAIR']).read_text());assert sha(__file__)==repair['auditor_sha256'] and sha(planpath)==repair['unchanged_producer_plan_sha256'] and not (O/'independent_audit.json').exists();inv=json.loads((O/'inventory.json').read_text());assert inv['plan_sha256']==sha(planpath)
    for n,h in inv['outputs'].items():assert sha(O/n)==h
    d=pd.read_csv(O/'all26994_full_control_CM_LR_associations.tsv.gz',sep='\t',keep_default_na=False);assert len(d)==26994;checked=0;cache={}
    for gid,g in d.groupby('group_id'):
        b=dict(np.load(L/'private_group_arrays'/(gid+'.npz'),allow_pickle=False));lr={(l,r):i for i,(l,r) in enumerate(zip(b['ligands'],b['receptors']))};cm={c:i for i,c in enumerate(b['CMs'])}
        for row in g.itertuples():
            j=lr[row.ligand,row.receptor];mask=np.isfinite(b['scores'][:,j,0]);Y=np.log1p(b['scores'][mask,j,0]);blocks=b['patient_blocks'][mask];activity=b['endpoint_excluded_activity' if row.variant=='endpoint_excluded_library' else 'original_H'][mask,cm[row.cm]];raw=np.column_stack([activity,b['source_within_major'][mask],b['target_within_major'][mask],np.log1p(b['source_n'][mask]),np.log1p(b['target_n'][mask])]);sd=raw.std(axis=0,ddof=1);z=np.divide(raw-raw.mean(axis=0),sd,out=np.zeros_like(raw),where=sd>1e-12);X=np.column_stack([np.ones(len(raw)),z]);rank=np.linalg.matrix_rank(X);pinv=np.linalg.pinv(X);diag=np.diag(pinv@X)
            if row.status!='ESTIMABLE_FULL_CONTROL_SPAN_ASSOCIATION':
                if row.status=='CM_EXPOSURE_NOT_IDENTIFIABLE':assert diag[1]<1-1e-6 or sd[0]<=1e-12
                elif row.status=='ILL_CONDITIONED_OR_HIGH_INFLUENCE_DESIGN':assert np.diag(X@pinv).max()>.95-1e-9
                continue
            assert diag[1]>1-1e-7 and exact_int(row.design_rank)==rank;yy=(Y-Y.mean())/Y.std(ddof=1);coef=pinv@yy;err=yy-X@coef;a=pinv[1];G=len(set(blocks));variance=sum((a[blocks==k]@err[blocks==k])**2 for k in set(blocks))*G/(G-1)*(len(yy)-1)/(len(yy)-rank);se=np.sqrt(variance);p=2*t.sf(abs(coef[1]/se),G-1);assert np.isclose(coef[1],float(row.beta),rtol=1e-6,atol=1e-8);assert np.isclose(se,float(row.CR1_se),rtol=1e-6,atol=1e-8);assert np.isclose(p,float(row.analytic_p),rtol=1e-6,atol=1e-10);checked+=1
    for variant,g in d.groupby('variant'):
        pp=pd.to_numeric(g.residual_block_permutation_p,errors='coerce');q=multipletests(pp.fillna(1),method='fdr_bh')[1];valid=pp.notna();assert np.allclose(pd.to_numeric(g.loc[valid,'full_family_BH_q']),q[valid],rtol=1e-10,atol=1e-10)
    # Fresh explicit block-resampled full-X fit for every saved example.
    examples=[]
    for file in sorted((O/'private_resampling_examples').glob('*.npz')):
        b=dict(np.load(file,allow_pickle=False));raw=b['raw_predictors'];Y=b['raw_response'];blocks=b['patient_blocks'];N=len(raw);sd=raw.std(axis=0,ddof=1);X=np.column_stack([np.ones(N),np.divide(raw-raw.mean(axis=0),sd,out=np.zeros_like(raw),where=sd>1e-12)]);yy=(Y-Y.mean())/Y.std(ddof=1);Z=X[:,[0,2,3,4,5]];fit0=Z@np.linalg.lstsq(Z,yy,rcond=None)[0];e0=yy-fit0;co=np.linalg.lstsq(X,yy,rcond=None)[0][1];order=list(dict.fromkeys(blocks));ids=[np.flatnonzero(blocks==k) for k in order];G=len(ids);bins={}
        for i,index in enumerate(ids):bins.setdefault(len(index),[]).append(i)
        rng=np.random.default_rng(int(b['resampling_seed']));perm=[]
        for draw in range(2000):
            vec=e0.copy()
            for members in bins.values():
                shuffle=rng.permutation(members)
                for target,source in zip(members,shuffle):vec[ids[target]]=e0[ids[source]]
            perm.append(np.linalg.lstsq(X,fit0+vec,rcond=None)[0][1])
        Qz,Rz,piv=qr(Z,mode='economic',pivoting=True);rz=np.linalg.matrix_rank(Z);basis=np.column_stack([Qz[:,:rz],X[:,1]])
        pp=(1+np.count_nonzero(abs(np.array(perm))>=abs(co)-1e-12))/2001;draws=rng.integers(0,G,size=(2000,G));bootstrap=[];rank=np.linalg.matrix_rank(X)
        for draw in draws:
            take=np.concatenate([ids[j] for j in draw]);xx=X[take];r=np.linalg.matrix_rank(xx)
            if r!=rank:continue
            pinv=np.linalg.pinv(xx);H=xx@pinv
            if np.diag(H).max()>.95:continue
            # Effective nonzero singular-value condition is a numerical guard.
            ev=np.linalg.eigvalsh(basis[take].T@basis[take]);
            if ev[0]<=1e-12*ev[-1] or ev[-1]/ev[0]>=1e16:continue
            bootstrap.append((pinv@yy[take])[1])
        row=d[(d.case_id==str(b['case_id']))&(d.cm==str(b['cm']))&(d.variant==str(b['variant']))].iloc[0];assert np.isclose(pp,float(row.residual_block_permutation_p));assert len(bootstrap)==exact_int(row.boot_valid),(file.name,len(bootstrap),row.boot_valid)
        if len(bootstrap)>=1800:assert np.allclose(np.quantile(bootstrap,[.025,.975]),[float(row.bootstrap_ci_low),float(row.bootstrap_ci_high)],rtol=1e-6,atol=1e-8)
        examples.append(dict(case_id=str(b['case_id']),cm=str(b['cm']),variant=str(b['variant']),status='PASS_EXPLICIT_FULLX2000_BLOCK_PERMUTATIONS_AND_CASE_BOOTSTRAPS',valid_bootstraps=len(bootstrap)))
    result=dict(audit_repair_plan_sha256=sha(os.environ['PANCANCER_CM_ASSOCIATION_AUDIT_REPAIR']),status='PASS_ALL26994_ROWS_FULL_CONTROL_TARGET_COEFFICIENT_CR1_BH_AND_EXPLICIT_RESAMPLING_EXAMPLES',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),estimated_models_independently_recomputed=checked,resampling_examples=examples,all_required_covariates_retained=True,scientific_or_causal_validation=False)
    (O/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
if __name__=='__main__':main()
