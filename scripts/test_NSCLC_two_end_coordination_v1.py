"""One complete author-subtype coordination family; no directed LR claim."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from threadpoolctl import threadpool_limits


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def correlation(x,y,Z=None,ranks=True):
    a=rankdata(x).astype(float) if ranks else np.array(x,dtype=float,copy=True)
    b=rankdata(y).astype(float) if ranks else np.array(y,dtype=float,copy=True)
    a-=a.mean();b-=b.mean()
    if Z is not None:
        Q=np.linalg.qr(Z,mode='reduced')[0];a-=Q@(Q.T@a);b-=Q@(Q.T@b)
    den=np.linalg.norm(a)*np.linalg.norm(b)
    return float(a@b/den) if den>1e-12 else np.nan


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);I=R/'NSCLC_two_end_coordination_preflight_20261016_v1'
    j=R/'HNSCC_receiver_followup_job_control_20261016_v2/prepare_NSCLC_two_end_coordination_v1/launch_receipt.json'
    launch=json.loads(j.read_text());deadline=time.monotonic()+3600
    while not (I/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists()
        st=p.read_text().rsplit(')',1)[1].split();assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z'
        assert time.monotonic()<deadline;time.sleep(10)
    terminal=json.loads((I/'private_terminal.json').read_text());plan=json.loads((I/'plan_before_any_two_end_effect.json').read_text())
    path=I/'private_all572_planned_pair_designs.tsv';assert sha(path)==terminal['scope_sha256']
    t=pd.read_csv(path,sep='\t');assert len(t)==572
    O=R/'NSCLC_two_end_coordination_results_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    (O/'execution_before_effects.json').write_text(json.dumps(dict(started_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        parent_plan_sha256=sha(I/'plan_before_any_two_end_effect.json'),preflight_terminal_sha256=sha(I/'private_terminal.json'),
        planned_family=len(t),order='Historical endpointA as explanatory rank, endpointB as response rank; no reverse-P selection and no direction-of-signalling claim.',
        caveat='Pathology adjustment is descriptive; conditioning on an outcome can induce associations. No-pathology sensitivity must be interpreted alongside the main result.'),indent=2)+'\n')
    result=[];loo=[]
    with threadpool_limits(limits=1):
        for r in t.itertuples():
            row=r._asdict();row.pop('Index',None)
            row.update(partial_rank_rho=np.nan,HC3_t=np.nan,wild_p=np.nan,CI_low=np.nan,CI_high=np.nan,valid_bootstraps=0,
                raw_spearman=np.nan,no_pathology_control_rho=np.nan,linear_signature_partial_rho=np.nan,old_residual_rank_partial_rho=np.nan,
                nonshared_gene_partial_rho=np.nan,LOO_estimable=0,LOO_same_sign=0,LOO_rho_minimum=np.nan,LOO_rho_maximum=np.nan)
            if r.status!='READY_FOR_FROZEN_TWO_END_TEST':result.append(row);continue
            z=np.load(I/'private_pair_designs'/(r.record+'.npz'),allow_pickle=False);x=z['X'];y=z['Y'];Z=z['Z'];n=len(x)
            assert Z.shape==(n,18) and np.linalg.matrix_rank(Z)==18 and len(set(z['cases']))==n
            xr=rankdata(x).astype(float);yr=rankdata(y).astype(float);xr=(xr-xr.mean())/xr.std();yr=(yr-yr.mean())/yr.std()
            D=np.column_stack([Z,xr]);G=np.linalg.pinv(D);beta=G@yr;e=yr-D@beta;h=np.sum(D*G.T,axis=1)
            se=float(np.sqrt(np.sum((G[-1]*e/(1-h))**2)));observed=float(beta[-1]/se)
            rho=correlation(x,y,Z);assert np.isfinite(rho) and np.isfinite(observed) and se>1e-12
            check=sm.OLS(yr,D).fit(cov_type='HC3',use_t=True);assert np.allclose([beta[-1],se,observed],[check.params[-1],check.bse[-1],check.tvalues[-1]],atol=1e-8,rtol=1e-8)
            seed=np.random.SeedSequence([int(plan['seed']),int(hashlib.sha256(r.record.encode()).hexdigest()[:12],16)])
            random,brng=[np.random.default_rng(s) for s in seed.spawn(2)]
            G0=np.linalg.pinv(Z);fit0=Z@(G0@yr);h0=np.sum(Z*G0.T,axis=1);err=(yr-fit0)/np.sqrt(1-h0);exceed=0
            for start in range(0,int(plan['wild_draws']),200):
                b=min(200,int(plan['wild_draws'])-start);sign=random.choice(np.array([-1.,1.]),size=(n,b));star=fit0[:,None]+err[:,None]*sign
                bs=G@star;res=star-D@bs;variance=np.sum((G[-1,:,None]*res/(1-h[:,None]))**2,axis=0)
                assert (variance>0).all() and np.isfinite(variance).all()
                null=bs[-1]/np.sqrt(variance);exceed+=int(np.sum(np.abs(null)>=abs(observed)-1e-12))
            boots=[]
            for _ in range(int(plan['bootstrap_draws'])):
                take=brng.integers(0,n,n);zz=Z[take]
                if np.linalg.matrix_rank(zz)!=18:continue
                xx=rankdata(x[take]);dd=np.column_stack([zz,xx])
                if np.linalg.matrix_rank(dd)!=19:continue
                v=correlation(x[take],y[take],zz)
                if np.isfinite(v):boots.append(v)
            ci=np.quantile(boots,[.025,.975]) if len(boots)>=int(plan['min_valid_bootstraps']) else [np.nan,np.nan]
            lv=[]
            for k in range(n):
                take=np.arange(n)!=k
                v=correlation(x[take],y[take],Z[take]) if np.linalg.matrix_rank(Z[take])==18 else np.nan
                loo.append(dict(record=r.record,removed_case_index=k,rho=v))
                if np.isfinite(v):lv.append(v)
            linear=correlation(x,y,Z,ranks=False);old_linear=correlation(z['oldA'],z['oldB'],Z,ranks=False)
            assert np.isclose(linear,old_linear,atol=1e-8,rtol=1e-8,equal_nan=True)
            masked=correlation(z['maskedA'],z['maskedB'],Z) if np.isfinite(z['maskedA']).all() and np.isfinite(z['maskedB']).all() else np.nan
            row.update(status='ESTIMABLE_EXPLORATORY_CASE_PROGRAMME_COORDINATION',partial_rank_rho=rho,HC3_t=observed,wild_p=(1+exceed)/(int(plan['wild_draws'])+1),
                CI_low=float(ci[0]),CI_high=float(ci[1]),valid_bootstraps=len(boots),raw_spearman=correlation(x,y),no_pathology_control_rho=correlation(x,y,np.delete(Z,9,axis=1)),
                linear_signature_partial_rho=linear,old_residual_rank_partial_rho=correlation(z['oldA'],z['oldB'],Z),nonshared_gene_partial_rho=masked,
                LOO_estimable=len(lv),LOO_same_sign=int(np.sum(np.sign(lv)==np.sign(rho))),LOO_rho_minimum=min(lv) if lv else np.nan,LOO_rho_maximum=max(lv) if lv else np.nan)
            result.append(row)
            (O/'private_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),planned_models=len(t),models_accounted=len(result),current_record=r.record)))
    out=pd.DataFrame(result);q=multipletests(out.wild_p.fillna(1.),method='fdr_bh')[1]
    out['BH_q_full572']=np.where(out.wild_p.notna(),q,np.nan)
    assert len(out)==572 and out.record.tolist()==t.record.tolist()
    out.to_csv(O/'private_all572_two_end_results.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(loo).to_csv(O/'private_all_estimable_pair_LOO.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='ALL572_TWO_END_COORDINATION_RECORDS_COMPLETE_INDEPENDENT_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),
        planned_family=len(out),estimable_models=int(out.wild_p.notna().sum()),BH_below05=int(out.BH_q_full572.lt(.05).sum()),status_counts=out.status.value_counts().to_dict(),
        result_sha256=sha(O/'private_all572_two_end_results.tsv'),parent_plan_sha256=sha(I/'plan_before_any_two_end_effect.json'),
        directed_LR_or_causal_interaction=False,independent_clinical_prediction=False,fine_identity_transfer_accepted=False)
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
