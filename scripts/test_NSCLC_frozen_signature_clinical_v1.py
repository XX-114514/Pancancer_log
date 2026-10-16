"""One frozen all-signature clinical family; case-block wild/bootstrap.

Association after treatment, not baseline prediction, intervention or causality.
"""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def hc3(y,D,G=None):
    if G is None:G=np.linalg.pinv(D)
    beta=G@y;res=y-D@beta;h=np.sum(D*G.T,axis=1)
    se=float(np.sqrt(np.sum((G[1]*res/(1-h))**2)))
    return float(beta[1]),se,float(beta[1]/se)


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);I=R/'NSCLC_frozen_signature_clinical_preflight_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/preflight_NSCLC_signature_clinical_v1/launch_receipt.json').read_text());deadline=time.monotonic()+3600
    while not (I/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists()
        s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(launch['proc_start_ticks']) and s[0]!='Z'
        assert time.monotonic()<deadline;time.sleep(10)
    terminal=json.loads((I/'private_terminal.json').read_text());plan=json.loads((I/'plan_before_any_clinical_effect.json').read_text())
    path=I/'private_all_planned_signature_subtype_designs.tsv';assert sha(path)==terminal['design_scope_sha256']
    t=pd.read_csv(path,sep='\t');O=R/'NSCLC_frozen_signature_clinical_results_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    (O/'execution_contract_before_effects.json').write_text(json.dumps(dict(started_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
        plan_sha256=sha(I/'plan_before_any_clinical_effect.json'),preflight_terminal_sha256=sha(I/'private_terminal.json'),
        planned_family=len(t),main_BH='All fixed planned signature/subtype slots; unestimable use1 for denominator and retain NA statistical output.',
        wild_bootstrap=int(plan['wild_bootstrap']),case_bootstrap=int(plan['case_bootstrap']),covariates_not_removed=True),indent=2)+'\n')
    result=[];loo=[]
    with threadpool_limits(limits=1):
        for r in t.itertuples():
            row=r._asdict();row.pop('Index',None)
            row.update(beta_MPR_minus_nonMPR=np.nan,HC3_se=np.nan,HC3_t=np.nan,wild_bootstrap_p=np.nan,CI_low=np.nan,CI_high=np.nan,
                valid_bootstraps=0,unadjusted_rank_biserial=np.nan,unadjusted_MannWhitney_p=np.nan,original_residual_score_beta=np.nan,original_residual_score_t=np.nan,
                LOO_estimable=0,LOO_same_sign=0,LOO_beta_minimum=np.nan,LOO_beta_maximum=np.nan)
            if r.status!='DESIGN_READY_NO_EFFECT_COMPUTED':result.append(row);continue
            z=np.load(I/'private_frozen_designs'/(r.record+'.npz'),allow_pickle=False);D=z['D'];Z=z['Z'];Y=z['Y'];old=z['Y_original_frozen_residual'];x=z['X'];n=len(Y)
            assert len(set(z['case_blocks']))==n and np.linalg.matrix_rank(D)==14 and np.linalg.matrix_rank(Z)==13
            G=np.linalg.pinv(D);h=np.sum(D*G.T,axis=1);beta,se,observed=hc3(Y,D,G)
            assert se>1e-12 and np.isfinite(observed)
            independent=sm.OLS(Y,D).fit(cov_type='HC3',use_t=True)
            assert np.allclose([beta,se,observed],[independent.params[1],independent.bse[1],independent.tvalues[1]],atol=1e-8,rtol=1e-8)
            old_beta,old_se,old_t=hc3(old,D,G)
            assert np.allclose([old_beta,old_se,old_t],[beta,se,observed],atol=1e-8,rtol=1e-8),'Changed score parameterization changes planned fully adjusted inference'
            seed=np.random.SeedSequence([int(plan['seed']),int(hashlib.sha256(r.record.encode()).hexdigest()[:12],16)])
            random,bootrng=[np.random.default_rng(s) for s in seed.spawn(2)]
            G0=np.linalg.pinv(Z);fitted=Z@(G0@Y);h0=np.sum(Z*G0.T,axis=1);error=(Y-fitted)/np.sqrt(1-h0)
            exceed=0;nulls=int(plan['wild_bootstrap'])
            for start in range(0,nulls,200):
                b=min(200,nulls-start);sign=random.choice(np.array([-1.,1.]),size=(n,b));star=fitted[:,None]+error[:,None]*sign
                bs=G@star;e=star-D@bs;variance=np.sum((G[1,:,None]*e/(1-h[:,None]))**2,axis=0)
                assert (variance>0).all() and np.isfinite(variance).all()
                null=bs[1]/np.sqrt(variance);exceed+=int(np.sum(np.abs(null)>=abs(observed)-1e-12))
            bp=[]
            for _ in range(int(plan['case_bootstrap'])):
                take=bootrng.integers(0,n,n)
                if np.linalg.matrix_rank(D[take])!=14:continue
                value=float(np.linalg.lstsq(D[take],Y[take],rcond=None)[0][1])
                if np.isfinite(value):bp.append(value)
            ci=np.quantile(bp,[.025,.975]) if len(bp)>=int(plan['minimum_valid_bootstrap']) else [np.nan,np.nan]
            a=Y[x==1];b=Y[x==0];u=mannwhitneyu(a,b,alternative='two-sided',method='asymptotic')
            vv=[]
            for k in range(n):
                take=np.arange(n)!=k;value=float(np.linalg.lstsq(D[take],Y[take],rcond=None)[0][1]) if np.linalg.matrix_rank(D[take])==14 else np.nan
                loo.append(dict(record=r.record,removed_case_index=k,beta=value))
                if np.isfinite(value):vv.append(value)
            row.update(status='ESTIMABLE_POSTTREATMENT_SIGNATURE_CASE_ASSOCIATION',beta_MPR_minus_nonMPR=beta,HC3_se=se,HC3_t=observed,wild_bootstrap_p=(1+exceed)/(nulls+1),
                CI_low=float(ci[0]),CI_high=float(ci[1]),valid_bootstraps=len(bp),unadjusted_rank_biserial=float(2*u.statistic/(len(a)*len(b))-1),unadjusted_MannWhitney_p=float(u.pvalue),
                original_residual_score_beta=old_beta,original_residual_score_t=old_t,LOO_estimable=len(vv),LOO_same_sign=int(np.sum(np.sign(vv)==np.sign(beta))),
                LOO_beta_minimum=min(vv) if vv else np.nan,LOO_beta_maximum=max(vv) if vv else np.nan)
            result.append(row)
            (O/'private_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),planned_models=len(t),models_accounted=len(result),current_record=r.record)))
    out=pd.DataFrame(result);p=out.wild_bootstrap_p.fillna(1.).to_numpy();q=multipletests(p,method='fdr_bh')[1]
    out['BH_q_full_planned_family']=np.where(out.wild_bootstrap_p.notna(),q,np.nan)
    assert len(out)==len(t) and out.record.tolist()==t.record.tolist()
    out.to_csv(O/'private_all_planned_signature_clinical_results.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(loo).to_csv(O/'private_all_estimable_clinical_LOO.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='ONE_FROZEN_ALL_SIGNATURE_CLINICAL_FAMILY_COMPLETE_INDEPENDENT_RESAMPLING_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),
        planned_family=len(out),estimated_models=int(out.wild_bootstrap_p.notna().sum()),BH_below05=int(out.BH_q_full_planned_family.lt(.05).sum()),
        status_counts=out.status.value_counts().to_dict(),result_sha256=sha(O/'private_all_planned_signature_clinical_results.tsv'),
        preflight_contract_sha256=sha(I/'plan_before_any_clinical_effect.json'),coefficient_t_raw_residual_invariance_verified=True,
        confirmatory_external_validation=False,baseline_prediction=False,LR_mechanism_verified=False)
    (O/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
