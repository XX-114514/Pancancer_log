"""Full independent QR/OLS replay of the matched-case clinical family."""
import hashlib,json,os,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy.linalg import qr,solve_triangular
from scipy.stats import beta as beta_dist


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def close(a,b):assert np.allclose(a,b,atol=1e-8,rtol=1e-8,equal_nan=True)


def coefficient(y,d):
    s=np.linalg.svd(d,compute_uv=False);rank=int((s>s[0]*1e-12).sum())
    if rank!=d.shape[1] or rank!=np.linalg.matrix_rank(d) or len(y)-rank<1 or s[0]/s[-1]>=1e8:return np.full(2,np.nan)
    q,rr=qr(d,mode='economic');return solve_triangular(rr,q.T@y)[-2:]


def exceedances(y,d,target,draws,seed):
    q,rr=qr(d,mode='economic');inverse=solve_triangular(rr,q.T);beta=inverse@y;e=y-d@beta;h=np.sum(q*q,axis=1)
    variance=np.sum((inverse[target]*e/(1-h))**2);observed=beta[target]/np.sqrt(variance)
    z=np.delete(d,target,axis=1);q0,_=qr(z,mode='economic');mu=q0@(q0.T@y);e0=(y-mu)/np.sqrt(1-np.sum(q0*q0,axis=1))
    rng=np.random.default_rng(seed);count=0
    for start in range(0,draws,500):
        signs=np.where(rng.random((min(500,draws-start),len(y)))<.5,-1.,1.).T
        yy=mu[:,None]+e0[:,None]*signs;b=inverse@yy;ee=yy-d@b
        variance=np.sum((inverse[target,:,None]*ee/(1-h[:,None]))**2,axis=0);tt=b[target]/np.sqrt(variance)
        assert np.isfinite(tt).all();count+=int((np.abs(tt)>=abs(observed)-1e-12).sum())
    return count


def bh(p):
    order=np.argsort(p,kind='mergesort');vals=np.asarray(p)[order];v=np.minimum.accumulate((vals*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    out=np.empty(len(p));out[order]=np.minimum(v,1);return out


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'NSCLC_paired_programme_clinical_increment_results_20261016_v1';pre=r/'NSCLC_paired_programme_clinical_increment_preflight_20261016_v1'
    job=r/'HNSCC_receiver_followup_job_control_20261016_v2/run_NSCLC_paired_clinical_increment_v1';receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+21600
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split()
        assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(10)
    terminal=json.loads((out/'private_terminal.json').read_text())
    for n,h in terminal['files'].items():assert sha(out/n)==h
    plan=json.loads((out/'plan_before_true_paired_clinical_effects.json').read_text());cfg=json.loads((pre/'plan_before_any_paired_clinical_increment_effect.json').read_text())
    assert plan['family']==cfg['family']==572 and cfg['wild_draws']==plan['wild_draws']==199999
    for path,key in [(pre/'private_independent_input_audit.json','input_audit_sha256'),(r/'paired_clinical_increment_calibration_20261016_v2/private_terminal.json','calibration_sha256')]:
        assert sha(path)==plan[key] and json.loads(path.read_text())['status'].startswith('PASS_')
    t=pd.read_csv(out/'private_all572_paired_clinical_results.tsv',sep='\t');scope=pd.read_csv(pre/'private_all572_scope.tsv',sep='\t')
    assert len(t)==572 and t.record.tolist()==scope.record.tolist();hashes=pd.read_csv(pre/'private_ready_design_SHA256.tsv',sep='\t').set_index('record').sha256
    rh=pd.read_csv(out/'private_resample_SHA256.tsv',sep='\t').set_index('record').sha256;count=0;maximum=0.;replayed_draws=0
    for i,row in enumerate(t.itertuples()):
        if row.status!='READY_NO_CLINICAL_EFFECT_COMPUTED':
            assert row.association_status=='INPUT_UNAVAILABLE' and pd.isna(row.pair_IUT_P) and pd.isna(row.BH_q_pair572);continue
        assert row.association_status=='ESTIMATED_PAIRED_CONDITIONAL_CLINICAL_EXPLORATORY';p=pre/'private_designs'/(row.record+'.npz');assert sha(p)==hashes.loc[row.record]
        a=np.load(p,allow_pickle=False);y=a['Y'];d=a['D'];ref=sm.OLS(y,d).fit(cov_type='HC3');seed=cfg['seed']+i*10
        for key,target,j in [('A',-2,0),('B',-1,1)]:
            close(getattr(row,key+'_beta'),ref.params[target]);close(getattr(row,key+'_HC3_t'),ref.tvalues[target]);close(getattr(row,key+'_HC3_se'),ref.bse[target])
            count_w=exceedances(y,d,target,199999,seed+j);assert count_w==getattr(row,key+'_wild_exceedances');replayed_draws+=199999
            close(getattr(row,key+'_P'),(1+count_w)/200000)
            upper=1. if count_w==199999 else beta_dist.ppf(1-.01/1144,count_w+1,199999-count_w)
            close(getattr(row,key+'_MC_upper_P'),min(1.,(1+199999*upper)/200000))
            ds=np.delete(d,-1 if j==0 else -2,axis=1);close(getattr(row,key+'_matched_single_end_beta'),sm.OLS(y,ds).fit().params[-1])
            reduced=np.delete(d,target,axis=1);e0=sm.OLS(y,reduced).fit().resid;close(getattr(row,key+'_incremental_partial_R2'),1-(ref.resid@ref.resid)/(e0@e0))
        close(row.pair_IUT_P,max(row.A_P,row.B_P));close(row.pair_MC_upper_IUT_P,max(row.A_MC_upper_P,row.B_MC_upper_P))
        p=out/'private_resamples'/(row.record+'.npz');assert sha(p)==rh.loc[row.record];s=np.load(p,allow_pickle=False)
        rng=np.random.default_rng(seed+2);boot=np.full((2000,2),np.nan)
        for k in range(2000):
            take=rng.integers(0,len(y),len(y));boot[k]=coefficient(y[take],d[take])
        leave=np.full((len(y),2),np.nan)
        for k in range(len(y)):
            use=np.arange(len(y))!=k;leave[k]=coefficient(y[use],d[use])
        for key,expected in [('bootstrap_coefficients',boot),('LOO_coefficients',leave)]:
            assert np.array_equal(np.isfinite(s[key]),np.isfinite(expected));close(s[key],expected)
            if np.isfinite(expected).any():maximum=max(maximum,float(np.nanmax(np.abs(s[key]-expected))))
        for j,key in enumerate(['A','B']):
            valid=np.isfinite(boot[:,j]);assert valid.sum()==getattr(row,key+'_valid_bootstraps')
            q=np.quantile(boot[valid,j],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan];close(q,[getattr(row,key+'_CI_low'),getattr(row,key+'_CI_high')])
            assert np.isfinite(leave[:,j]).sum()==getattr(row,key+'_LOO_valid')
            stable=np.isfinite(leave[:,j]).all() and (np.sign(leave[:,j])==np.sign(ref.params[-2+j])).all();assert bool(stable)==getattr(row,key+'_LOO_same_sign')
            close([np.nanmin(leave[:,j]),np.nanmax(leave[:,j])],[getattr(row,key+'_LOO_min'),getattr(row,key+'_LOO_max')])
        mask=np.column_stack([a['masked_A'],a['masked_B']])
        if row.nonshared_conditional_status=='DESCRIPTIVE_ONLY_NO_NEW_P':
            dm=np.column_stack([d[:,:-2],(mask-mask.mean(0))/mask.std(0)]);close(coefficient(y,dm),[row.nonshared_A_beta,row.nonshared_B_beta])
        count+=1;(out/'private_independent_audit_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),models_checked=count,wild_draws_independently_replayed=replayed_draws))+'\n')
    valid=t.pair_IUT_P.notna().to_numpy();p=t.pair_IUT_P.fillna(1).to_numpy();upper=t.pair_MC_upper_IUT_P.fillna(1).to_numpy()
    for key,values in [('BH_q_pair572',bh(p)),('BY_q_pair572',np.minimum(1,bh(p)*np.sum(1/np.arange(1,573)))),('MC_upper_BH_q_pair572',bh(upper))]:
        close(t.loc[valid,key],values[valid]);assert t.loc[~valid,key].isna().all()
    assert count==terminal['estimated_models']==len(rh)
    proof=dict(status='PASS_ALL572_CONDITIONAL_POINTS_ALL_FIXED_WILD_BOOTSTRAP_LOO_SENSITIVITIES_AND_IUT_BH_BY_MC_REPLAY',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_code_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),
               estimated_models=count,all_conditional_wild_draws_independently_replayed=replayed_draws,maximum_resample_coefficient_difference=maximum,clinical_prediction_or_causal_validation=False)
    p=out/'private_independent_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
