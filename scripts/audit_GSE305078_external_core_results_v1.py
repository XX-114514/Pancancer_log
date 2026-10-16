"""Independent OLS/pseudoinverse replay of every frozen external result."""
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import rankdata, beta as beta_distribution


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def independent_point(x, y, c, main=False):
    n = len(x)
    if n < (20 if main else 1) or not np.isfinite(c).all() or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('input')
    sd = c.std(0); z = np.column_stack([np.ones(n), (c-c.mean(0))/np.where(sd > 1e-12, sd, 1)])
    s = np.linalg.svd(z, compute_uv=False); rz = int((s > s[0]*1e-12).sum())
    if rz != int((s > s[0]*np.finfo(float).eps*max(z.shape)).sum()):
        raise ValueError('ambiguous')
    xx = rankdata(x).astype(float); yy = rankdata(y).astype(float)
    if xx.std() <= 1e-12 or np.ptp(y) <= 1e-12:
        raise ValueError('constant')
    xx = (xx-xx.mean())/xx.std(); yy = (yy-yy.mean())/yy.std()
    d = np.column_stack([z, xx]); ss = np.linalg.svd(d, compute_uv=False); rd = int((ss > ss[0]*1e-12).sum())
    if rd != int((ss > ss[0]*np.finfo(float).eps*max(d.shape)).sum()):
        raise ValueError('ambiguous')
    if rd != rz+1 or n-rd < (10 if main else 1) or ss[0]/ss[rd-1] >= 1e8:
        raise ValueError('rank_df_condition')
    rx = xx-z@np.linalg.lstsq(z,xx,rcond=1e-12)[0]
    ry = yy-z@np.linalg.lstsq(z,yy,rcond=1e-12)[0]
    if rx@rx <= 1e-16*(xx@xx) or ry@ry <= 1e-16*(yy@yy):
        raise ValueError('residual')
    rho = float(rx@ry/(np.linalg.norm(rx)*np.linalg.norm(ry)))
    result = dict(rho=rho, z=z, d=d, y=yy)
    if main:
        inverse = np.linalg.pinv(d, rcond=1e-12); hat = np.sum(d*inverse.T,axis=1)
        if hat.max() >= .999999:
            raise ValueError('leverage')
        fit = sm.OLS(yy,d).fit(cov_type='HC3')
        if not np.isfinite(fit.tvalues[-1]) or fit.cov_params()[-1,-1] <= 0:
            raise ValueError('variance')
        result.update(inverse=inverse, hat=hat, coefficient=float(fit.params[-1]), t=float(fit.tvalues[-1]), variance=float(fit.cov_params()[-1,-1]))
    return result


def bootstrap(x,y,c,b,seed):
    rng=np.random.default_rng(seed); values=np.full(b,np.nan)
    for i in range(b):
        take=rng.integers(0,len(x),size=len(x))
        try: values[i]=independent_point(x[take],y[take],c[take])['rho']
        except ValueError: pass
    return values


def wild_exceed(g,draws,seed):
    z=g['z']; d=g['d']; inverse=g['inverse']; h=g['hat']
    inverse_z=np.linalg.pinv(z,rcond=1e-12)
    mu=z@np.linalg.lstsq(z,g['y'],rcond=1e-12)[0]
    h0=np.sum(z*inverse_z.T,axis=1); adjusted=(g['y']-mu)/np.sqrt(1-h0)
    rng=np.random.default_rng(seed); exceed=0
    for start in range(0,draws,500):
        signs=np.where(rng.random((min(500,draws-start),len(mu)))<.5,-1.,1.).T
        y=mu[:,None]+adjusted[:,None]*signs
        coef=inverse@y; error=y-d@coef
        variance=np.sum((inverse[-1,:,None]*error/(1-h[:,None]))**2,axis=0)
        t=coef[-1]/np.sqrt(variance)
        assert np.isfinite(t).all()
        exceed+=int((np.abs(t)>=abs(g['t'])-1e-12).sum())
    return exceed


def bh(p):
    order=np.argsort(p,kind='mergesort'); ordered=np.asarray(p)[order]
    q=np.minimum.accumulate((ordered*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    result=np.empty(len(p));result[order]=np.minimum(q,1);return result


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out=root/'GSE305078_external_ligand_programme_core_results_20261016_v1'
    job=root/'HNSCC_receiver_followup_job_control_20261016_v2/run_GSE305078_external_core_v1'
    receipt=json.loads((job/'launch_receipt.json').read_text()); deadline=time.monotonic()+21600
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists()
        s=p.read_text().rsplit(')',1)[1].split()
        assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline
        time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text())
    for name,digest in terminal['files'].items():assert sha(out/name)==digest
    pre=root/'GSE305078_external_ligand_programme_core_20261016_v1'
    cfg=json.loads((root/'external_ligand_programme_core_followup_config_20261016_v1.json').read_text())
    data=pd.read_csv(out/'private_all676_external_core_results.tsv',sep='\t')
    assert len(data)==676 and data.model_id.is_unique
    design_hash=pd.read_csv(pre/'private_ready_design_SHA256.tsv',sep='\t').set_index('model_id').sha256
    resample_hash=pd.read_csv(out/'private_resample_SHA256.tsv',sep='\t').set_index('model_id').sha256
    counts=0; maximum_point=0.; maximum_resample=0.
    for i,row in enumerate(data.itertuples()):
        if row.status!='READY_FOR_INDEPENDENT_INPUT_AUDIT':
            assert row.association_status=='INPUT_UNAVAILABLE' and pd.isna(row.wild_P)
            continue
        p=pre/'private_designs'/(row.model_id+'.npz');assert sha(p)==design_hash.loc[row.model_id]
        a=np.load(p,allow_pickle=False);x=a['X'];y=a['Y'];c=a['C_raw']
        try:g=independent_point(x,y,c,main=True)
        except ValueError:
            assert row.association_status.startswith('UNESTIMABLE_') and pd.isna(row.wild_P);continue
        assert row.association_status=='ESTIMATED_FIXED_FAMILY_WILD_APPROXIMATION'
        for name,expected in [('partial_rank_rho',g['rho']),('coefficient',g['coefficient']),('HC3_t',g['t']),('HC3_variance',g['variance'])]:
            got=getattr(row,name);assert np.isclose(got,expected,rtol=1e-8,atol=1e-8),name
            maximum_point=max(maximum_point,abs(got-expected))
        seed=cfg['seed']+i*10
        exceed=wild_exceed(g,cfg['wild_draws'],seed)
        assert exceed==int(row.wild_exceedances) and row.wild_draws==cfg['wild_draws']
        assert np.isclose(row.wild_P,(1+exceed)/(cfg['wild_draws']+1),atol=1e-14)
        upper=1. if exceed==cfg['wild_draws'] else beta_distribution.ppf(1-.01/676,exceed+1,cfg['wild_draws']-exceed)
        upper=(1+cfg['wild_draws']*upper)/(cfg['wild_draws']+1)
        assert np.isclose(upper,row.MC_upper_P,atol=1e-12)
        p=out/'private_resamples'/(row.model_id+'.npz');assert sha(p)==resample_hash.loc[row.model_id]
        original=np.load(p,allow_pickle=False)
        boot=bootstrap(x,y,c,cfg['bootstrap_draws'],seed+1)
        loo=np.full(len(x),np.nan)
        for k in range(len(x)):
            use=np.arange(len(x))!=k
            try:loo[k]=independent_point(x[use],y[use],c[use])['rho']
            except ValueError:pass
        for key,expected in [('bootstrap_rho',boot),('LOO_rho',loo)]:
            got=original[key];assert np.array_equal(np.isfinite(got),np.isfinite(expected))
            assert np.allclose(got,expected,rtol=1e-8,atol=1e-8,equal_nan=True),key
            if np.isfinite(expected).any():maximum_resample=max(maximum_resample,float(np.nanmax(np.abs(got-expected))))
        valid=int(np.isfinite(boot).sum());assert valid==row.bootstrap_valid
        if valid>=1800:
            assert np.allclose(np.quantile(boot[np.isfinite(boot)],[.025,.975]),[row.bootstrap_low,row.bootstrap_high],atol=1e-8)
        else:assert pd.isna(row.bootstrap_low) and pd.isna(row.bootstrap_high)
        assert row.LOO_valid==np.isfinite(loo).sum()
        assert bool(row.LOO_direction_stable)==bool(np.isfinite(loo).all() and np.all(np.sign(loo)==np.sign(g['rho'])))
        for key,extra in [('sender_programme_adjusted_rho',a['sender_programme']),('sender_same_regulator_adjusted_rho',a['sender_same_regulator'])]:
            expected=np.nan
            if np.isfinite(extra).all():
                try:expected=independent_point(x,y,np.column_stack([c,extra]))['rho']
                except ValueError:pass
            assert np.isclose(getattr(row,key),expected,rtol=1e-8,atol=1e-8,equal_nan=True)
        raw=float(np.corrcoef(rankdata(x),rankdata(y))[0,1]);assert np.isclose(raw,row.unadjusted_spearman)
        counts+=1
        (out/'private_independent_audit_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),models_checked=counts,last_model=row.model_id))+'\n')
    p=data.wild_P.fillna(1).to_numpy();upper=data.MC_upper_P.fillna(1).to_numpy()
    for name,expected in [('BH_q_complete676',bh(p)),('BY_q_complete676',np.minimum(1,bh(p)*np.sum(1/np.arange(1,677)))),('MC_upper_BH_q_complete676',bh(upper))]:
        valid=data.wild_P.notna().to_numpy();assert np.allclose(data[name].to_numpy()[valid],expected[valid],atol=1e-12)
        assert data.loc[~valid,name].isna().all()
    assert counts==terminal['estimated_models']==len(resample_hash)
    # Verify every same-core discovery bridge point without new P/BH.
    bridge=pd.read_csv(out/'private_all_discovery_core10_bridge_points.tsv',sep='\t')
    tfroot=root/'NSCLC_gene_disjoint_receiver_TF_measurement_20261016_v1'
    tf_scope=pd.read_csv(tfroot/'private_all_TF_measurement_scope.tsv',sep='\t')
    for row in bridge.itertuples():
        p=root/'NSCLC_automatic_ligand_response_preflight_20261016_v1/private_frozen_designs'/(row.original_record+'.npz')
        original=np.load(p,allow_pickle=False);y=original['Y']
        if row.readout!='masked_programme':
            rec=tf_scope[tf_scope.record.eq(row.original_record)&tf_scope.TF.eq(cfg['secondary_regulator'])].iloc[0]
            y=np.load(tfroot/'private_frozen_designs'/(rec.TF_record+'.npz'),allow_pickle=False)['Y']
        assert np.isclose(independent_point(original['X'],y,original['covariates'],True)['rho'],row.full18_partial_rank_rho,atol=1e-9)
        assert np.isclose(independent_point(original['X'],y,original['covariates'][:,8:],True)['rho'],row.core10_partial_rank_rho,atol=1e-9)
    proof=dict(status='PASS_ALL676_EXTERNAL_POINTS_HC3_COMPLETE_FIXED_WILD_REPLAYS_ALL_BOOTSTRAP_LOO_SENSITIVITIES_BH_BY_MC_AND_DISCOVERY_BRIDGE',
               finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),
               complete_family=len(data),estimated_models=counts,maximum_point_difference=maximum_point,maximum_resample_difference=maximum_resample,
               all_wild_draws_replayed=True,all_bootstrap_and_LOO_replayed=True,clinical_prediction_or_causal_validation=False)
    p=out/'private_independent_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
