"""Pipeline every private numeric result through independent QR/random-stream replay."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy.stats import beta as beta_dist


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rngseed(identifier,stream):return int.from_bytes(hashlib.sha256((identifier+'|'+stream+'|20261016').encode()).digest()[:8],'little')


def same(a,b):assert np.isclose(a,b,rtol=1e-8,atol=1e-9,equal_nan=True),(a,b)


def adjusted(p,by=False):
    p=np.asarray(p,float);order=np.argsort(p,kind='stable');n=len(p);ranked=p[order]*n/np.arange(1,n+1);q=np.minimum(1.,np.minimum.accumulate(ranked[::-1])[::-1]);result=np.empty(n);result[order]=q
    return np.minimum(1.,result*sum(1/k for k in range(1,n+1))) if by else result


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);pre=r/'HNSCC_function_change_clinical_preflight_20261016_v1';source=r/'HNSCC_function_change_clinical_results_20261016_v1';out=r/'HNSCC_function_change_clinical_audit_20261016_v1';assert not out.exists();out.mkdir(mode=0o700);(out/'private_model_audits').mkdir(mode=0o700)
    rc=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/run_HNSCC_function_change_clinical_v1/launch_receipt.json').read_text());deadline=time.monotonic()+43200
    def wait(path):
        while not path.exists():
            stat=Path('/proc')/str(rc['pid'])/'stat';assert stat.exists();a=stat.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(3)
    wait(source/'plan_before_observed_clinical_effects.json');plan=json.loads((source/'plan_before_observed_clinical_effects.json').read_text());assert plan['code_sha256']==rc['source_sha256'];model=pd.read_csv(pre/'private_unique_numeric_designs.tsv',sep='\t').sort_values('model_id').reset_index(drop=True);assert len(model)==plan['unique_numeric_designs'];B=int(plan['parameters']['wild_draws']);draws=int(plan['parameters']['bootstrap_draws']);total_wild=0;total_boot=0;maximum=0.;record_cache={}
    for number,row in enumerate(model.itertuples()):
        identifier=row.model_id;point_path=source/'private_models'/(identifier+'.json');wait(point_path);point=json.loads(point_path.read_text());p=pre/'private_designs'/(identifier+'.npz');assert sha(p)==row.design_sha256==point['design_sha256'];archive=np.load(p,allow_pickle=False);d=archive['D'];y=archive['Y'];strata=archive['strata'];assert point['n']==len(y) and point['residual_df']==len(y)-11
        saved=source/'private_models'/(identifier+'.npz');assert sha(saved)==point['resamples_sha256'];resamples=np.load(saved,allow_pickle=False);observed_boot=resamples['bootstrap_coefficients'];observed_loo=resamples['LOO_coefficients'];assert observed_boot.shape==(draws,) and observed_loo.shape==(len(y),)
        if point['status'].startswith('ESTIMATED_'):
            m=sm.OLS(y,d).fit(cov_type='HC3');same(point['beta'],m.params[-1]);same(point['HC3_t'],m.tvalues[-1]);same(point['HC3_se'],m.bse[-1]);q,R=np.linalg.qr(d,mode='reduced');solver=np.linalg.solve(R,q.T);h=np.square(q).sum(axis=1);z=d[:,:-1];q0,_=np.linalg.qr(z,mode='reduced');mu=q0@(q0.T@y);e0=(y-mu)/np.sqrt(1-np.square(q0).sum(axis=1));rng=np.random.default_rng(rngseed(identifier,'wild'));exceed=0
            for start in range(0,B,1000):
                signs=np.where(rng.random((min(1000,B-start),len(y)))<.5,-1.,1.).T;yy=mu[:,None]+e0[:,None]*signs;coefficients=solver@yy;errors=yy-d@coefficients;variance=np.sum(np.square(solver[-1,:,None]*errors/(1-h[:,None])),axis=0);statistic=coefficients[-1]/np.sqrt(variance);assert np.isfinite(statistic).all();exceed+=int((abs(statistic)>=abs(m.tvalues[-1])-1e-12).sum())
            assert point['wild_exceedances']==exceed;pv=(1+exceed)/(1+B);upper=1. if exceed==B else float(beta_dist.isf(.01/len(model),exceed+1,B-exceed));upper=min(1.,(1+B*upper)/(1+B));same(point['P'],pv);same(point['MC_upper_P'],upper);total_wild+=B
            rng=np.random.default_rng(rngseed(identifier,'bootstrap'));take=np.empty((draws,len(y)),int);offset=0
            for level in sorted(set(strata)):
                ids=np.flatnonzero(strata==level);take[:,offset:offset+len(ids)]=ids[rng.integers(0,len(ids),(draws,len(ids)))];offset+=len(ids)
            expected_boot=np.full(draws,np.nan)
            for start in range(0,draws,128):
                end=min(draws,start+128);indices=take[start:end];xx=d[indices];s=np.linalg.svd(xx,compute_uv=False);rank=(s>s[:,:1]*1e-12).sum(axis=1);machine=(s>s[:,:1]*np.finfo(float).eps*max(d.shape)).sum(axis=1);good=(rank==11)&(rank==machine)&(s[:,0]/np.maximum(s[:,-1],np.finfo(float).tiny)<1e8)
                qq,rr=np.linalg.qr(xx[good],mode='reduced');rhs=np.einsum('bni,bn->bi',qq,y[indices[good]]);coefficients=np.linalg.solve(rr,rhs[...,None])[...,0];expected_boot[np.arange(start,end)[good]]=coefficients[:,-1]
            assert np.array_equal(np.isfinite(expected_boot),np.isfinite(observed_boot));assert np.allclose(expected_boot,observed_boot,rtol=1e-8,atol=1e-8,equal_nan=True);maximum=max(maximum,float(np.nanmax(abs(expected_boot-observed_boot))));total_boot+=draws;valid=np.isfinite(expected_boot);assert int(valid.sum())==point['valid_bootstraps'];limits=np.quantile(expected_boot[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan];same(point['CI_low'],limits[0]);same(point['CI_high'],limits[1]);leave=np.full(len(y),np.nan)
            for j in range(len(y)):
                keep=np.arange(len(y))!=j;dd=d[keep];s=np.linalg.svd(dd,compute_uv=False);r1=int(np.sum(s>s[0]*1e-12));r2=int(np.sum(s>s[0]*np.finfo(float).eps*max(dd.shape)))
                if r1==r2==11 and len(dd)-11>=1 and s[0]/s[-1]<1e8:leave[j]=sm.OLS(y[keep],dd).fit().params[-1]
            assert np.allclose(leave,observed_loo,rtol=1e-8,atol=1e-8,equal_nan=True);assert point['LOO_valid']==int(np.isfinite(leave).sum()) and point['LOO_same_sign']==bool(np.isfinite(leave).all() and (np.sign(leave)==np.sign(m.params[-1])).all());same(point['LOO_min'],np.nanmin(leave));same(point['LOO_max'],np.nanmax(leave));reduced=sm.OLS(y,z).fit().resid;full=sm.OLS(y,d).fit().resid;same(point['incremental_partial_R2'],1-(full@full)/(reduced@reduced))
        else:
            m=sm.OLS(y,d).fit(cov_type='HC3');assert not np.isfinite(m.bse[-1]) or m.bse[-1]<=0;assert np.isnan(observed_boot).all() and np.isnan(point['P'])
        record_cache[identifier]=point;proof=dict(status='PASS_NUMERIC_MODEL_POINT_ALL_WILD_BOOTSTRAP_LOO',model_id=identifier,producer_point_sha256=sha(point_path),resamples_sha256=sha(saved),auditor_sha256=sha(__file__));(out/'private_model_audits'/(identifier+'.json')).write_text(json.dumps(proof)+'\n');(out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),unique_models_audited=number+1,unique_models_planned=len(model),wild_draws_replayed=total_wild,bootstrap_draws_replayed=total_boot))+'\n')
    wait(source/'private_terminal.json');terminal=json.loads((source/'private_terminal.json').read_text());assert terminal['code_sha256']==plan['code_sha256']
    for name,h in terminal['files'].items():assert sha(source/name)==h
    original=pd.read_csv(pre/'private_all8000_clinical_scope.tsv',sep='\t',keep_default_na=False).set_index('slot');result=pd.read_csv(source/'private_all8000_longitudinal_function_clinical_results.tsv',sep='\t',keep_default_na=False,na_values=['NA']);assert len(result)==8000 and result.slot.is_unique and set(result.slot)==set(original.index)
    fields=['beta','HC3_t','HC3_se','P','MC_upper_P','CI_low','CI_high','valid_bootstraps','LOO_valid','LOO_same_sign','LOO_min','LOO_max','incremental_partial_R2']
    for row in result.itertuples():
        old=original.loc[row.slot];assert row.family==old.family and row.pool_A==old.pool_A and row.pool_B==old.pool_B and row.mask_id==old.mask_id and row.reference==old.reference and row.status==old.status
        endpoints=['A'] if row.family=='single3200' else ['A','B'];all_ready=old.status=='READY_NO_OBSERVED_CLINICAL_EFFECT'
        if all_ready:
            for label in endpoints:
                identifier=old[label+'_model'];point=record_cache[identifier];all_ready &= point['status'].startswith('ESTIMATED_')
                for field in fields:
                    got=getattr(row,label+'_'+field)
                    if isinstance(point[field],bool):assert bool(got)==point[field]
                    else:same(got,point[field])
            if all_ready:same(row.clinical_P,max(record_cache[old[label+'_model']]['P'] for label in endpoints));same(row.MC_upper_clinical_P,max(record_cache[old[label+'_model']]['MC_upper_P'] for label in endpoints));assert row.association_status.startswith('ESTIMATED_')
        if not all_ready:assert not row.association_status.startswith('ESTIMATED_') and np.isnan(row.clinical_P)
    for family,g in result.groupby('family'):
        assert len(g)==(4800 if family=='pair4800' else 3200);valid=g.association_status.str.startswith('ESTIMATED_').to_numpy();p=g.clinical_P.fillna(1).to_numpy(float);upper=g.MC_upper_clinical_P.fillna(1).to_numpy(float)
        for column,wanted in [('BH_q',adjusted(p)),('BY_q',adjusted(p,True)),('MC_upper_BH_q',adjusted(upper))]:wanted[~valid]=np.nan;assert np.allclose(g[column],wanted,rtol=1e-9,atol=1e-10,equal_nan=True)
    proof=dict(status='PASS_ALL_UNIQUE_MODELS_FULL_RANDOM_STREAMS_BOOTSTRAPS_LOO_AND_ALL4800_3200_FAMILY_MAPPINGS',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),unique_models_audited=len(model),logical_slots_checked=len(result),wild_draws_replayed=total_wild,bootstrap_draws_replayed=total_boot,maximum_bootstrap_coefficient_difference=maximum,no_direct_communication_independent_validation_or_prediction_claim=True)
    (out/'private_terminal.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
