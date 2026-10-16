"""Independent complete numerical audit, with three full random-stream replays."""
import hashlib,json,math,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def close(a,b):assert np.allclose(a,b,rtol=1e-7,atol=1e-7,equal_nan=True),(a,b)

def rho(x,y,Z):
    a=sm.OLS(rankdata(x),Z).fit().resid;b=sm.OLS(rankdata(y),Z).fit().resid
    return float(np.corrcoef(a,b)[0,1]) if np.linalg.norm(a)*np.linalg.norm(b)>1e-12 else np.nan

def manual_q(p,by=False):
    p=np.asarray(p,float);n=len(p);order=np.argsort(p,kind='stable');factor=math.fsum(1/i for i in range(1,n+1)) if by else 1.
    out=np.empty(n);v=p[order]*n*factor/np.arange(1,n+1);out[order]=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1]);return out

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'NSCLC_automatic_ligand_response_results_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/test_NSCLC_automatic_ligand_response_v1/launch_receipt.json').read_text());deadline=time.monotonic()+86400
    while not (O/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split()
        assert s[19]==str(launch['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(20)
    output=O/'private_independent_audit.json';assert not output.exists()
    terminal=json.loads((O/'private_terminal.json').read_text());plan=json.loads((O/'plan_before_any_ligand_response_effect.json').read_text());settings=plan['settings']
    assert settings==dict(sequential_h=100,maximum_null_draws=999999,null_batch=256,bootstrap_draws=2000,min_valid_bootstraps=1900,seed=2026101611)
    I=R/'NSCLC_automatic_ligand_response_preflight_20261016_v1';ap=I/'private_independent_audit.json';assert sha(ap)==plan['input_audit_sha256'];audit=json.loads(ap.read_text())
    scope_path=I/'private_all_unique_ligand_response_designs.tsv';assert sha(scope_path)==audit['unique_scope_sha256']==plan['unique_scope_sha256']
    scope=pd.read_csv(scope_path,sep='\t');scope_lookup=scope.set_index('record');result_path=O/'private_all_unique_ligand_response_results.tsv';assert sha(result_path)==terminal['result_sha256']
    t=pd.read_csv(result_path,sep='\t');assert len(t)==len(scope)==plan['complete_unique_hypothesis_family'] and t.record.tolist()==scope.record.tolist()
    llpath=O/'private_all_estimable_LOO.tsv.gz';assert sha(llpath)==terminal['LOO_sha256'];loo=pd.read_csv(llpath,sep='\t')
    loo_groups={name:f.set_index('removed_case_index') for name,f in loo.groupby('record')}
    hashes_path=I/'private_all_ready_design_SHA256.tsv';assert sha(hashes_path)==audit['ready_design_hash_table_sha256'];hashes=pd.read_csv(hashes_path,sep='\t').set_index('record').sha256
    observed=t.sequential_wild_p.notna();p=t.sequential_wild_p.fillna(1.).to_numpy(float)
    close(t.loc[observed,'BH_q_complete_unique_family'],manual_q(p)[observed]);close(t.loc[observed,'BY_q_dependence_sensitivity'],manual_q(p,True)[observed])
    assert t.loc[~observed,'BH_q_complete_unique_family'].isna().all() and t.loc[~observed,'BY_q_dependence_sensitivity'].isna().all()
    estimates=t[observed];selected=[]
    if len(estimates):
        selected=[estimates.record.iloc[0],estimates.record.iloc[len(estimates)//2]]
        cap=estimates[estimates.stopping_reason.eq('FIXED_MAXIMUM_STOP')]
        selected.append(cap.record.iloc[0] if len(cap) else estimates.record.iloc[-1]);selected=list(dict.fromkeys(selected))
    replays=[];points=0;LOOcount=0;maximum_rho_error=0.
    with threadpool_limits(limits=1):
        for row in t.itertuples():
            parent=scope_lookup.loc[row.record]
            if parent.status!='READY_FOR_INDEPENDENT_DESIGN_AUDIT':assert row.status==parent.status and pd.isna(row.sequential_wild_p);continue
            path=I/'private_frozen_designs'/(row.record+'.npz');assert sha(path)==hashes.loc[row.record];z=np.load(path,allow_pickle=False)
            x=z['X'];y=z['Y'];Z=z['Z'];n=len(x);xr=rankdata(x);yr=rankdata(y);xr=(xr-xr.mean())/xr.std();yr=(yr-yr.mean())/yr.std();D=np.column_stack([Z,xr])
            fit=sm.OLS(yr,D).fit(cov_type='HC3',use_t=True)
            if pd.isna(row.sequential_wild_p):
                assert row.status=='NO_NONZERO_ERROR_VARIANCE_FOR_STUDENTIZED_TEST' and (not np.isfinite(fit.bse[-1]) or fit.bse[-1]<=1e-12);continue
            assert row.status=='ESTIMABLE_EXPLORATORY_LIGAND_TO_GENE_DISJOINT_RECEIVER_PROGRAMME' and row.cases==n
            observed_rho=rho(x,y,Z);close(observed_rho,row.partial_rank_rho);maximum_rho_error=max(maximum_rho_error,abs(observed_rho-row.partial_rank_rho))
            close([fit.params[-1],fit.bse[-1],fit.tvalues[-1]],[row.rank_beta,row.HC3_se,row.HC3_t])
            close(np.corrcoef(pd.Series(x).rank(),pd.Series(y).rank())[0,1],row.raw_spearman);close(rho(x,z['receiver_full_programme'],Z),row.full_receiver_signature_rho)
            h=int(settings['sequential_h']);N=int(settings['maximum_null_draws']);draws=int(row.null_draws);g=int(row.null_exceedances)
            if row.stopping_reason=='FIXED_H_EXCEEDANCE_STOP':assert g==h and h<=draws<=N;close(h/draws,row.sequential_wild_p)
            else:assert row.stopping_reason=='FIXED_MAXIMUM_STOP' and draws==N and g<h;close((g+1)/(N+1),row.sequential_wild_p)
            assert 0<row.sequential_wild_p<=1
            for name,extra in [('sender_programme',z['sender_programme']),('pathology',z['response'])]:
                extra=np.asarray(extra,float);sd=extra.std();zz=np.column_stack([Z,(extra-extra.mean())/sd]) if sd>1e-12 else None
                ok=zz is not None and np.linalg.matrix_rank(zz)==20 and np.linalg.matrix_rank(np.column_stack([zz,rankdata(x)]))==21
                status=getattr(row,name+'_sensitivity_status')
                assert status==('DESCRIPTIVE_FULL_ADDITIONAL_CONTROL' if ok else 'ADDITIONAL_CONTROL_NOT_IDENTIFIABLE')
                close(rho(x,y,zz) if ok else np.nan,getattr(row,name+'_adjusted_rho'))
            ll=loo_groups[row.record];assert ll.index.tolist()==list(range(n));values=[]
            for k in range(n):
                use=np.arange(n)!=k;zz=Z[use];xx=x[use]
                valid=np.linalg.matrix_rank(zz)==19 and np.linalg.matrix_rank(np.column_stack([zz,rankdata(xx)]))==20
                value=rho(xx,y[use],zz) if valid else np.nan;close(value,ll.loc[k,'rho']);LOOcount+=1
                if np.isfinite(value):values.append(value)
            assert row.LOO_estimable==len(values) and row.LOO_same_sign==int(np.sum(np.sign(values)==np.sign(observed_rho)))
            close([min(values) if values else np.nan,max(values) if values else np.nan],[row.LOO_rho_minimum,row.LOO_rho_maximum])
            if row.record in selected:
                seed=np.random.SeedSequence([int(settings['seed']),int(hashlib.sha256(row.record.encode()).hexdigest()[:12],16)])
                rng,brng=[np.random.default_rng(s) for s in seed.spawn(2)]
                G=np.linalg.pinv(D);H=np.sum(D*G.T,axis=1);G0=np.linalg.pinv(Z);H0=np.sum(Z*G0.T,axis=1);fit0=Z@(G0@yr);e0=(yr-fit0)/np.sqrt(1-H0)
                consumed=0;exceed=0;pv=None
                while pv is None:
                    b=min(int(settings['null_batch']),N-consumed);sign=rng.choice(np.array([-1.,1.]),size=(n,b));star=fit0[:,None]+e0[:,None]*sign
                    coef=G@star;res=star-D@coef;se=np.sqrt(np.sum((G[-1,:,None]*res/(1-H[:,None]))**2,axis=0));null=np.divide(coef[-1],se,out=np.full(b,np.inf),where=se>0)
                    which=np.flatnonzero(np.abs(null)>=abs(row.HC3_t)-1e-12);needed=h-exceed
                    if len(which)>=needed:
                        consumed+=int(which[needed-1])+1;exceed=h;pv=h/consumed;reason='FIXED_H_EXCEEDANCE_STOP'
                    else:
                        consumed+=b;exceed+=len(which)
                        if consumed==N:pv=(exceed+1)/(N+1);reason='FIXED_MAXIMUM_STOP'
                assert consumed==draws and exceed==g and reason==row.stopping_reason;close(pv,row.sequential_wild_p)
                boots=[]
                for _ in range(int(settings['bootstrap_draws'])):
                    take=brng.integers(0,n,n);zz=Z[take]
                    if np.linalg.matrix_rank(zz)!=19 or np.linalg.matrix_rank(np.column_stack([zz,rankdata(x[take])]))!=20:continue
                    value=rho(x[take],y[take],zz)
                    if np.isfinite(value):boots.append(value)
                assert len(boots)==row.valid_bootstraps
                ci=np.quantile(boots,[.025,.975]) if len(boots)>=int(settings['min_valid_bootstraps']) else [np.nan,np.nan];close(ci,[row.CI_low,row.CI_high])
                replays.append(dict(record=row.record,null_draws=consumed,stopping_reason=reason,bootstrap_draws=int(settings['bootstrap_draws']),valid_bootstraps=len(boots)))
            points+=1
            if points%50==0:(O/'private_audit_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),points=points,all_LOO_records=LOOcount)))
    assert LOOcount==len(loo)
    result=dict(status='PASS_COMPLETE_UNIQUE_FAMILY_ALL_HC3_POINTS_SENSITIVITIES_BH_BY_LOO_AND_PRESPECIFIED_FULL_SEQUENTIAL_BOOTSTRAP_REPLAYS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(O/'private_terminal.json'),
        result_sha256=sha(result_path),complete_hypothesis_family=len(t),all_estimable_points=points,all_LOO_records=LOOcount,maximum_partial_rho_error=maximum_rho_error,
        complete_random_stream_replays=replays,all_random_streams_independently_replayed=False,directed_receptor_specific_or_causal_interaction_validated=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
