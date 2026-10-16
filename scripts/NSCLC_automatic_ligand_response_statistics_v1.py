"""Finite sequential wild inference and pre-data synthetic calibration.

The sequential Monte Carlo rule is not evidence that wild resampling is an
exact experiment. Its statistical interpretation retains wild-model assumptions.
"""
import hashlib,itertools,json,math,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.linalg import qr,solve_triangular
from scipy.stats import rankdata,binom
from threadpoolctl import threadpool_limits

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def consume(exceedances,previous_draws,previous_exceeds,h,maximum):
    """Stop at the FIRST h-th exceedance, not the end of its vector batch."""
    flags=np.asarray(exceedances,bool);remaining=maximum-previous_draws;flags=flags[:remaining]
    accumulated=previous_exceeds+np.cumsum(flags);stops=np.flatnonzero(accumulated>=h)
    if len(stops):
        draws=previous_draws+int(stops[0])+1
        return draws,h,h/draws,'FIXED_H_EXCEEDANCE_STOP'
    draws=previous_draws+len(flags);count=previous_exceeds+int(flags.sum())
    if draws==maximum:return draws,count,(count+1)/(maximum+1),'FIXED_MAXIMUM_STOP'
    return draws,count,None,None

def geometry(x,y,Z):
    xx=rankdata(x).astype(float);yy=rankdata(y).astype(float)
    if xx.std()<=1e-12 or yy.std()<=1e-12:return None
    xx=(xx-xx.mean())/xx.std();yy=(yy-yy.mean())/yy.std();D=np.column_stack([Z,xx])
    if np.linalg.matrix_rank(D)!=D.shape[1]:return None
    Q,T=qr(D,mode='economic');G=solve_triangular(T,Q.T);h=np.sum(Q*Q,axis=1)
    residual=yy-Q@(Q.T@yy);coefficient=G[-1]@yy;se=float(np.sqrt(np.sum((G[-1]*residual/(1-h))**2)))
    if not np.isfinite(se) or se<=1e-12:return None
    return dict(D=D,Q=Q,T=T,G=G,h=h,yy=yy,beta=float(coefficient),se=se,t=float(coefficient/se))

def partial_rho(x,y,Z):
    a=rankdata(x).astype(float);b=rankdata(y).astype(float);Q=qr(Z,mode='economic')[0]
    a-=Q@(Q.T@a);b-=Q@(Q.T@b);den=float(np.linalg.norm(a)*np.linalg.norm(b))
    return float(a@b/den) if den>1e-12 else np.nan

def sequential_wild(x,y,Z,rng,h=100,maximum=999999,batch=256):
    fit=geometry(x,y,Z)
    if fit is None:return None
    Q0=qr(Z,mode='economic')[0];yy=fit['yy'];f0=Q0@(Q0.T@yy);e0=(yy-f0)/np.sqrt(1-np.sum(Q0*Q0,axis=1))
    draws=0;count=0;p=None;reason=None
    while p is None:
        b=min(batch,maximum-draws);sign=rng.choice(np.array([-1.,1.]),size=(len(x),b));star=f0[:,None]+e0[:,None]*sign
        coef=solve_triangular(fit['T'],fit['Q'].T@star);res=star-fit['Q']@(fit['Q'].T@star)
        se=np.sqrt(np.sum((fit['G'][-1,:,None]*res/(1-fit['h'][:,None]))**2,axis=0))
        null=np.divide(coef[-1],se,out=np.full(b,np.inf),where=se>0)
        flags=np.abs(null)>=abs(fit['t'])-1e-12
        draws,count,p,reason=consume(flags,draws,count,h,maximum)
    return dict(HC3_t=fit['t'],rank_beta=fit['beta'],HC3_se=fit['se'],sequential_wild_p=float(p),null_draws=draws,
        null_exceedances=count,stopping_reason=reason,sequential_h=h,maximum_null_draws=maximum)

def model_statistics(record,path,expected_sha,settings):
    assert sha(path)==expected_sha;z=np.load(path,allow_pickle=False);x=z['X'];y=z['Y'];Z=z['Z'];n=len(x)
    seed=np.random.SeedSequence([int(settings['seed']),int(hashlib.sha256(record.encode()).hexdigest()[:12],16)])
    rng,brng=[np.random.default_rng(s) for s in seed.spawn(2)]
    result=sequential_wild(x,y,Z,rng,h=int(settings['sequential_h']),maximum=int(settings['maximum_null_draws']),batch=int(settings['null_batch']))
    if result is None:return dict(record=record,status='NO_NONZERO_ERROR_VARIANCE_FOR_STUDENTIZED_TEST'),[]
    rho=partial_rho(x,y,Z);boots=[]
    for _ in range(int(settings['bootstrap_draws'])):
        take=brng.integers(0,n,n);zz=Z[take]
        if np.linalg.matrix_rank(zz)!=19 or np.linalg.matrix_rank(np.column_stack([zz,rankdata(x[take])]))!=20:continue
        v=partial_rho(x[take],y[take],zz)
        if np.isfinite(v):boots.append(v)
    ci=np.quantile(boots,[.025,.975]) if len(boots)>=int(settings['min_valid_bootstraps']) else [np.nan,np.nan]
    loo=[];values=[]
    for k in range(n):
        take=np.arange(n)!=k;zz=Z[take];xx=x[take]
        valid=np.linalg.matrix_rank(zz)==19 and np.linalg.matrix_rank(np.column_stack([zz,rankdata(xx)]))==20
        value=partial_rho(xx,y[take],zz) if valid else np.nan
        loo.append(dict(record=record,removed_case_index=k,rho=value))
        if np.isfinite(value):values.append(value)
    sensitivity={}
    for name,extra in [('sender_programme',z['sender_programme']),('pathology',z['response'])]:
        extra=np.asarray(extra,float);sd=extra.std();zz=np.column_stack([Z,(extra-extra.mean())/sd]) if sd>1e-12 else None
        ok=zz is not None and np.linalg.matrix_rank(zz)==20 and np.linalg.matrix_rank(np.column_stack([zz,rankdata(x)]))==21
        sensitivity[name+'_adjusted_rho']=partial_rho(x,y,zz) if ok else np.nan
        sensitivity[name+'_sensitivity_status']='DESCRIPTIVE_FULL_ADDITIONAL_CONTROL' if ok else 'ADDITIONAL_CONTROL_NOT_IDENTIFIABLE'
    result.update(record=record,status='ESTIMABLE_EXPLORATORY_LIGAND_TO_GENE_DISJOINT_RECEIVER_PROGRAMME',cases=n,partial_rank_rho=rho,CI_low=float(ci[0]),CI_high=float(ci[1]),
        valid_bootstraps=len(boots),raw_spearman=float(np.corrcoef(rankdata(x),rankdata(y))[0,1]),full_receiver_signature_rho=partial_rho(x,z['receiver_full_programme'],Z),
        LOO_estimable=len(values),LOO_same_sign=int(np.sum(np.sign(values)==np.sign(rho))),LOO_rho_minimum=min(values) if values else np.nan,LOO_rho_maximum=max(values) if values else np.nan,
        **sensitivity)
    return result,loo

def calibrate():
    import statsmodels.api as sm
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'NSCLC_ligand_response_statistics_calibration_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),statistics_source_sha256=sha(__file__),
        no_native_or_clinical_effects=True,exact_small_cap_enumeration='All2^9 Bernoulli paths, uniform true-tail mixing, h2/3, every support threshold, batch1/2/4/7.',
        real_inference_settings={'sequential_h':100,'maximum_null_draws':999999,'null_batch':256,'bootstrap_draws':2000,'min_valid_bootstraps':1900,'seed':2026101611},
        synthetic='500 null n96 full20-column designs, half explanatory ligands have50% zero ties; heterogeneous errors; 80 injected alternatives; h20 maximum4999 for calibration only.',
        gate='Exact finite-rule superuniformity; statsmodels HC3 parity; each250-null stratum rejects no more than one-sided99.5% binomial.05 bound; injected alternatives reject>=40/80.',
        limitations='Finite synthetic configurations diagnose gross errors, do not prove universal finite-sample wild validity, patient independence, causality, or BH assumptions.')
    (O/'plan_before_synthetic_calibration.json').write_text(json.dumps(plan,indent=2)+'\n')
    enumerated=[]
    for h in [2,3]:
        probabilities=[]
        for sequence in itertools.product([False,True],repeat=9):
            sequence=np.asarray(sequence,bool);g=int(sequence.sum());mass=1/(10*math.comb(9,g));reference=None
            for b in [1,2,4,7]:
                draws=0;exceed=0;p=None;reason=None
                while p is None:draws,exceed,p,reason=consume(sequence[draws:draws+b],draws,exceed,h,9)
                value=(draws,exceed,p,reason)
                if reference is None:reference=value
                else:assert value==reference
            probabilities.append((reference[2],mass))
        assert np.isclose(sum(mass for _,mass in probabilities),1.)
        error=max(sum(mass for p,mass in probabilities if p<=a+1e-14)-a for a in sorted({p for p,_ in probabilities}))
        assert error<1e-12;enumerated.append(dict(h=h,maximum=9,all_paths=len(probabilities),maximum_superuniform_error=error))
    rng=np.random.default_rng(2026101612);null=[];alternative=[];maximum_t_error=0.
    with threadpool_limits(limits=1):
        for k in range(580):
            n=96;cov=rng.normal(size=(n,18));Z=np.column_stack([np.ones(n),(cov-cov.mean(0))/cov.std(0)])
            x=rng.normal(size=n)
            tied=k%2==1
            if tied:x=np.where(rng.random(n)<.5,0.,np.exp(x/2))
            errors=rng.normal(size=n)*np.exp(.3*cov[:,0]);y=.3*cov[:,0]+.3*cov[:,1]+errors
            if k>=500:y+=.75*(x-x.mean())/x.std()
            fit=geometry(x,y,Z);assert fit is not None
            if k<8:
                sf=sm.OLS(fit['yy'],fit['D']).fit(cov_type='HC3',use_t=True);err=abs(float(sf.tvalues[-1])-fit['t']);assert err<1e-8;maximum_t_error=max(maximum_t_error,err)
            result=sequential_wild(x,y,Z,rng,h=20,maximum=4999,batch=256);assert result is not None
            row=dict(index=k,tied=tied,p=result['sequential_wild_p'],draws=result['null_draws'],stop=result['stopping_reason'])
            (null if k<500 else alternative).append(row)
    groups=[]
    for tied in [False,True]:
        f=[x for x in null if x['tied']==tied];reject=sum(x['p']<=.05 for x in f);bound=int(binom.ppf(.995,len(f),.05));assert reject<=bound
        groups.append(dict(tied=tied,replicates=len(f),reject05=reject,one_sided995pct_binomial_bound=bound))
    power=sum(x['p']<=.05 for x in alternative);assert power>=40
    (O/'private_synthetic_replications.json').write_text(json.dumps(dict(null=null,alternative=alternative),indent=2)+'\n')
    result=dict(status='PASS_EXACT_FINITE_SEQUENTIAL_RULE_BATCH_STOP_HC3_PARITY_AND_PRESPECIFIED_SYNTHETIC_NULL_ALTERNATIVE_CHECKS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),statistics_source_sha256=sha(__file__),finite_enumeration=enumerated,
        synthetic_null_strata=groups,alternative_replicates=len(alternative),alternative_reject05=power,maximum_HC3_t_error=maximum_t_error,
        all_configurations_or_scientific_effects_validated=False,real_inference_settings=plan['real_inference_settings'])
    (O/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':calibrate()
