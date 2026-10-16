"""All paired inputs, points, sensitivities, BH/LOO and three full replays."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from scipy.linalg import qr,solve_triangular,svdvals
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary','response']
UNIT=['fraction','logN','mean_logRNA','mean_loggenes']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def close(a,b):assert np.allclose(a,b,atol=1e-7,rtol=1e-8,equal_nan=True),(a,b)


def independent_rho(x,y,Z=None,ranks=True):
    a=rankdata(x).astype(float) if ranks else np.asarray(x,float)
    b=rankdata(y).astype(float) if ranks else np.asarray(y,float)
    if Z is not None:a=sm.OLS(a,Z).fit().resid;b=sm.OLS(b,Z).fit().resid
    else:a=a-a.mean();b=b-b.mean()
    return float(np.corrcoef(a,b)[0,1]) if np.linalg.norm(a)*np.linalg.norm(b)>1e-12 else np.nan


def full_rank(A):
    v=svdvals(A);return int(np.sum(v>max(A.shape)*np.finfo(float).eps*v[0]))


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);I=R/'NSCLC_two_end_coordination_preflight_20261016_v1';O=R/'NSCLC_two_end_coordination_results_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/test_NSCLC_two_end_coordination_v1/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (O/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists()
        st=p.read_text().rsplit(')',1)[1].split();assert st[19]==str(launch['proc_start_ticks']) and st[0]!='Z'
        assert time.monotonic()<deadline;time.sleep(15)
    output=O/'private_independent_audit.json';assert not output.exists()
    terminal=json.loads((O/'private_terminal.json').read_text());plan=json.loads((I/'plan_before_any_two_end_effect.json').read_text())
    path=O/'private_all572_two_end_results.tsv';assert sha(path)==terminal['result_sha256'];assert sha(I/'plan_before_any_two_end_effect.json')==terminal['parent_plan_sha256']
    t=pd.read_csv(path,sep='\t');scope=pd.read_csv(I/'private_all572_planned_pair_designs.tsv',sep='\t');assert len(t)==572 and t.record.tolist()==scope.record.tolist()
    loo=pd.read_csv(O/'private_all_estimable_pair_LOO.tsv',sep='\t')
    q=multipletests(t.wild_p.fillna(1.),method='fdr_bh')[1];close(t.loc[t.wild_p.notna(),'BH_q_full572'],q[t.wild_p.notna()])
    P=R/'NSCLC_context_native_programme_projection_20261016_v1';sp=P/'private_native_author_subtype_programme_scores.tsv.gz';assert sha(sp)==plan['input_sha256']['scores']
    scores=pd.read_csv(sp,sep='\t');cp=R/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv'
    assert sha(cp)==plan['input_sha256']['case_eligibility'];c=pd.read_csv(cp,sep='\t').set_index('sampleID')
    pp=R/'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv';assert sha(pp)==plan['input_sha256']['parameters'];params=pd.read_csv(pp,sep='\t')
    N=R/'NSCLC_full_native_gene_units_20261016_v1';genes=pd.Index(pd.read_csv(N/'private_native_gene_axis.tsv',sep='\t').gene);M=np.load(N/'mean_logCP10K.npy',mmap_mode='r')
    masks=pd.read_csv(I/'private_all16_parent_shared_gene_masks.tsv',sep='\t').set_index(['context','program','type']);masked={}
    for context,program in t[['context','program']].drop_duplicates().itertuples(index=False):
        a=params[params.context.eq(context)&params.program.eq(program)&params.type.eq('typeA')];b=params[params.context.eq(context)&params.program.eq(program)&params.type.eq('typeB')]
        shared=set(a.gene)&set(b.gene)
        for typ,p in [('typeA',a),('typeB',b)]:
            keep=~p.gene.isin(shared);provided=genes.get_indexer(p.gene)>=0;absw=np.abs(p.weight.to_numpy());retained=float(absw[keep].sum()/absw.sum())
            measured=float(absw[keep&provided].sum()/absw[keep].sum()) if keep.any() else 0.
            available=retained>=.8 and measured>=.95 and np.sum(keep&provided)>=20
            row=masks.loc[(context,program,typ)];close([retained,measured],[row.nonshared_original_L1_fraction,row.nonshared_native_weight_fraction]);assert row.nonshared_score_available==available
            if available:
                p=p.loc[keep&provided];f=scores[scores.context.eq(context)&scores.program.eq(program)&scores.type.eq(typ)]
                values=np.asarray(M[np.ix_(f.unit_index.to_numpy(int),genes.get_indexer(p.gene))])@(p.weight/p.gene_sd).to_numpy()-float((p.weight*p.gene_mean/p.gene_sd).sum())
                masked[(context,program,typ)]=pd.Series(values,index=f.unit_index)
    estimates=t[t.wild_p.notna()];selected=set(estimates.record.iloc[np.unique(np.array([0,len(estimates)//2,len(estimates)-1]))]);checks=[];maxerr=0.;LOOcount=0
    with threadpool_limits(limits=1):
        for row in t.itertuples():
            a=scores[scores.context.eq(row.context)&scores.program.eq(row.program)&scores.type.eq('typeA')&scores.sub_cell_type.eq(row.author_subtype_A)].set_index('sampleID')
            b=scores[scores.context.eq(row.context)&scores.program.eq(row.program)&scores.type.eq('typeB')&scores.sub_cell_type.eq(row.author_subtype_B)].set_index('sampleID')
            cases=sorted(set(a.index)&set(b.index));aa=a.loc[cases];bb=b.loc[cases];cc=c.loc[cases]
            use=aa.n_cells.ge(20).to_numpy()&bb.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool)
            aa=aa.loc[use];bb=bb.loc[use];cc=cc.loc[use];n=len(aa);assert n==row.eligible_paired_cases and len(cases)==row.shared_cases_any
            if pd.isna(row.wild_p):assert row.status!='ESTIMABLE_EXPLORATORY_CASE_PROGRAMME_COORDINATION' and pd.isna(row.BH_q_full572);continue
            z=np.load(I/'private_pair_designs'/(row.record+'.npz'),allow_pickle=False);x=z['X'];y=z['Y'];Z=z['Z'];assert z['cases'].tolist()==aa.index.tolist()==bb.index.tolist()
            close(x,aa.raw_weighted_score);close(y,bb.raw_weighted_score);close(z['oldA'],aa.score);close(z['oldB'],bb.score)
            cov=np.column_stack([cc[COMMON].to_numpy(float),aa[UNIT].to_numpy(float),bb[UNIT].to_numpy(float)])
            close(z['covariates'],cov);expected=np.column_stack([np.ones(n),(cov-cov.mean(0))/cov.std(0)]);close(Z,expected)
            xr=rankdata(x);yr=rankdata(y);xr=(xr-xr.mean())/xr.std();yr=(yr-yr.mean())/yr.std();D=np.column_stack([Z,xr]);assert full_rank(D)==19
            observed=independent_rho(x,y,Z);close(observed,row.partial_rank_rho);maxerr=max(maxerr,abs(observed-row.partial_rank_rho))
            fit=sm.OLS(yr,D).fit(cov_type='HC3',use_t=True);close(fit.tvalues[-1],row.HC3_t)
            close(independent_rho(x,y),row.raw_spearman);close(independent_rho(x,y,np.delete(Z,1+COMMON.index('response'),axis=1)),row.no_pathology_control_rho)
            close(independent_rho(x,y,Z,ranks=False),row.linear_signature_partial_rho);close(independent_rho(z['oldA'],z['oldB'],Z),row.old_residual_rank_partial_rho)
            ma=masked.get((row.context,row.program,'typeA'));mb=masked.get((row.context,row.program,'typeB'))
            if ma is not None and mb is not None:
                av=ma.reindex(aa.unit_index).to_numpy();bv=mb.reindex(bb.unit_index).to_numpy();close(av,z['maskedA']);close(bv,z['maskedB']);close(independent_rho(av,bv,Z),row.nonshared_gene_partial_rho)
            else:assert np.isnan(z['maskedA']).all() and np.isnan(z['maskedB']).all() and pd.isna(row.nonshared_gene_partial_rho)
            ll=loo[loo.record.eq(row.record)].set_index('removed_case_index');assert len(ll)==n;values=[]
            for k in range(n):
                take=np.arange(n)!=k;v=independent_rho(x[take],y[take],Z[take]) if full_rank(Z[take])==18 else np.nan
                close(v,ll.loc[k,'rho']);LOOcount+=1
                if np.isfinite(v):values.append(v)
            assert len(values)==row.LOO_estimable and np.sum(np.sign(values)==np.sign(observed))==row.LOO_same_sign
            close([min(values),max(values)],[row.LOO_rho_minimum,row.LOO_rho_maximum])
            if row.record in selected:
                seed=np.random.SeedSequence([int(plan['seed']),int(hashlib.sha256(row.record.encode()).hexdigest()[:12],16)])
                rng,brng=[np.random.default_rng(s) for s in seed.spawn(2)];Q,T=qr(D,mode='economic');G=solve_triangular(T,Q.T);h=np.sum(Q*Q,axis=1)
                Q0,_=qr(Z,mode='economic');fit0=Q0@(Q0.T@yr);e0=(yr-fit0)/np.sqrt(1-np.sum(Q0*Q0,axis=1));exceed=0
                for start in range(0,int(plan['wild_draws']),200):
                    b=min(200,int(plan['wild_draws'])-start);sign=rng.choice(np.array([-1.,1.]),size=(n,b));star=fit0[:,None]+e0[:,None]*sign
                    coef=solve_triangular(T,Q.T@star);res=star-Q@(Q.T@star);se=np.sqrt(np.sum((G[-1,:,None]*res/(1-h[:,None]))**2,axis=0));null=coef[-1]/se
                    exceed+=int(np.sum(np.abs(null)>=abs(fit.tvalues[-1])-1e-12))
                close((1+exceed)/(int(plan['wild_draws'])+1),row.wild_p)
                boot=[]
                for _ in range(int(plan['bootstrap_draws'])):
                    take=brng.integers(0,n,n);zz=Z[take]
                    if full_rank(zz)!=18 or full_rank(np.column_stack([zz,rankdata(x[take])]))!=19:continue
                    v=independent_rho(x[take],y[take],zz)
                    if np.isfinite(v):boot.append(v)
                assert len(boot)==row.valid_bootstraps
                ci=np.quantile(boot,[.025,.975]) if len(boot)>=int(plan['min_valid_bootstraps']) else [np.nan,np.nan];close(ci,[row.CI_low,row.CI_high])
                checks.append(dict(record=row.record,wild_draws=int(plan['wild_draws']),bootstrap_draws=int(plan['bootstrap_draws']),valid_bootstraps=len(boot)))
            (O/'private_audit_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),current_record=row.record,LOO_records=LOOcount)))
    assert len(checks)==len(selected)
    report=dict(status='PASS_ALL572_SCOPE_ALL_ESTIMABLE_INPUTS_HC3_SENSITIVITIES_BH_LOO_AND_THREE_FULL_RANDOMIZATION_REPLAYS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        result_sha256=sha(path),audit_code_sha256=sha(__file__),planned_models=len(t),estimable_points=len(estimates),all_LOO_records=LOOcount,maximum_partial_rho_error=maxerr,
        complete_randomization_replays=checks,all_models_randomization_independently_replayed=False,directed_or_causal_interaction_validated=False)
    output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
