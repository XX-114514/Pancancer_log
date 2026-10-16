"""All clinical designs/HC3/BH/LOO, three full independent randomization replays."""
import hashlib,json,os,time,re
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.linalg import qr,solve_triangular,svdvals
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from threadpoolctl import threadpool_limits


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def close(a,b):assert np.allclose(a,b,atol=1e-7,rtol=1e-8,equal_nan=True),(a,b)


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);I=R/'NSCLC_frozen_signature_clinical_preflight_20261016_v1';O=R/'NSCLC_frozen_signature_clinical_results_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/test_NSCLC_frozen_signature_clinical_v1/launch_receipt.json').read_text());deadline=time.monotonic()+7200
    while not (O/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists()
        s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(launch['proc_start_ticks']) and s[0]!='Z'
        assert time.monotonic()<deadline;time.sleep(10)
    out=O/'private_independent_audit.json';assert not out.exists()
    producer=json.loads((O/'private_terminal.json').read_text());path=O/'private_all_planned_signature_clinical_results.tsv';assert sha(path)==producer['result_sha256']
    plan=json.loads((I/'plan_before_any_clinical_effect.json').read_text());assert sha(I/'plan_before_any_clinical_effect.json')==producer['preflight_contract_sha256']
    t=pd.read_csv(path,sep='\t');planned=pd.read_csv(I/'private_all_planned_signature_subtype_designs.tsv',sep='\t')
    assert t.record.tolist()==planned.record.tolist() and len(t)==producer['planned_family']
    q=multipletests(t.wild_bootstrap_p.fillna(1),method='fdr_bh')[1];close(t.loc[t.wild_bootstrap_p.notna(),'BH_q_full_planned_family'],q[t.wild_bootstrap_p.notna()])
    loo=pd.read_csv(O/'private_all_estimable_clinical_LOO.tsv',sep='\t')
    c=pd.read_csv(I/'private_all243_clinical_eligibility.tsv',sep='\t').set_index('sampleID')
    original=pd.read_csv(R/'NSCLC_all_case_clinical_scope_20261016_v1/private_all243_author_case_fields.tsv',sep='\t',keep_default_na=False).set_index('sampleID')
    assert c.index.tolist()==original.index.tolist()
    expected_pd1=original['anti-PD1_therapy'].str.contains('Tislelizumab|Pembrolizumab|Camrelizumab|Sintilimab|Nivolumab|Toripalimab')|original['anti-PD1_therapy'].eq('Yes')
    expected_pd1&=~original['anti-PD1_therapy'].str.contains('Placebo');assert np.array_equal(c.PD1_recorded,expected_pd1)
    close(c.response,original.pathological_response.map({'pCR':1.,'MPR':1.,'non-MPR':0.}))
    for col,expected in [('age_numeric',pd.to_numeric(original.age,errors='coerce')),('male',original.gender.map({'M':1.,'F':0.})),('smoker',original.smoking_history.map({'Y':1.,'N':0.})),
        ('LUSC',original.cancer_type.map({'LUSC':1.,'LUAD':0.})),('cycles_numeric',pd.to_numeric(original.cycles,errors='coerce'))]:close(c[col],expected)
    stage=original.pre_treatment_staging.str.extract(r'^(IV|III|II|I)',expand=False).map({'I':1.,'II':2.,'III':3.,'IV':4.});close(c.stage_ordinal,stage)
    chemo=pd.Series(np.nan,index=original.index);chemo.loc[original.chemotherapy.eq('No')]=0.
    chemo.loc[original.chemotherapy.str.contains('Carboplatin|Cisplatin|Nedaplatin|Lobaplatin|Abraxane|Pemetrexed|Paclitaxel|Gemcitabine|Docetaxel')|original.chemotherapy.eq('Yes')]=1.;close(c.chemo_binary,chemo)
    targeted=pd.Series(np.nan,index=original.index);targeted.loc[original.targeted_therapy.eq('No')]=0.
    targeted.loc[original.targeted_therapy.str.contains('Apatinib|Bevacizumab|Avastin|Endostar|Anlotinib')]=1.;close(c.targeted_binary,targeted)
    controls=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
    expected_eligible=expected_pd1&c[controls+['response']].notna().all(axis=1);assert np.array_equal(c.clinical_eligible,expected_eligible)
    scorepath=R/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz';assert sha(scorepath)==plan['source_score_sha256']
    scores=pd.read_csv(scorepath,sep='\t');ready=t[t.wild_bootstrap_p.notna()];check_keys=set(ready.record.iloc[np.unique(np.array([0,len(ready)//2,len(ready)-1]))])
    maximum=0.;LOOcount=0;replays=[]
    with threadpool_limits(limits=1):
        for row in t.itertuples():
            f=scores[scores.context.eq(row.context)&scores.type.eq(row.type)&scores.program.eq(row.program)&scores.sub_cell_type.eq(row.sub_cell_type)]
            cc=c.reindex(f.sampleID);mask=f.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool);f=f.loc[mask];cc=cc.loc[mask]
            assert len(f)==row.eligible_cases and f.sampleID.is_unique
            if pd.isna(row.wild_bootstrap_p):assert row.status!='ESTIMABLE_POSTTREATMENT_SIGNATURE_CASE_ASSOCIATION' and pd.isna(row.BH_q_full_planned_family);continue
            z=np.load(I/'private_frozen_designs'/(row.record+'.npz'),allow_pickle=False);D=z['D'];Z=z['Z'];Y=z['Y'];n=len(Y)
            assert z['case_blocks'].tolist()==f.sampleID.tolist() and np.array_equal(z['X'],cc.response)
            cov=np.column_stack([cc[controls].to_numpy(),f[['fraction','logN','mean_logRNA','mean_loggenes']].to_numpy()])
            expected=np.column_stack([np.ones(n),cc.response,(cov-cov.mean(0))/cov.std(0)]);close(D,expected);close(Z,np.delete(D,1,axis=1))
            y=f.raw_weighted_score.to_numpy();close(Y,(y-y.mean())/y.std(ddof=1))
            fit=sm.OLS(Y,D).fit(cov_type='HC3',use_t=True);close([row.beta_MPR_minus_nonMPR,row.HC3_se,row.HC3_t],[fit.params[1],fit.bse[1],fit.tvalues[1]])
            maximum=max(maximum,abs(row.HC3_t-fit.tvalues[1]))
            old=(f.score.to_numpy()-f.score.mean())/y.std(ddof=1);close(old,z['Y_original_frozen_residual']);ff=sm.OLS(old,D).fit(cov_type='HC3',use_t=True)
            close([row.original_residual_score_beta,row.original_residual_score_t],[ff.params[1],ff.tvalues[1]]);close(ff.tvalues[1],fit.tvalues[1])
            ranks=rankdata(Y);n1=int(cc.response.sum());n0=n-n1;u=ranks[cc.response.to_numpy()==1].sum()-n1*(n1+1)/2
            close(row.unadjusted_rank_biserial,2*u/(n1*n0)-1)
            inv=np.linalg.pinv(D);h=np.sum(D*inv.T,axis=1);res=fit.resid
            expected_LOO=fit.params[1]-inv[1]*res/(1-h)
            ll=loo[loo.record.eq(row.record)].sort_values('removed_case_index');assert len(ll)==n
            finite=ll.beta.notna().to_numpy();close(ll.beta.to_numpy()[finite],expected_LOO[finite]);LOOcount+=n
            assert finite.sum()==row.LOO_estimable and np.sum(np.sign(expected_LOO[finite])==np.sign(fit.params[1]))==row.LOO_same_sign
            close(row.LOO_beta_minimum,np.min(expected_LOO[finite]));close(row.LOO_beta_maximum,np.max(expected_LOO[finite]))
            if row.record not in check_keys:continue
            seed=np.random.SeedSequence([int(plan['seed']),int(hashlib.sha256(row.record.encode()).hexdigest()[:12],16)])
            rng,brng=[np.random.default_rng(s) for s in seed.spawn(2)]
            Q,T=qr(D,mode='economic');G=solve_triangular(T,Q.T);hz=(Q*Q).sum(1)
            Q0,T0=qr(Z,mode='economic');fitted=Q0@(Q0.T@Y);error=(Y-fitted)/np.sqrt(1-(Q0*Q0).sum(1));exceed=0
            for start in range(0,int(plan['wild_bootstrap']),200):
                b=min(200,int(plan['wild_bootstrap'])-start);sign=rng.choice(np.array([-1.,1.]),size=(n,b));star=fitted[:,None]+error[:,None]*sign
                coef=solve_triangular(T,Q.T@star);e=star-Q@(Q.T@star);se=np.sqrt(np.sum((G[1,:,None]*e/(1-hz[:,None]))**2,axis=0))
                null=coef[1]/se;exceed+=int(np.sum(np.abs(null)>=abs(fit.tvalues[1])-1e-12))
            pv=(1+exceed)/(int(plan['wild_bootstrap'])+1);assert pv==row.wild_bootstrap_p or abs(pv-row.wild_bootstrap_p)<1e-14
            boots=[]
            for _ in range(int(plan['case_bootstrap'])):
                take=brng.integers(0,n,n);a=D[take];sv=svdvals(a);rank=int(np.sum(sv>max(a.shape)*np.finfo(float).eps*sv[0]))
                if rank!=14:continue
                qq,tt=qr(a,mode='economic');boots.append(float(solve_triangular(tt,qq.T@Y[take])[1]))
            assert len(boots)==row.valid_bootstraps
            ci=np.quantile(boots,[.025,.975]) if len(boots)>=int(plan['minimum_valid_bootstrap']) else [np.nan,np.nan];close(ci,[row.CI_low,row.CI_high])
            replays.append(dict(record=row.record,wild_draws=int(plan['wild_bootstrap']),case_draws=int(plan['case_bootstrap']),valid_cases=len(boots)))
    assert len(replays)==len(check_keys)
    a=dict(status='PASS_ALL_PLANNED_CLINICAL_SCOPE_EXPOSURE_DESIGNS_HC3_BH_LOO_AND_THREE_FULL_WILD_BOOTSTRAP_REPLAYS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        result_sha256=sha(path),audit_code_sha256=sha(__file__),planned_models=len(t),all_estimable_points=len(ready),all_LOO_records=LOOcount,
        maximum_HC3_t_error=maximum,full_randomization_replay_cases=replays,all_models_randomization_replayed=False,
        raw_residual_adjusted_inference_invariance_verified=True,patient_identity_linkage_reviewed=False,clinical_prediction_or_LR_mechanism_verified=False)
    out.write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a))


if __name__=='__main__':main()
