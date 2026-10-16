"""Independent OLS reconstruction of every functional description and unavailable slot."""
import hashlib,json,os,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy.linalg import svdvals
from scipy.stats import rankdata

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
QUALITY=['fraction','logN','mean_logRNA','mean_loggenes']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def full(d,rank):
    s=svdvals(d)
    return int((s>s[0]*1e-12).sum())==rank==np.linalg.matrix_rank(d) and len(d)-rank>=10 and s[0]/s[-1]<1e8


def rho(x,y,z=None):
    if not np.isfinite(x).all() or not np.isfinite(y).all() or np.ptp(x)<=1e-12 or np.ptp(y)<=1e-12:return np.nan
    a=rankdata(x).astype(float);b=rankdata(y).astype(float);a-=a.mean();b-=b.mean();ssx=a@a;ssy=b@b
    if z is not None:a=sm.OLS(a,z).fit().resid;b=sm.OLS(b,z).fit().resid
    if a@a<=1e-16*ssx or b@b<=1e-16*ssy:return np.nan
    return float(np.corrcoef(a,b)[0,1])


def close(a,b):assert np.allclose(a,b,atol=1e-8,rtol=1e-8,equal_nan=True)


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'NSCLC_disjoint_hallmark_functional_descriptions_20261016_v1';job=r/'HNSCC_receiver_followup_job_control_20261016_v2/describe_NSCLC_disjoint_hallmark_v1'
    receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text())
    for n,h in terminal['files'].items():assert sha(out/n)==h
    plan=json.loads((out/'plan_before_any_functional_association_description.json').read_text());assert plan['code_sha256']==receipt['source_sha256']
    measurement=r/'NSCLC_disjoint_hallmark_measurement_20261016_v1';proof=measurement/'private_independent_measurement_audit.json';assert sha(proof)==plan['measurement_audit_sha256'] and json.loads(proof.read_text())['status'].startswith('PASS_ALL1200')
    pre=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1';scope=pd.read_csv(pre/'private_all10600_description_scope.tsv',sep='\t');t=pd.read_csv(out/'private_all10600_functional_descriptions.tsv',sep='\t')
    assert len(t)==len(scope)==10600 and t.description_id.tolist()==scope.description_id.tolist()
    pd.testing.assert_frame_equal(t[scope.columns],scope,check_dtype=False)
    coverage=pd.read_csv(pre/'private_all1200_functional_scope.tsv',sep='\t').set_index(['signature','reference'])
    scores=pd.read_csv(r/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz',sep='\t')
    c=pd.read_csv(r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',sep='\t').set_index('sampleID')
    v=np.load(measurement/'private_all_native_functional_scores_and_QA.npz',allow_pickle=False);arrays={k:v[k] for k in ['ULM_score','target_mean','nontarget_mean','expressed_target_genes']}
    checked=0;maximum=0.;positive=0;negative=0
    for rec,f in t.groupby('original_clinical_record',sort=True):
        first=f.iloc[0];s=scores[scores.context.eq(first.context)&scores.type.eq(first.type)&scores.program.eq(first.program)&scores.sub_cell_type.eq(first.sub_cell_type)]
        cc=c.reindex(s.sampleID);use=s.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool);s=s.loc[use];cc=cc.loc[use];n=len(s);ids=s.index.to_numpy(int);x=s.raw_weighted_score.to_numpy(float);response=cc.response.to_numpy(float)
        cov=np.column_stack([cc[COMMON].to_numpy(float),s[QUALITY].to_numpy(float)]);counts=np.bincount(response.astype(int),minlength=2);status='DESCRIBABLE_NO_NEW_P';rank=np.nan;condition=np.nan;hmax=np.nan;Z=None;ZR=None;zclinical=None
        if n<50 or counts.min()<10:status='INSUFFICIENT_ORIGINAL_CASE_OR_RESPONSE_GROUP_COVERAGE'
        elif not np.isfinite(cov).all() or np.any(cov.std(0)<=1e-12):status='CONSTANT_OR_NONFINITE_FIXED_CONTROL'
        elif np.ptp(x)<=1e-12:status='CONSTANT_ORIGINAL_PROGRAMME'
        else:
            Z=np.column_stack([np.ones(n),(cov-cov.mean(0))/cov.std(0)]);rx=rankdata(x);rx=(rx-rx.mean())/rx.std();D=np.column_stack([Z,rx]);rank=int(np.linalg.matrix_rank(D));condition=float(np.linalg.cond(D))
            if not full(D,14):status='NONIDENTIFIABLE_COMPLETE_FUNCTIONAL_DESIGN'
            else:
                hmax=float(sm.OLS(x,D).fit().get_influence().hat_matrix_diag.max());candidate=np.column_stack([Z,response])
                if full(candidate,14):
                    zclinical=candidate
                    if full(np.column_stack([candidate,rx]),15):ZR=candidate
        for row in f.itertuples():
            assert row.eligible_cases==n==row.old_eligible_cases and row.response0==counts[0] and row.response1==counts[1]
            close([row.design_rank,row.condition,row.maximum_leverage],[rank,condition,hmax]);assert bool(row.extra_response_control_available)==(ZR is not None)
            signature=f'{row.context}_{row.type}_{row.program}';cv=coverage.loc[(signature,row.reference)];j=int(cv.reference_column);y=arrays['ULM_score'][ids,j];fg=arrays['target_mean'][ids,j];bg=arrays['nontarget_mean'][ids,j];det=arrays['expressed_target_genes'][ids,j]
            assert row.native_target_genes==cv.native_measured_targets and np.isclose(row.reference_retained_fraction,cv.original_reference_fraction_retained)
            assert row.positive_foreground_cases==int((fg>0).sum()) and row.minimum_expressed_targets==(int(det.min()) if n else -1)
            close(row.median_expressed_targets,float(np.median(det)) if n else np.nan)
            assert bool(row.foreground_expression_variable)==bool(n>1 and np.isfinite(fg).all() and np.ptp(fg)>1e-12)
            expected_status=status;main=np.nan;raw=np.nan;extra=np.nan;foreground=np.nan;background=np.nan;clinical=np.nan
            if cv.status!='MEASURABLE_INPUT_NO_RNA_SCORE':expected_status='FUNCTIONAL_INPUT_COVERAGE_UNAVAILABLE'
            elif status=='DESCRIBABLE_NO_NEW_P':
                if not np.isfinite(y).all():expected_status='NONFINITE_READOUT_ORIGINAL_CASES_NOT_REMOVED'
                elif not np.any(fg>0):expected_status='BACKGROUND_ONLY_SCORE_ALL_FOREGROUND_RNA_ZERO'
                elif np.ptp(y)<=1e-12:expected_status='CONSTANT_HELDOUT_RNA_READOUT'
                else:
                    main=rho(x,y,Z)
                    if not np.isfinite(main):expected_status='NO_READOUT_VARIATION_AFTER_FULL_CONTROLS'
                    else:
                        raw=rho(x,y);extra=rho(x,y,ZR) if ZR is not None else np.nan;foreground=rho(x,fg,Z);background=rho(x,bg,Z)
                        if zclinical is not None:clinical=float(sm.OLS((y-y.mean())/y.std(ddof=1),zclinical).fit().params[-1])
            assert row.status==expected_status
            got=[row.partial_rank_rho,row.raw_spearman,row.response_controlled_rho,row.target_mean_partial_rho,row.nontarget_mean_partial_rho,row.functional_MPR_minus_nonMPR_SD_beta]
            expected=[main,raw,extra,foreground,background,clinical];close(got,expected)
            assert bool(row.clinical_direction_available)==bool(np.isfinite(clinical))
            if np.isfinite(main):maximum=max(maximum,abs(row.partial_rank_rho-main));positive+=int(main>0);negative+=int(main<0)
            checked+=1
        (out/'private_independent_audit_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),description_rows_checked=checked))+'\n')
    summary=pd.read_csv(out/'private_all1200_author_subtype_description_summary.tsv',sep='\t');assert len(summary)==1200
    for row in summary.itertuples():
        f=t[t.signature.eq(row.signature)&t.reference.eq(row.reference)];a=f.partial_rank_rho.dropna()
        assert row.planned_author_subtypes==len(f) and row.describable_author_subtypes==len(a)
        close([row.subtype_median_rho,row.minimum_rho,row.maximum_rho],[np.median(a) if len(a) else np.nan,np.min(a) if len(a) else np.nan,np.max(a) if len(a) else np.nan])
    assert checked==10600 and positive>0 and negative>0 and t.partial_rank_rho.isna().any()
    result=dict(status='PASS_ALL10600_CASE_CONTROL_KEYS_STATUSES_PRIMARY_AND_ALL_SENSITIVITY_OLS_POINTS_ALL1200_SUMMARIES_NO_NEW_P',finished_UTC=datetime.now(timezone.utc).isoformat(),
                audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),description_rows_checked=checked,maximum_main_rho_difference=maximum,positive_negative_and_unavailable_cases_checked=True,new_association_P=False)
    p=out/'private_independent_description_audit.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
