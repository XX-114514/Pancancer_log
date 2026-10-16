"""Full prespecified functional descriptions, no new association P or BH."""
import hashlib,json,os,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.stats import rankdata

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
QUALITY=['fraction','logN','mean_logRNA','mean_loggenes']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def correlation(x,y,q=None):
    if not np.isfinite(x).all() or not np.isfinite(y).all() or np.ptp(x)<=1e-12 or np.ptp(y)<=1e-12:return np.nan
    a=rankdata(x).astype(float);b=rankdata(y).astype(float);a-=a.mean();b-=b.mean();sa=a@a;sb=b@b
    if q is not None:a-=q@(q.T@a);b-=q@(q.T@b)
    if a@a<=1e-16*sa or b@b<=1e-16*sb:return np.nan
    return float(a@b/(np.linalg.norm(a)*np.linalg.norm(b)))


def identifiable(d,expected_rank):
    s=np.linalg.svd(d,compute_uv=False);rank=int((s>s[0]*1e-12).sum())
    return rank==expected_rank==np.linalg.matrix_rank(d) and len(d)-rank>=10 and s[0]/s[-1]<1e8


def synthetic_check():
    import statsmodels.api as sm
    rng=np.random.default_rng(2026101615);n=73;c=rng.normal(size=(n,12));z=np.column_stack([np.ones(n),c]);q=np.linalg.qr(z,mode='reduced')[0]
    x=rng.normal(size=n);x[:8]=0;response=rng.binomial(1,.5,n);maximum=0.
    for sign in [-1,1]:
        y=sign*x+c[:,0]+rng.normal(size=n)
        expected=np.corrcoef(sm.OLS(rankdata(x),z).fit().resid,sm.OLS(rankdata(y),z).fit().resid)[0,1]
        got=correlation(x,y,q);assert np.isclose(got,expected,atol=1e-10);maximum=max(maximum,abs(got-expected))
        ys=(y-y.mean())/y.std(ddof=1);rx=response-q@(q.T@response);b=rx@ys/(rx@rx)
        assert np.isclose(b,sm.OLS(ys,np.column_stack([z,response])).fit().params[-1],atol=1e-10)
    assert np.isnan(correlation(x,np.ones(n),q));assert not identifiable(np.column_stack([z,c[:,0]]),14)
    return dict(status='PASS_SYNTHETIC_RANK_RESIDUAL_AND_BINARY_DIRECTION_OLS_NOT_REAL_EFFECTS',maximum_rho_difference=maximum)


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);measurement=r/'NSCLC_disjoint_hallmark_measurement_20261016_v1';proof=measurement/'private_independent_measurement_audit.json'
    receipt=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_disjoint_hallmark_measurement_v1/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not proof.exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    assert json.loads(proof.read_text())['status'].startswith('PASS_ALL1200')
    terminal=json.loads((measurement/'private_terminal.json').read_text());assert sha(measurement/'private_terminal.json')==json.loads(proof.read_text())['producer_terminal_sha256']
    for n,h in terminal['files'].items():assert sha(measurement/n)==h
    pre=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1';before=json.loads((pre/'plan_before_functional_RNA_measurements.json').read_text())
    score_path=r/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz';clinical_path=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv'
    assert sha(score_path)==before['input_sha256']['scores'] and sha(clinical_path)==before['input_sha256']['clinical']
    out=r/'NSCLC_disjoint_hallmark_functional_descriptions_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),measurement_audit_sha256=sha(proof),input_scope_plan_sha256=sha(pre/'plan_before_functional_RNA_measurements.json'),
              synthetic_check=synthetic_check(),scope='All10600 frozen author-subtype/signature/reference positions, no old clinical/LR/function significance filter.',
              main='Rank raw frozen programme and heldout ULM; residual correlation under all8 clinical and4 own composition/quality controls plus intercept. >=50 cases,10 each old pathology group,14-column fullrank/condition<1e8/df10.',
              descriptions='Unadjusted rank correlation; extra-response-controlled rank; target/non-target means versus programme rank; raw RNA support; standardized functional-score MPR coefficient under same12 controls (no P or CI).',
              unavailable='No removal of cases to fit a readout. Missing/constant/nonidentifiable marked; all foreground RNA zero marked background-only rather than a functional response.',
              no_new_P_BH_or_cell_state_or_independent_validation=True)
    (out/'plan_before_any_functional_association_description.json').write_text(json.dumps(plan,indent=2)+'\n')
    scope=pd.read_csv(pre/'private_all10600_description_scope.tsv',sep='\t');coverage=pd.read_csv(pre/'private_all1200_functional_scope.tsv',sep='\t').set_index(['signature','reference'])
    scores=pd.read_csv(score_path,sep='\t');c=pd.read_csv(clinical_path,sep='\t').set_index('sampleID');z=np.load(measurement/'private_all_native_functional_scores_and_QA.npz',allow_pickle=False)
    arrays={k:z[k] for k in ['ULM_score','target_mean','nontarget_mean','expressed_target_genes']};rows=[]
    for clinical_record,g in scope.groupby('original_clinical_record',sort=True):
        first=g.iloc[0];f=scores[scores.context.eq(first.context)&scores.type.eq(first.type)&scores.program.eq(first.program)&scores.sub_cell_type.eq(first.sub_cell_type)]
        cc=c.reindex(f.sampleID);use=f.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool);f=f.loc[use];cc=cc.loc[use];n=len(f)
        assert n==int(first.old_eligible_cases) and f.sampleID.is_unique and np.array_equal(f.sampleID,cc.index)
        ids=f.index.to_numpy(int);x=f.raw_weighted_score.to_numpy(float);response=cc.response.to_numpy(float);counts=np.bincount(response.astype(int),minlength=2)
        cov=np.column_stack([cc[COMMON].to_numpy(float),f[QUALITY].to_numpy(float)]);status='DESCRIBABLE_NO_NEW_P';rank=np.nan;condition=np.nan;hmax=np.nan;q=None;qresponse=None;clinical_residual=None
        if n<50 or counts.min()<10:status='INSUFFICIENT_ORIGINAL_CASE_OR_RESPONSE_GROUP_COVERAGE'
        elif not np.isfinite(cov).all() or np.any(cov.std(0)<=1e-12):status='CONSTANT_OR_NONFINITE_FIXED_CONTROL'
        elif np.ptp(x)<=1e-12:status='CONSTANT_ORIGINAL_PROGRAMME'
        else:
            Z=np.column_stack([np.ones(n),(cov-cov.mean(0))/cov.std(0)]);xx=rankdata(x);xx=(xx-xx.mean())/xx.std();D=np.column_stack([Z,xx])
            rank=int(np.linalg.matrix_rank(D));condition=float(np.linalg.cond(D))
            if not identifiable(D,14):status='NONIDENTIFIABLE_COMPLETE_FUNCTIONAL_DESIGN'
            else:
                q=np.linalg.qr(Z,mode='reduced')[0];hmax=float(np.max(np.sum(D*np.linalg.pinv(D).T,axis=1)));zr=np.column_stack([Z,response])
                if identifiable(zr,14):
                    clinical_residual=response-q@(q.T@response)
                    if identifiable(np.column_stack([zr,xx]),15):qresponse=np.linalg.qr(zr,mode='reduced')[0]
        for row in g.itertuples():
            a=row._asdict();signature=f'{row.context}_{row.type}_{row.program}';cv=coverage.loc[(signature,row.reference)];j=int(cv.reference_column)
            a.update(signature=signature,eligible_cases=n,response0=int(counts[0]),response1=int(counts[1]),design_rank=rank,condition=condition,maximum_leverage=hmax,status=status,
                     native_target_genes=int(cv.native_measured_targets),reference_retained_fraction=float(cv.original_reference_fraction_retained),partial_rank_rho=np.nan,raw_spearman=np.nan,
                     response_controlled_rho=np.nan,target_mean_partial_rho=np.nan,nontarget_mean_partial_rho=np.nan,functional_MPR_minus_nonMPR_SD_beta=np.nan,
                     extra_response_control_available=qresponse is not None,clinical_direction_available=False)
            y=arrays['ULM_score'][ids,j];foreground=arrays['target_mean'][ids,j];background=arrays['nontarget_mean'][ids,j];detected=arrays['expressed_target_genes'][ids,j]
            a.update(positive_foreground_cases=int(np.sum(foreground>0)),median_expressed_targets=float(np.median(detected)) if n else np.nan,
                     minimum_expressed_targets=int(detected.min()) if n else -1,foreground_expression_variable=bool(n>1 and np.isfinite(foreground).all() and np.ptp(foreground)>1e-12))
            if cv.status!='MEASURABLE_INPUT_NO_RNA_SCORE':a['status']='FUNCTIONAL_INPUT_COVERAGE_UNAVAILABLE'
            elif status=='DESCRIBABLE_NO_NEW_P':
                if not np.isfinite(y).all():a['status']='NONFINITE_READOUT_ORIGINAL_CASES_NOT_REMOVED'
                elif not np.any(foreground>0):a['status']='BACKGROUND_ONLY_SCORE_ALL_FOREGROUND_RNA_ZERO'
                elif np.ptp(y)<=1e-12:a['status']='CONSTANT_HELDOUT_RNA_READOUT'
                else:
                    a['partial_rank_rho']=correlation(x,y,q)
                    if not np.isfinite(a['partial_rank_rho']):a['status']='NO_READOUT_VARIATION_AFTER_FULL_CONTROLS'
                    else:
                        a['raw_spearman']=correlation(x,y);a['response_controlled_rho']=correlation(x,y,qresponse) if qresponse is not None else np.nan
                        a['target_mean_partial_rho']=correlation(x,foreground,q);a['nontarget_mean_partial_rho']=correlation(x,background,q)
                        if clinical_residual is not None and clinical_residual@clinical_residual>1e-12:
                            ys=(y-y.mean())/y.std(ddof=1);a['functional_MPR_minus_nonMPR_SD_beta']=float(clinical_residual@ys/(clinical_residual@clinical_residual));a['clinical_direction_available']=True
            rows.append(a)
        (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),description_rows=len(rows),original_clinical_slot=clinical_record))+'\n')
    table=pd.DataFrame(rows);assert len(table)==10600 and table.description_id.is_unique
    table.to_csv(out/'private_all10600_functional_descriptions.tsv',sep='\t',index=False,na_rep='NA')
    summary=[]
    for (signature,reference),f in table.groupby(['signature','reference'],sort=True):
        v=f.partial_rank_rho.dropna();summary.append(dict(signature=signature,reference=reference,planned_author_subtypes=len(f),describable_author_subtypes=len(v),
                                                        subtype_median_rho=float(v.median()) if len(v) else np.nan,minimum_rho=float(v.min()) if len(v) else np.nan,maximum_rho=float(v.max()) if len(v) else np.nan,
                                                        same_cases_across_subtypes_not_independent_replications=True))
    pd.DataFrame(summary).to_csv(out/'private_all1200_author_subtype_description_summary.tsv',sep='\t',index=False,na_rep='NA')
    terminal=dict(status='COMPLETE_ALL10600_TRAINING_DISJOINT_FUNCTION_DESCRIPTIONS_NO_NEW_P_INDEPENDENT_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
                  description_positions=len(table),status_counts=table.status.value_counts().to_dict(),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_P_BH_or_independent_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
