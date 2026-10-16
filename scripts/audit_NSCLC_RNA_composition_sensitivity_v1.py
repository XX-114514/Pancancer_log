"""Independent complete OLS replay of author RNA-composition sensitivities."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm
from scipy.linalg import svdvals
from scipy.stats import rankdata

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
QUALITY=['fraction','logN','mean_logRNA','mean_loggenes'];QC=['mean_author_pct_mito','mean_author_pct_ribo']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def full(d,rank):
    s=svdvals(d)
    return int((s>s[0]*1e-12).sum())==rank==np.linalg.matrix_rank(d) and len(d)-rank>=10 and s[0]/s[-1]<1e8


def rho(x,y,z):
    if not np.isfinite(x).all() or not np.isfinite(y).all() or np.ptp(x)<=1e-12 or np.ptp(y)<=1e-12:return np.nan
    a=rankdata(x).astype(float);b=rankdata(y).astype(float);a-=a.mean();b-=b.mean();ssx=a@a;ssy=b@b
    a=sm.OLS(a,z).fit().resid;b=sm.OLS(b,z).fit().resid
    if a@a<=1e-16*ssx or b@b<=1e-16*ssy:return np.nan
    return float(np.corrcoef(a,b)[0,1])


def make(cov,x):
    if not np.isfinite(cov).all() or np.any(cov.std(0)<=1e-12):return None,'CONSTANT_OR_NONFINITE_COMPLETE_CONTROL'
    z=np.column_stack([np.ones(len(x)),(cov-cov.mean(0))/cov.std(0)]);a=rankdata(x);a=(a-a.mean())/a.std();d=np.column_stack([z,a])
    return (z,'DESCRIBABLE_NO_NEW_P') if full(d,d.shape[1]) else (None,'NONIDENTIFIABLE_COMPLETE_DESIGN')


def beta(y,response,z):
    if z is None or not np.isfinite(y).all() or np.ptp(y)<=1e-12 or not full(np.column_stack([z,response]),z.shape[1]+1):return np.nan
    rr=sm.OLS(response,z).fit().resid
    if rr@rr<=1e-12:return np.nan
    return float(sm.OLS((y-y.mean())/y.std(ddof=1),np.column_stack([z,response])).fit().params[-1])


def close(got,expected):assert np.allclose(got,expected,rtol=1e-8,atol=1e-8,equal_nan=True)


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'NSCLC_RNA_composition_sensitivity_20261016_v1';job=r/'HNSCC_receiver_followup_job_control_20261016_v2/describe_NSCLC_RNA_composition_sensitivity_v1'
    receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text());plan=json.loads((out/'plan_before_any_QC_association.json').read_text());assert plan['code_sha256']==receipt['source_sha256']
    for f,h in terminal['files'].items():assert sha(out/f)==h
    paths=dict(qc=r/'NSCLC_author_RNA_composition_20261016_v1/private_all_native_unit_RNA_composition.tsv',old=r/'NSCLC_disjoint_hallmark_functional_descriptions_20261016_v1/private_all10600_functional_descriptions.tsv',scores=r/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz',clinical=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',readouts=r/'NSCLC_disjoint_hallmark_measurement_20261016_v1/private_all_native_functional_scores_and_QA.npz',references=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1/private_reference_axis.tsv')
    for key,p in paths.items():assert sha(p)==plan['input_sha256'][key]
    proof=r/'NSCLC_author_RNA_composition_20261016_v1/private_independent_QC_audit.json';assert sha(proof)==plan['qc_audit_sha256'] and json.loads(proof.read_text())['status'].startswith('PASS_ALL_SOURCE_CELL_ORDER')
    t=pd.read_csv(out/'private_all10600_functional_QC_sensitivities.tsv',sep='\t');slots=pd.read_csv(out/'private_all212_programme_composition_descriptions.tsv',sep='\t');old=pd.read_csv(paths['old'],sep='\t');scores=pd.read_csv(paths['scores'],sep='\t');c=pd.read_csv(paths['clinical'],sep='\t').set_index('sampleID');qc=pd.read_csv(paths['qc'],sep='\t').set_index('unit_index');ref=pd.read_csv(paths['references'],sep='\t').set_index('reference')
    assert len(t)==len(old)==10600 and len(slots)==212 and not t.description_id.duplicated().any() and not slots.original_clinical_record.duplicated().any()
    assert t.description_id.tolist()==old.description_id.tolist();old=old.set_index('description_id');slots=slots.set_index('original_clinical_record')
    zfile=np.load(paths['readouts'],allow_pickle=False);arrays={k:zfile[k] for k in ['ULM_score','target_mean','nontarget_mean']};checked=0;maximum=0.;positive=0;negative=0
    for rec,g in t.groupby('original_clinical_record',sort=True):
        first=g.iloc[0];f=scores[scores.context.eq(first.context)&scores.type.eq(first.type)&scores.program.eq(first.program)&scores.sub_cell_type.eq(first.sub_cell_type)]
        cc=c.reindex(f.sampleID);use=f.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool);f=f.loc[use];cc=cc.loc[use];n=len(f);e=qc.reindex(f.unit_index);assert np.array_equal(e.n_cells.to_numpy(int),f.n_cells.to_numpy(int))
        x=f.raw_weighted_score.to_numpy(float);response=cc.response.to_numpy(float);ids=f.index.to_numpy(int);cov=np.column_stack([cc[COMMON].to_numpy(float),f[QUALITY].to_numpy(float)]);extended=np.column_stack([cov,e[QC].to_numpy(float)]);counts=np.bincount(response.astype(int),minlength=2)
        status='DESCRIBABLE_NO_NEW_P';z=None;z0=None;zr=None;rank=np.nan;cond=np.nan
        if n<50 or counts.min()<10:status='INSUFFICIENT_ORIGINAL_CASE_OR_RESPONSE_GROUP_COVERAGE'
        elif np.ptp(x)<=1e-12:status='CONSTANT_ORIGINAL_PROGRAMME'
        else:
            z0,base=make(cov,x)
            if base!='DESCRIBABLE_NO_NEW_P':status='ORIGINAL_'+base
            elif not np.isfinite(e[QC].to_numpy(float)).all():status='INCOMPLETE_QC_ORIGINAL_CASES_NOT_REMOVED'
            else:
                z,status=make(extended,x)
                if z is not None:
                    a=rankdata(x);a=(a-a.mean())/a.std();d=np.column_stack([z,a]);rank=int(np.linalg.matrix_rank(d));cond=float(np.linalg.cond(d));candidate=np.column_stack([z,response])
                    if full(np.column_stack([candidate,a]),17):zr=candidate
        s=slots.loc[rec];assert s.QC_design_status==status and s.eligible_cases==n and s.response0==counts[0] and s.response1==counts[1]
        for field in ['context','type','program','sub_cell_type','signature']:assert s[field]==first[field]
        close([s.design_rank,s.condition],[rank,cond])
        expected=[rho(x,e[QC[0]].to_numpy(float),z0) if z0 is not None else np.nan,rho(x,e[QC[1]].to_numpy(float),z0) if z0 is not None else np.nan,beta(x,response,z0),beta(x,response,z)]
        close([s.programme_mito_rho,s.programme_ribo_rho,s.programme_MPR_SD_beta_original,s.programme_MPR_SD_beta_QC],expected)
        for row in g.itertuples():
            oldrow=old.loc[row.description_id];assert row.eligible_cases==n==oldrow.eligible_cases and row.original_status==oldrow.status
            for field in ['original_clinical_record','context','type','program','sub_cell_type','signature','reference']:assert getattr(row,field)==oldrow[field]
            close([row.original_partial_rank_rho,row.original_functional_MPR_SD_beta],[oldrow.partial_rank_rho,oldrow.functional_MPR_minus_nonMPR_SD_beta])
            expected_status=status;main=np.nan;response_rho=np.nan;fg_rho=np.nan;bg_rho=np.nan;clinical=np.nan;change=np.nan;j=int(ref.loc[row.reference,'reference_column']);y=arrays['ULM_score'][ids,j]
            if oldrow.status!='DESCRIBABLE_NO_NEW_P':expected_status='ORIGINAL_'+oldrow.status
            elif z is not None:
                main=rho(x,y,z)
                if not np.isfinite(main):expected_status='NO_READOUT_VARIATION_AFTER_COMPLETE_QC_CONTROLS'
                else:
                    response_rho=rho(x,y,zr) if zr is not None else np.nan;fg_rho=rho(x,arrays['target_mean'][ids,j],z);bg_rho=rho(x,arrays['nontarget_mean'][ids,j],z);clinical=beta(y,response,z);change=main-oldrow.partial_rank_rho
            assert row.QC_status==expected_status
            close([row.QC_partial_rank_rho,row.QC_response_controlled_rho,row.QC_target_mean_rho,row.QC_nontarget_mean_rho,row.QC_functional_MPR_SD_beta,row.QC_minus_original_rho],[main,response_rho,fg_rho,bg_rho,clinical,change])
            if np.isfinite(main):maximum=max(maximum,abs(row.QC_partial_rank_rho-main));positive+=int(main>0);negative+=int(main<0)
            checked+=1
        (out/'private_independent_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),rows_checked=checked))+'\n')
    summary=pd.read_csv(out/'private_all1200_functional_QC_summary.tsv',sep='\t');assert len(summary)==1200 and not summary.duplicated(['signature','reference']).any()
    for row in summary.itertuples():
        f=t[t.signature.eq(row.signature)&t.reference.eq(row.reference)];v=f[f.QC_partial_rank_rho.notna()];assert row.planned_author_subtypes==len(f) and row.QC_describable_author_subtypes==len(v)
        close([row.old_median_same_subtypes,row.QC_median_same_subtypes,row.median_change_same_subtypes],[v.original_partial_rank_rho.median(),v.QC_partial_rank_rho.median(),v.QC_minus_original_rho.median()])
    assert checked==10600 and positive>0 and negative>0 and t.QC_partial_rank_rho.isna().any()
    proof=dict(status='PASS_ALL212_PROGRAMME_ALL10600_FUNCTIONAL_COMPLETE_CASE_CONTROL_STATUS_AND_OLS_POINTS_ALL1200_SUMMARIES_NO_NEW_P',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),description_rows_checked=checked,maximum_main_rho_difference=maximum,positive_negative_and_unavailable_cases_checked=True,new_P_BH=False)
    p=out/'private_independent_sensitivity_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
