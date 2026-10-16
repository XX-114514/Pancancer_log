"""Independent all-map summaries and4800 same-patient OLS residual checks."""
import hashlib,itertools,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def correlation(x,y):
    return pd.Series(x).rank().corr(pd.Series(y).rank()) if len(x)>=10 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan


def conditional(x,y,controls):
    if len(x)<20:return np.nan,'INSUFFICIENT_SAME_PATIENTS'
    if len(x)-controls.shape[1]-2<10:return np.nan,'INSUFFICIENT_DF_WITH_COMPLETE_CONTROLS'
    if not np.isfinite(controls).all() or any(np.ptp(c)==0 for c in controls.T):return np.nan,'MISSING_OR_CONSTANT_REQUIRED_CONTROL'
    ranked=pd.DataFrame(controls).rank().to_numpy();ranked=(ranked-ranked.mean(axis=0))/ranked.std(axis=0);design=sm.add_constant(ranked,has_constant='add');rx=pd.Series(x).rank().to_numpy();ry=pd.Series(y).rank().to_numpy();full=np.column_stack([design,rx])
    if np.ptp(x)==0 or np.ptp(y)==0 or np.linalg.matrix_rank(full)!=full.shape[1] or np.linalg.cond(full)>1e8:return np.nan,'CONSTANT_OR_COLLINEAR_COMPLETE_DESIGN'
    ex=sm.OLS(rx,design).fit().resid;ey=sm.OLS(ry,design).fit().resid
    if np.linalg.norm(ex)<1e-10*np.linalg.norm(rx-rx.mean()) or np.linalg.norm(ey)<1e-10*np.linalg.norm(ry-ry.mean()):return np.nan,'FUNCTION_FULLY_EXPLAINED_BY_REQUIRED_CONTROLS'
    return float(np.dot(ex-ex.mean(),ey-ey.mean())/np.sqrt(np.dot(ex-ex.mean(),ex-ex.mean())*np.dot(ey-ey.mean(),ey-ey.mean()))),'DESCRIPTIVE_COMPLETE_CONTROL_RANK_RESIDUAL_NO_P'


def check(got,wanted):assert np.isclose(got,wanted,rtol=1e-9,atol=1e-10,equal_nan=True),(got,wanted)


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'HNSCC_complete_multicell_function_review_20261016_v3';rc=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/review_HNSCC_complete_multicell_function_changes_v3/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(rc['pid'])/'stat';assert p.exists();a=p.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(5)
    t=json.loads((out/'private_terminal.json').read_text());assert t['code_sha256']==rc['source_sha256']
    for name,h in t['files'].items():assert sha(out/name)==h
    source=r/'HNSCC_paired_multicell_function_changes_20261016_v2';assert sha(source/'private_independent_paired_change_audit.json')==t['input_audit_sha256'];all_results=pd.read_csv(source/'private_all3200_paired_functional_change_results.tsv',sep='\t');summary=pd.read_csv(out/'private_all200_function_pool_summaries.tsv',sep='\t');pairs=pd.read_csv(out/'private_all4800_same_patient_change_coordination.tsv',sep='\t');cases=pd.read_csv(source/'private_all_slot_exact_patient_pair_differences.tsv.gz',sep='\t',dtype={'patient_block':str});native=r/'HNSCC_native_multicell_functions_20261016_v1';axis=pd.read_csv(native/'private_all228_pool_library_axis_and_quality.tsv',sep='\t').set_index('pool_index');archive=np.load(native/'private_all_multicell_function_scores_and_QA.npz',allow_pickle=False);scores=archive['ULM_score'];masks=pd.read_csv(native/'private_all16_mask_axis.tsv',sep='\t');refs=pd.read_csv(native/'private_all50_reference_axis.tsv',sep='\t');mi={x:i for i,x in enumerate(masks.mask_id)};ri={x:i for i,x in enumerate(refs.reference)}
    assert len(summary)==200 and not summary.duplicated(['reference','pool']).any();assert len(pairs)==4800 and not pairs.duplicated(['mask_id','reference','pool_A','pool_B']).any()
    assert set(map(tuple,pairs[['mask_id','reference','pool_A','pool_B']].to_numpy()))=={(m,f,a,b) for m in masks.mask_id for f in refs.reference for a,b in itertools.combinations(sorted(axis.pool.unique()),2)}
    for row in summary.itertuples():
        g=all_results[all_results.reference.eq(row.reference)&all_results.pool.eq(row.pool)];assert len(g)==row.mask_count==16 and row.estimable_masks==g.status.str.startswith('ESTIMATED_').sum() and row.min_paired_patients==g.patients.min() and row.max_paired_patients==g.patients.max()
        for prefix,field in [('primary_BH','primary_sign_fdr_bh'),('primary_BY','primary_sign_fdr_by'),('companion_BH','companion_Wilcoxon_fdr_bh'),('companion_BY','companion_Wilcoxon_fdr_by')]:assert getattr(row,prefix+'_masks')==int((g[field]<.05).sum())
        check(row.median_patient_signed_fraction,np.median(2*g.positive_fraction-g.nonzero_patients/g.patients));check(row.median_raw_ULM_change,g.median_change.median());check(row.min_primary_BH,g.primary_sign_fdr_bh.min());check(row.max_primary_BH,g.primary_sign_fdr_bh.max());check(row.median_foreground_RNA_change,g.median_change_target_mean.median());check(row.median_nontarget_RNA_change,g.median_change_nontarget_mean.median());assert row.all_mask_change_directions_equal==(g.median_change.notna().all() and len(set(np.sign(g.median_change)))==1) and not row.independent_mask_replications
    groups={k:g.set_index('patient_block') for k,g in cases.groupby(['mask_id','reference','pool'],sort=False)};fields=['n_cells','fraction_of_all_source_library_cells','mean_logRNA','mean_loggenes','mean_symbol_pct_mito','mean_symbol_pct_ribo'];third='extension::B_plus_Plasma';estimated=0;unavailable=0
    def quality(g):
        pre=axis.loc[g.pre_pool_index,fields].to_numpy(float);post=axis.loc[g.post_pool_index,fields].to_numpy(float);pre[:,0]=np.log1p(pre[:,0]);post[:,0]=np.log1p(post[:,0]);return post-pre
    for row in pairs.itertuples():
        a=groups[(row.mask_id,row.reference,row.pool_A)];b=groups[(row.mask_id,row.reference,row.pool_B)];people=sorted(set(a.index)&set(b.index));aa=a.loc[people];bb=b.loc[people];x=aa.rounded_difference.to_numpy();y=bb.rounded_difference.to_numpy();check(row.raw_change_rho,correlation(x,y));assert row.same_patients==len(people) and row.quality12_required_controls==12 and not row.new_P_or_clinical_claim
        controls=np.column_stack([quality(aa),quality(bb)]);value,status=conditional(x,y,controls);check(row.quality12_adjusted_change_rho,value);assert row.quality12_status==status;estimated+=np.isfinite(value);unavailable+=not np.isfinite(value)
        i=mi[row.mask_id];j=ri[row.reference];check(row.baseline_same_patient_rho,correlation(scores[i,aa.pre_pool_index.to_numpy(int),j],scores[i,bb.pre_pool_index.to_numpy(int),j]));check(row.post_same_patient_rho,correlation(scores[i,aa.post_pool_index.to_numpy(int),j],scores[i,bb.post_pool_index.to_numpy(int),j]))
        if third not in [row.pool_A,row.pool_B]:
            bg=groups[(row.mask_id,row.reference,third)];people=sorted(set(people)&set(bg.index));aa=a.loc[people];bb=b.loc[people];x=aa.rounded_difference.to_numpy();y=bb.rounded_difference.to_numpy();z=bg.loc[people].rounded_difference.to_numpy();check(row.raw_rho_on_B_subset,correlation(x,y));assert row.third_B_same_patients==len(people)
            value,status=conditional(x,y,z[:,None]);check(row.B_function_only_adjusted_rho,value);assert row.B_only_status==status;value,status=conditional(x,y,np.column_stack([quality(aa),quality(bb),z]));check(row.quality12_plus_B_rho,value);assert row.quality12_plus_B_status==status
        else:assert row.third_B_same_patients==0 and row.B_only_status==row.quality12_plus_B_status=='B_ALREADY_ENDPOINT' and np.isnan(row.B_function_only_adjusted_rho)
    proof=dict(status='PASS_ALL200_MAP_SUMMARIES_ALL4800_EXACT_PATIENT_CHANGE_RAW_AND_COMPLETE_QC_OLS_CORRELATIONS',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),summary_rows_checked=len(summary),same_patient_pair_rows_checked=len(pairs),complete_QC_estimable=int(estimated),complete_QC_unavailable=int(unavailable),independent_statsmodels_OLS_not_residual_reranking=True,new_P_BH=False)
    p=out/'private_independent_complete_review_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
