"""Independent complete patient-key, exact-test, interval and correction replay."""
import hashlib,json,os,time
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import binomtest


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def signed_probability(delta,cache):
    d=[float(x) for x in delta if x!=0];magnitudes=sorted(set(abs(x) for x in d));weights={}
    for v in magnitudes:
        below=sum(abs(x)<v for x in d);equal=sum(abs(x)==v for x in d);weights[v]=2*below+equal+1
    ranks=[weights[abs(x)] for x in d];key=tuple(sorted(ranks))
    if key not in cache:
        histogram=Counter({0:1})
        for rank in key:
            updated=histogram.copy()
            for total,count in histogram.items():updated[total+rank]+=count
            histogram=updated
        assert sum(histogram.values())==2**len(d);cache[key]=histogram
    histogram=cache[key];observed=sum(rank for rank,value in zip(ranks,d) if value>0);left=sum(n for s,n in histogram.items() if s<=observed);right=sum(n for s,n in histogram.items() if s>=observed)
    return min(1.,2*min(left,right)/2**len(d)),observed/2


def adjusted(p):
    n=len(p);order=np.argsort(p,kind='stable');ranked=p[order]*n/np.arange(1,n+1);q=np.minimum(1.,np.minimum.accumulate(ranked[::-1])[::-1]);out=np.empty(n);out[order]=q;return out,np.minimum(1.,out*sum(1/k for k in range(1,n+1)))


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'HNSCC_paired_multicell_function_changes_20261016_v1';rc=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/test_HNSCC_paired_multicell_function_changes_v1/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(rc['pid'])/'stat';assert p.exists();a=p.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(5)
    t=json.loads((out/'private_terminal.json').read_text());plan=json.loads((out/'plan_before_any_paired_function_effect.json').read_text());assert t['code_sha256']==rc['source_sha256']==plan['code_sha256']
    for n,h in t['files'].items():assert sha(out/n)==h
    source=r/'HNSCC_native_multicell_functions_20261016_v1'
    for name,h in plan['input_sha256'].items():assert sha(source/name)==h
    assert sha(source/'private_independent_multicell_function_audit.json')==plan['measurement_audit_sha256']
    axis=pd.read_csv(source/'private_all228_pool_library_axis_and_quality.tsv',sep='\t',dtype={'library_id':str,'patient_block':str}).set_index('pool_index');masks=pd.read_csv(source/'private_all16_mask_axis.tsv',sep='\t');refs=pd.read_csv(source/'private_all50_reference_axis.tsv',sep='\t');npz=np.load(source/'private_all_multicell_function_scores_and_QA.npz',allow_pickle=False);result=pd.read_csv(out/'private_all3200_paired_functional_change_results.tsv',sep='\t');cases=pd.read_csv(out/'private_all_slot_exact_patient_pair_differences.tsv.gz',sep='\t',dtype={'patient_block':str});exclusions=pd.read_csv(out/'private_all_patient_pool_capture_exclusions.tsv',sep='\t',dtype={'patient_block':str})
    key=['mask_id','reference','pool'];assert len(result)==3200 and not result.duplicated(key).any();assert set(map(tuple,result[key].to_numpy()))=={(m,ref,pool) for m in masks.mask_id for ref in refs.reference for pool in axis.pool.unique()};assert not cases.duplicated(key+['patient_block']).any()
    mask_index={x:i for i,x in enumerate(masks.mask_id)};reference_index={x:i for i,x in enumerate(refs.reference)};candidates={};expected_exclusions=[]
    for (pool,patient),g in axis.groupby(['pool','patient_block'],sort=True):
        phases={phase:g[g.timepoint.eq(phase)] for phase in ['pre','post']};reason='READY'
        if any(len(v)!=1 for v in phases.values()):reason='NO_UNIQUE_REVIEWED_PRE_POST'
        else:
            pre=phases['pre'].iloc[0];post=phases['post'].iloc[0]
            if min(pre.n_cells,post.n_cells)<20:reason='CAPTURE_LT20'
            elif pre.tissue!=post.tissue or pre.treatment_arm!=post.treatment_arm:reason='INCOMPARABLE_RECORDED_TISSUE_OR_ARM'
        if reason=='READY':candidates.setdefault(pool,[]).append((patient,int(phases['pre'].index[0]),int(phases['post'].index[0])))
        expected_exclusions.append((pool,patient,reason,int((~g.timepoint.isin(['pre','post'])).sum())))
    assert set(expected_exclusions)==set(map(tuple,exclusions[['pool','patient_block','status','other_phase_libraries']].to_numpy()))
    case_index={k:g for k,g in cases.groupby(key,sort=False)};cache={};primary=[];companion=[];tested=0;bootstraps=0;checked_cases=0
    for row in result.itertuples():
        mi=mask_index[row.mask_id];ri=reference_index[row.reference];m=masks.iloc[mi];assert row.context==m.context and row.program==m.program
        eligible=[]
        for patient,pre,post in candidates.get(row.pool,[]):
            y=npz['ULM_score'][mi,[pre,post],ri];support=npz['expressed_target_genes'][mi,[pre,post],ri]
            if np.isfinite(y).all() and min(support)>0:eligible.append((patient,pre,post,float(y[1]-y[0]),float(np.round(y[1]-y[0],8))))
        g=case_index.get((row.mask_id,row.reference,row.pool),pd.DataFrame(columns=['patient_block','pre_pool_index','post_pool_index','raw_difference','rounded_difference']));assert set(map(tuple,g[['patient_block','pre_pool_index','post_pool_index']].to_numpy()))=={v[:3] for v in eligible}
        gg=g.set_index('patient_block')
        for patient,pre,post,raw,rounded in eligible:assert abs(gg.loc[patient,'raw_difference']-raw)<1e-12 and gg.loc[patient,'rounded_difference']==rounded
        checked_cases+=len(eligible);delta=np.array([x[4] for x in eligible]);nonzero=delta[delta!=0];assert row.patients==len(eligible) and row.paired_capture_patients==len(candidates.get(row.pool,[])) and row.RNA_or_nonfinite_exclusions==row.paired_capture_patients-len(eligible) and row.nonzero_patients==len(nonzero)
        if len(delta)>=20 and len(nonzero)>=10:
            sign=float(binomtest(int((nonzero>0).sum()),len(nonzero),.5,alternative='two-sided').pvalue);signed,rank=signed_probability(delta,cache);assert row.status=='ESTIMATED_LONGITUDINAL_FUNCTIONAL_SHIFT_NOT_COMMUNICATION'
            assert abs(row.primary_sign_P-sign)<1e-13 and abs(row.companion_Wilcoxon_P-signed)<1e-13 and row.Wilcoxon_positive_rank==rank
            assert abs(row.median_change-float(np.median(delta)))<1e-12 and abs(row.positive_fraction-float((delta>0).mean()))<1e-13 and abs(row.positive_nonzero_fraction-float((nonzero>0).mean()))<1e-13
            slot=f'{row.mask_id}|{row.reference}|{row.pool}';seed=int.from_bytes(hashlib.sha256(slot.encode()).digest()[:8],'little');rng=np.random.default_rng(seed);idx=rng.integers(0,len(delta),(2000,len(delta)));medians=np.median(delta[idx],axis=1);limits=np.percentile(medians,[2.5,97.5]);assert np.allclose([row.median_change_low,row.median_change_high],limits,rtol=1e-12,atol=1e-12);tested+=1;bootstraps+=2000
        else:assert row.status=='INSUFFICIENT_PAIRED_OR_NONZERO_PATIENTS' and row.primary_sign_P==row.companion_Wilcoxon_P==1 and np.isnan(row.median_change)
        if eligible:
            before=np.array([x[1] for x in eligible]);after=np.array([x[2] for x in eligible])
            for k in ['target_mean','nontarget_mean','expressed_target_genes']:assert np.isclose(getattr(row,'median_change_'+k),np.median(npz[k][mi,after,ri]-npz[k][mi,before,ri]),rtol=1e-12,atol=1e-12)
            for k in ['n_cells','fraction_of_all_source_library_cells','mean_logRNA','mean_loggenes','mean_symbol_pct_mito','mean_symbol_pct_ribo']:assert np.isclose(getattr(row,'median_change_'+k),np.median(axis.loc[after,k].to_numpy()-axis.loc[before,k].to_numpy()),rtol=1e-12,atol=1e-12)
        assert not row.new_LR_or_clinical_prediction_claim;primary.append(row.primary_sign_P);companion.append(row.companion_Wilcoxon_P)
    for prefix,ps in [('primary_sign',primary),('companion_Wilcoxon',companion)]:
        bh,by=adjusted(np.array(ps));assert np.allclose(result[prefix+'_fdr_bh'],bh,rtol=1e-12,atol=1e-14) and np.allclose(result[prefix+'_fdr_by'],by,rtol=1e-12,atol=1e-14)
    assert checked_cases==len(cases)
    proof=dict(status='PASS_ALL3200_EXACT_PATIENT_PAIRS_SIGN_WILCOXON_BOOTSTRAP_QC_AND_BOTH_FULL_FAMILY_BH_BY',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),slots_checked=len(result),estimated_slots=tested,patient_pair_rows_checked=checked_cases,bootstrap_draws_replayed=bootstraps,old_clinical_and_LR_families_unchanged=True,new_direct_communication_or_causal_claim=False)
    p=out/'private_independent_paired_change_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
