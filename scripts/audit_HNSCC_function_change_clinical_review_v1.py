"""Recheck every joined clinical row, source case set, RNA QA and map summary."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,statsmodels.api as sm


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def same(a,b):assert np.isclose(a,b,rtol=1e-9,atol=1e-10,equal_nan=True),(a,b)


def truth(v):return str(v).lower() in {'true','1'}


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'HNSCC_function_change_clinical_complete_review_20261016_v1';rc=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/review_HNSCC_function_change_clinical_v1/launch_receipt.json').read_text());deadline=time.monotonic()+43200
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(rc['pid'])/'stat';assert p.exists();a=p.read_text().rsplit(')',1)[1].split();assert a[0]!='Z' and a[19]==str(rc['proc_start_ticks']) and time.monotonic()<deadline;time.sleep(5)
    t=json.loads((out/'private_terminal.json').read_text());assert t['code_sha256']==rc['source_sha256']
    for name,h in t['files'].items():assert sha(out/name)==h
    source=r/'HNSCC_function_change_clinical_results_20261016_v1';pre=r/'HNSCC_function_change_clinical_preflight_20261016_v1';original=pd.read_csv(source/'private_all8000_longitudinal_function_clinical_results.tsv',sep='\t',keep_default_na=False,na_values=['NA']).set_index('slot');evidence=pd.read_csv(out/'private_all8000_exact_case_multilevel_clinical_evidence.tsv',sep='\t',keep_default_na=False,na_values=['NA']).set_index('slot');assert len(evidence)==8000 and evidence.index.is_unique and set(evidence.index)==set(original.index)
    for col in original:pd.testing.assert_series_equal(evidence.loc[original.index,col],original[col],check_names=False,check_dtype=False,rtol=1e-12,atol=1e-12)
    meta=pd.read_csv(pre/'private_unique_numeric_designs.tsv',sep='\t');models={}
    for row in meta.itertuples():
        a=np.load(pre/'private_designs'/(row.model_id+'.npz'),allow_pickle=False);x=a['D'][:,-1];z=a['D'][:,:-1];y=a['Y'];ex=sm.OLS(x,z).fit().resid;ey=sm.OLS(y,z).fit().resid;models[row.model_id]=(a['patients'].astype(str).tolist(),float(pd.Series(x).corr(pd.Series(y))),float(pd.Series(ex).corr(pd.Series(ey))))
    native=r/'HNSCC_native_multicell_functions_20261016_v1';masks=pd.read_csv(native/'private_all16_mask_axis.tsv',sep='\t');refs=pd.read_csv(native/'private_all50_reference_axis.tsv',sep='\t');mi={v:i for i,v in enumerate(masks.mask_id)};ri={v:i for i,v in enumerate(refs.reference)};a=np.load(native/'private_all_multicell_function_scores_and_QA.npz',allow_pickle=False);values={k:a[k] for k in a.files};cases=pd.read_csv(r/'HNSCC_paired_multicell_function_changes_20261016_v2/private_all_slot_exact_patient_pair_differences.tsv.gz',sep='\t',dtype={'patient_block':str});groups={k:g.set_index('patient_block') for k,g in cases.groupby(['mask_id','reference','pool'],sort=False)};turn=pd.read_csv(pre/'private_all_patient_pool_fine_composition_turnover.tsv',sep='\t',dtype={'patient_block':str}).set_index(['patient_block','pool']);coord=pd.read_csv(r/'HNSCC_complete_multicell_function_review_20261016_v3/private_all4800_same_patient_change_coordination.tsv',sep='\t').set_index(['mask_id','reference','pool_A','pool_B'])
    for slot,row in evidence.iterrows():
        endpoints=['A'] if row.family=='single3200' else ['A','B'];patients=None;all_ci=True;all_loo=True
        for label in endpoints:
            if row[label+'_model'] not in models:all_ci=False;all_loo=False;continue
            people,raw,partial=models[row[label+'_model']]
            if patients is None:patients=people
            else:assert patients==people
            same(row[label+'_raw_clinical_rank_rho'],raw);same(row[label+'_complete_model_partial_rank_rho'],partial);pool=row['pool_'+label];g=groups[(row.mask_id,row.reference,pool)].loc[people];same(row[label+'_median_pool_composition_TV'],np.median([turn.loc[(p,pool),'composition_total_variation'] for p in people]));same(row[label+'_median_change_functional_RNA'],np.median(g.rounded_difference))
            i=mi[row.mask_id];j=ri[row.reference]
            for key in ['target_mean','nontarget_mean','expressed_target_genes']:same(row[label+'_same_patient_median_change_'+key],np.median(values[key][i,g.post_pool_index.to_numpy(int),j]-values[key][i,g.pre_pool_index.to_numpy(int),j]))
            lo=row[label+'_CI_low'];hi=row[label+'_CI_high'];all_ci &= np.isfinite([lo,hi]).all() and lo*hi>0;all_loo &= truth(row[label+'_LOO_same_sign']) and row[label+'_LOO_valid']==len(people)
        assert truth(row.all_endpoint_CI_exclude_zero)==bool(all_ci) and truth(row.all_endpoint_LOO_valid_same_sign)==bool(all_loo) and not truth(row.new_direct_communication_or_independent_validation_accepted);assert truth(row.same_exact_patient_sets_confirmed)==(patients is not None and len(patients)==row.clinical_complete_patients)
        if row.family=='pair4800':
            expected=sorted(set(groups[(row.mask_id,row.reference,row.pool_A)].index)&set(groups[(row.mask_id,row.reference,row.pool_B)].index));assert truth(row.historical_coordination_same_exact_patients)==(patients==expected);assert truth(row.clinical_endpoint_betas_same_sign)==bool(np.sign(row.A_beta)==np.sign(row.B_beta))
            if patients==expected:
                prior=coord.loc[(row.mask_id,row.reference,row.pool_A,row.pool_B)]
                for key in ['raw_change_rho','quality12_adjusted_change_rho','quality12_status','B_function_only_adjusted_rho','B_only_status','quality12_plus_B_rho','quality12_plus_B_status']:
                    if 'status' in key:assert row['matched_coordination_'+key]==prior[key]
                    else:same(row['matched_coordination_'+key],prior[key])
        if not str(row.association_status).startswith('ESTIMATED_'):category='INSUFFICIENT_CLINICAL_INFORMATION'
        elif row.BH_q>=.05:category='NO_FULL_FAMILY_SHARED_CLINICAL_SUPPORT' if len(endpoints)==2 else 'NO_FULL_FAMILY_SINGLE_POOL_CLINICAL_SUPPORT'
        elif row.BY_q<.05 and row.MC_upper_BH_q<.05 and all_ci and all_loo:category='SOURCE_LIMITED_SHARED_CLINICAL_CANDIDATE' if len(endpoints)==2 else 'SOURCE_LIMITED_SINGLE_POOL_CLINICAL_CANDIDATE'
        else:category='EXPLORATORY_CLINICAL_ASSOCIATION_WITH_WEAK_OR_SENSITIVE_DIAGNOSTICS'
        assert row.evidence_category==category
    for family,file,keys in [('pair4800','private_all300_pair_function_summaries.tsv',['reference','pool_A','pool_B']),('single3200','private_all200_single_function_summaries.tsv',['reference','pool_A'])]:
        summary=pd.read_csv(out/file,sep='\t');assert len(summary)==(300 if family=='pair4800' else 200) and not summary.duplicated(keys).any()
        for row in summary.itertuples():
            g=evidence[evidence.family.eq(family)]
            for key in keys:g=g[g[key].eq(getattr(row,key))]
            assert len(g)==row.mask_count==16 and row.primary_BH_masks==(g.BH_q<.05).sum() and row.BY_masks==(g.BY_q<.05).sum() and row.MCupper_BH_masks==(g.MC_upper_BH_q<.05).sum() and row.robust_source_masks==g.evidence_category.str.contains('SOURCE_LIMITED').sum();same(row.min_BH_q,g.BH_q.min());same(row.max_BH_q,g.BH_q.max());assert row.min_patients==g.clinical_complete_patients.min() and row.max_patients==g.clinical_complete_patients.max() and not row.independent_mask_replications
            for label in ['A'] if family=='single3200' else ['A','B']:
                for field in ['beta','raw_clinical_rank_rho','complete_model_partial_rank_rho','median_pool_composition_TV']:same(getattr(row,label+'_median_'+field),g[label+'_'+field].median())
    proof=dict(status='PASS_ALL8000_EXACT_CASE_CLINICAL_COORDINATION_RNA_QA_COMPOSITION_AND300_200_MAP_SUMMARIES',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(out/'private_terminal.json'),clinical_rows_checked=len(evidence),pair_function_summaries_checked=300,single_function_summaries_checked=200,original_scope_denominators_and_case_sets_checked=True,new_P_or_direct_communication_or_independent_validation=False)
    p=out/'private_independent_clinical_review_audit.json';assert not p.exists();p.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))


if __name__=='__main__':main()
