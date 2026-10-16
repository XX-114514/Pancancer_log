"""Independent full OLS audit of third-lineage conditioning and capture denominators."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from audit_NSCLC_RNA_composition_sensitivity_v1 import COMMON,QUALITY,QC,rho,full,close

FIELDS=['function_AB_rho20_same_subset','function_AB_rho26_same_subset','function_AB_rho26_plus_third','function_AB_rho26_plus_third_response','function_A_third_rho26','function_B_third_rho26','programme_AB_rho26','programme_AB_rho26_plus_third','programme_A_third_rho26','programme_B_third_rho26','target_mean_AB_rho26','target_mean_AB_rho26_plus_third','function_change_after_third']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def standardized(x):return (x-x.mean(0))/x.std(0)


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=r/'NSCLC_third_lineage_background_descriptions_20261016_v1';job=r/'HNSCC_receiver_followup_job_control_20261016_v2/describe_NSCLC_third_lineage_background_v1';broot=r/'NSCLC_third_lineage_B_functions_20261016_v1'
    receipt=json.loads((job/'launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    terminal=json.loads((out/'private_terminal.json').read_text());plan=json.loads((out/'plan_before_any_third_lineage_association.json').read_text());assert plan['code_sha256']==receipt['source_sha256']
    for n,h in terminal['files'].items():assert sha(out/n)==h
    paths=dict(scores=r/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz',clinical=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',QC=r/'NSCLC_author_RNA_composition_20261016_v1/private_all_native_unit_RNA_composition.tsv',readouts=r/'NSCLC_disjoint_hallmark_measurement_20261016_v1/private_all_native_functional_scores_and_QA.npz',scope=r/'NSCLC_two_end_coordination_preflight_20261016_v1/private_all572_planned_pair_designs.tsv',references=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1/private_reference_axis.tsv',B_axis=broot/'private_all_B_case_axis_and_quality.tsv',B_masks=broot/'private_all16_mask_axis.tsv',B_readouts=broot/'private_all_B_function_scores_and_QA.npz')
    for k,p in paths.items():assert sha(p)==plan['input_sha256'][k]
    proof=broot/'private_independent_B_function_audit.json';assert sha(proof)==plan['B_measurement_audit_sha256'] and json.loads(proof.read_text())['status'].startswith('PASS_ALL_AUTHOR_B')
    t=pd.read_csv(out/'private_all28600_third_lineage_descriptions.tsv',sep='\t');capt=pd.read_csv(out/'private_all572_explicit_B_capture_opportunities.tsv',sep='\t').set_index('record');scope=pd.read_csv(paths['scope'],sep='\t');s=pd.read_csv(paths['scores'],sep='\t').assign(score_row=lambda f:np.arange(len(f)));c=pd.read_csv(paths['clinical'],sep='\t').set_index('sampleID');qc=pd.read_csv(paths['QC'],sep='\t').set_index('unit_index');B=pd.read_csv(paths['B_axis'],sep='\t').set_index('sampleID');masks=pd.read_csv(paths['B_masks'],sep='\t').set_index('mask_id');refs=pd.read_csv(paths['references'],sep='\t').sort_values('reference_column');ref_index=refs.set_index('reference')
    assert len(t)==28600 and len(capt)==572 and not t.duplicated(['record','reference']).any() and t.record.drop_duplicates().tolist()==scope.record.tolist();groups={k:f for k,f in t.groupby('record',sort=False)}
    z=np.load(paths['readouts'],allow_pickle=False);bz=np.load(paths['B_readouts'],allow_pickle=False);arrays={k:z[k] for k in ['ULM_score','target_mean','expressed_target_genes']};barrays={k:bz[k] for k in arrays};checked=0;maximum=0.;positive=0;negative=0;different_subset_cases_checked=False
    for pair in scope.itertuples():
        g=groups[pair.record];assert g.reference.tolist()==refs.reference.tolist()
        a=s[s.context.eq(pair.context)&s.program.eq(pair.program)&s.type.eq('typeA')&s.sub_cell_type.eq(pair.author_subtype_A)].set_index('sampleID');b=s[s.context.eq(pair.context)&s.program.eq(pair.program)&s.type.eq('typeB')&s.sub_cell_type.eq(pair.author_subtype_B)].set_index('sampleID')
        cases=sorted(set(a.index)&set(b.index));a=a.loc[cases];b=b.loc[cases];cc=c.loc[cases];use=a.n_cells.ge(20).to_numpy()&b.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool);a=a.loc[use];b=b.loc[use];cc=cc.loc[use];original_n=len(a);assert original_n==pair.eligible_paired_cases
        eligible_B=set(B.index[B.n_cells.ge(20)]);subset=[case for case in a.index if case in eligible_B];a=a.loc[subset];b=b.loc[subset];cc=cc.loc[subset];third=B.loc[subset];n=len(a);different_subset_cases_checked|=n!=original_n
        assert a.index.tolist()==b.index.tolist()==cc.index.tolist()==third.index.tolist()
        qa=qc.reindex(a.unit_index);qb=qc.reindex(b.unit_index);assert np.array_equal(qa.n_cells,a.n_cells) and np.array_equal(qb.n_cells,b.n_cells)
        cov20=np.column_stack([cc[COMMON].to_numpy(float),a[QUALITY].to_numpy(float),b[QUALITY].to_numpy(float),qa[QC].to_numpy(float),qb[QC].to_numpy(float)]);cov26=np.column_stack([cov20,third[QUALITY+QC].to_numpy(float)]);response=cc.response.to_numpy(float);counts=np.bincount(response.astype(int),minlength=2)
        x=a.raw_weighted_score.to_numpy(float);v=b.raw_weighted_score.to_numpy(float);ia=a.score_row.to_numpy(int);ib=b.score_row.to_numpy(int);it=third.pool_index.to_numpy(int);mi=int(masks.loc[f'C{pair.context}_{pair.program}','mask_index']);status='DESCRIBABLE_NO_NEW_P';Z=None;Z20=None;rank=np.nan;condition=np.nan
        if n<50 or counts.min()<10:status='INSUFFICIENT_EXPLICIT_B_CAPTURE_CASE_OR_RESPONSE_GROUP_COVERAGE'
        elif not np.isfinite(cov26).all() or np.any(cov26.std(0)<=1e-12):status='CONSTANT_OR_NONFINITE_COMPLETE26_CONTROL'
        else:
            Z=np.column_stack([np.ones(n),standardized(cov26)]);rank=int(np.linalg.matrix_rank(Z));condition=float(np.linalg.cond(Z))
            if not full(Z,27):status='NONIDENTIFIABLE_COMPLETE26_NUISANCE_SPACE';Z=None
            else:Z20=np.column_stack([np.ones(n),standardized(cov20)]);assert full(Z20,21)
        cv=capt.loc[pair.record];assert cv.context==pair.context and cv.program==pair.program and cv.base_status==status and cv.original_paired_cases==original_n and cv.B_captured_paired_cases==n and cv.excluded_for_B_capture==original_n-n and cv.response0==counts[0] and cv.response1==counts[1];close([cv.nuisance_rank,cv.condition],[rank,condition])
        for row in g.itertuples():
            for field in ['context','program','author_subtype_A','author_subtype_B','family_A','family_B']:assert getattr(row,field)==getattr(pair,field)
            assert row.original_paired_cases==original_n and row.B_captured_paired_cases==n and row.excluded_for_B_capture==original_n-n
            j=int(ref_index.loc[row.reference,'reference_column']);fa=arrays['ULM_score'][ia,j];fb=arrays['ULM_score'][ib,j];fc=barrays['ULM_score'][mi,it,j];fgA=arrays['target_mean'][ia,j];fgB=arrays['target_mean'][ib,j];fgC=barrays['target_mean'][mi,it,j]
            close([row.median_target_support_A,row.median_target_support_B,row.median_target_support_third],[np.median(arrays['expressed_target_genes'][ia,j]) if n else np.nan,np.median(arrays['expressed_target_genes'][ib,j]) if n else np.nan,np.median(barrays['expressed_target_genes'][mi,it,j]) if n else np.nan])
            expected={k:np.nan for k in FIELDS};expected_status=status;extra=False
            if Z is not None:
                if not all(np.isfinite(y).all() for y in [fa,fb,fc]):expected_status='NONFINITE_FUNCTION_ORIGINAL_B_CAPTURE_CASES_NOT_REMOVED'
                elif not all(np.any(y>0) for y in [fgA,fgB,fgC]):expected_status='BACKGROUND_ONLY_ALL_ZERO_FOREGROUND_IN_ONE_OR_MORE_LINEAGES'
                elif any(np.ptp(y)<=1e-12 for y in [fa,fb,fc]):expected_status='CONSTANT_FUNCTION_IN_ONE_OR_MORE_LINEAGES'
                else:
                    rc=standardized(rankdata(fc));ra=standardized(rankdata(fa));Zc=np.column_stack([Z,rc]);D=np.column_stack([Zc,ra])
                    if not full(D,29):expected_status='NONIDENTIFIABLE_COMPLETE26_PLUS_THIRD_FUNCTION_DESIGN'
                    else:
                        main=rho(fa,fb,Zc)
                        if not np.isfinite(main):expected_status='NO_TARGET_FUNCTION_VARIATION_AFTER_COMPLETE_CONTROLS'
                        else:
                            expected['function_AB_rho20_same_subset']=rho(fa,fb,Z20);expected['function_AB_rho26_same_subset']=rho(fa,fb,Z);expected['function_AB_rho26_plus_third']=main;expected['function_A_third_rho26']=rho(fa,fc,Z);expected['function_B_third_rho26']=rho(fb,fc,Z);expected['programme_AB_rho26']=rho(x,v,Z);expected['programme_AB_rho26_plus_third']=rho(x,v,Zc);expected['programme_A_third_rho26']=rho(x,fc,Z);expected['programme_B_third_rho26']=rho(v,fc,Z);expected['target_mean_AB_rho26']=rho(fgA,fgB,Z);expected['target_mean_AB_rho26_plus_third']=rho(fgA,fgB,Zc);expected['function_change_after_third']=main-expected['function_AB_rho26_same_subset']
                            Zcr=np.column_stack([Zc,response])
                            if full(np.column_stack([Zcr,ra]),30):extra=True;expected['function_AB_rho26_plus_third_response']=rho(fa,fb,Zcr)
            assert row.status==expected_status and bool(row.extra_response_estimable)==extra;close([getattr(row,k) for k in FIELDS],[expected[k] for k in FIELDS])
            main=expected['function_AB_rho26_plus_third']
            if np.isfinite(main):maximum=max(maximum,abs(row.function_AB_rho26_plus_third-main));positive+=int(main>0);negative+=int(main<0)
            checked+=1
        (out/'private_independent_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),description_rows_checked=checked))+'\n')
    assert checked==28600 and positive>0 and negative>0 and t.function_AB_rho26_plus_third.isna().any() and different_subset_cases_checked
    result=dict(status='PASS_ALL28600_THIRD_LINEAGE_CASE_CONTROL_STATUS_AND_ALL_OLS_POINTS_ALL572_CAPTURE_DENOMINATORS',finished_UTC=datetime.now(timezone.utc).isoformat(),audit_source_sha256=sha(__file__),OLS_helper_sha256=sha(Path(__file__).with_name('audit_NSCLC_RNA_composition_sensitivity_v1.py')),producer_terminal_sha256=sha(out/'private_terminal.json'),rows_checked=checked,maximum_main_rho_difference=maximum,positive_negative_unavailable_and_changed_case_subset_checked=True,new_P_BH=False)
    p=out/'private_independent_third_lineage_audit.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
