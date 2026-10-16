"""All original pairs/functions, matched-case conditioning on third-lineage RNA."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from describe_NSCLC_disjoint_hallmark_v1 import COMMON,QUALITY,correlation,identifiable

QC=['mean_author_pct_mito','mean_author_pct_ribo']
METRICS=['function_AB_rho20_same_subset','function_AB_rho26_same_subset','function_AB_rho26_plus_third','function_AB_rho26_plus_third_response','function_A_third_rho26','function_B_third_rho26','programme_AB_rho26','programme_AB_rho26_plus_third','programme_A_third_rho26','programme_B_third_rho26','target_mean_AB_rho26','target_mean_AB_rho26_plus_third','function_change_after_third']


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def standardize(v):return (v-v.mean(0))/v.std(0)


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);broot=r/'NSCLC_third_lineage_B_functions_20261016_v1';proof=broot/'private_independent_B_function_audit.json'
    receipt=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_third_lineage_B_functions_v1/launch_receipt.json').read_text());deadline=time.monotonic()+14400
    while not proof.exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    assert json.loads(proof.read_text())['status'].startswith('PASS_ALL_AUTHOR_B') and sha(broot/'private_terminal.json')==json.loads(proof.read_text())['producer_terminal_sha256']
    paths=dict(scores=r/'NSCLC_context_native_programme_projection_20261016_v1/private_native_author_subtype_programme_scores.tsv.gz',clinical=r/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv',QC=r/'NSCLC_author_RNA_composition_20261016_v1/private_all_native_unit_RNA_composition.tsv',readouts=r/'NSCLC_disjoint_hallmark_measurement_20261016_v1/private_all_native_functional_scores_and_QA.npz',scope=r/'NSCLC_two_end_coordination_preflight_20261016_v1/private_all572_planned_pair_designs.tsv',references=r/'NSCLC_disjoint_hallmark_preflight_20261016_v1/private_reference_axis.tsv',B_axis=broot/'private_all_B_case_axis_and_quality.tsv',B_masks=broot/'private_all16_mask_axis.tsv',B_readouts=broot/'private_all_B_function_scores_and_QA.npz')
    hashes={k:sha(p) for k,p in paths.items()};previous=json.loads((r/'NSCLC_matched_function_coordination_20261016_v1/plan_before_any_matched_function_effect.json').read_text())
    for k in ['scores','clinical','QC','readouts','scope','references']:assert hashes[k]==previous['input_sha256'][k]
    bt=json.loads((broot/'private_terminal.json').read_text())
    for k in ['B_axis','B_masks','B_readouts']:assert hashes[k]==bt['files'][paths[k].name]
    out=r/'NSCLC_third_lineage_background_descriptions_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),input_sha256=hashes,B_measurement_audit_sha256=sha(proof),previous_question_scope_sha256=sha(r/'NSCLC_third_lineage_background_plan_20261016_v1/PRIVATE_NEXT_QUESTION_AND_SCOPE.json'),
              scope='Every original572 pair x50 functions,28600 slots. Original paired eligibility then explicit third-lineage B>=20 subset; no effect selection.',
              before_after='Exactly same B-captured cases: first original20 controls, then all26 including B fraction/logN/logRNA/loggenes/author mito/ribo, then add ranked B functional readout.',
              gates='n>=50 and10 per old pathology group; all26 controls complete/nonconstant,27 nuisance columns and full29 columns with B-functional rank and A-functional rank; df>=10,condition<1e8. Do not drop controls/cases to repair a readout.',
              descriptions='Function A-B before/after B functional rank; function and original programme A/B with third functional rank; target-RNA means; extra pathology sensitivity. No new P/BH/CI; not causal mediation or a noninteracting negative control.',
              zero='An all-zero target foreground in any of three lineages is background-only and unavailable for interpretation. Actual zero retained, missing not imputed.',
              unchanged='All old clinical/LR/coordination families, original labels and frozen weights remain unchanged.')
    (out/'plan_before_any_third_lineage_association.json').write_text(json.dumps(plan,indent=2)+'\n')
    s=pd.read_csv(paths['scores'],sep='\t').assign(score_row=lambda f:np.arange(len(f)));c=pd.read_csv(paths['clinical'],sep='\t').set_index('sampleID');qc=pd.read_csv(paths['QC'],sep='\t').set_index('unit_index');scope=pd.read_csv(paths['scope'],sep='\t');refs=pd.read_csv(paths['references'],sep='\t').sort_values('reference_column');B=pd.read_csv(paths['B_axis'],sep='\t').set_index('sampleID');masks=pd.read_csv(paths['B_masks'],sep='\t').set_index('mask_id')
    zz=np.load(paths['readouts'],allow_pickle=False);bz=np.load(paths['B_readouts'],allow_pickle=False);arrays={k:zz[k] for k in ['ULM_score','target_mean','expressed_target_genes']};barrays={k:bz[k] for k in arrays};assert B.index.is_unique and len(scope)==572 and len(refs)==50
    rows=[];coverage=[]
    for pair in scope.itertuples():
        a=s[s.context.eq(pair.context)&s.program.eq(pair.program)&s.type.eq('typeA')&s.sub_cell_type.eq(pair.author_subtype_A)].set_index('sampleID');b=s[s.context.eq(pair.context)&s.program.eq(pair.program)&s.type.eq('typeB')&s.sub_cell_type.eq(pair.author_subtype_B)].set_index('sampleID');assert a.index.is_unique and b.index.is_unique
        cases=sorted(set(a.index)&set(b.index));a=a.loc[cases];b=b.loc[cases];cc=c.loc[cases];use=a.n_cells.ge(20).to_numpy()&b.n_cells.ge(20).to_numpy()&cc.clinical_eligible.to_numpy(bool);a=a.loc[use];b=b.loc[use];cc=cc.loc[use];original_n=len(a);assert original_n==pair.eligible_paired_cases
        third=B.reindex(a.index);captured=third.n_cells.ge(20).to_numpy();a=a.loc[captured];b=b.loc[captured];cc=cc.loc[captured];third=third.loc[captured];n=len(a);assert a.index.tolist()==b.index.tolist()==cc.index.tolist()==third.index.tolist()
        x=a.raw_weighted_score.to_numpy(float);v=b.raw_weighted_score.to_numpy(float);response=cc.response.to_numpy(float);ia=a.score_row.to_numpy(int);ib=b.score_row.to_numpy(int);it=third.pool_index.to_numpy(int);mi=int(masks.loc[f'C{pair.context}_{pair.program}','mask_index']);counts=np.bincount(response.astype(int),minlength=2)
        qa=qc.reindex(a.unit_index);qb=qc.reindex(b.unit_index);assert np.array_equal(qa.n_cells.to_numpy(int),a.n_cells.to_numpy(int)) and np.array_equal(qb.n_cells.to_numpy(int),b.n_cells.to_numpy(int))
        cov20=np.column_stack([cc[COMMON].to_numpy(float),a[QUALITY].to_numpy(float),b[QUALITY].to_numpy(float),qa[QC].to_numpy(float),qb[QC].to_numpy(float)]);cov26=np.column_stack([cov20,third[QUALITY+QC].to_numpy(float)])
        base='DESCRIBABLE_NO_NEW_P';Z=None;q26=None;q20=None;rank=np.nan;condition=np.nan
        if n<50 or counts.min()<10:base='INSUFFICIENT_EXPLICIT_B_CAPTURE_CASE_OR_RESPONSE_GROUP_COVERAGE'
        elif not np.isfinite(cov26).all() or np.any(cov26.std(0)<=1e-12):base='CONSTANT_OR_NONFINITE_COMPLETE26_CONTROL'
        else:
            Z=np.column_stack([np.ones(n),standardize(cov26)]);rank=int(np.linalg.matrix_rank(Z));condition=float(np.linalg.cond(Z))
            if not identifiable(Z,27):base='NONIDENTIFIABLE_COMPLETE26_NUISANCE_SPACE'
            else:q26=np.linalg.qr(Z,mode='reduced')[0];Z20=np.column_stack([np.ones(n),standardize(cov20)]);assert identifiable(Z20,21);q20=np.linalg.qr(Z20,mode='reduced')[0]
        coverage.append(dict(record=pair.record,context=pair.context,program=pair.program,original_paired_cases=original_n,B_captured_paired_cases=n,excluded_for_B_capture=original_n-n,response0=int(counts[0]),response1=int(counts[1]),base_status=base,nuisance_rank=rank,condition=condition))
        for ref in refs.itertuples():
            j=int(ref.reference_column);fa=arrays['ULM_score'][ia,j];fb=arrays['ULM_score'][ib,j];fc=barrays['ULM_score'][mi,it,j];fgA=arrays['target_mean'][ia,j];fgB=arrays['target_mean'][ib,j];fgC=barrays['target_mean'][mi,it,j]
            row=dict(record=pair.record,context=pair.context,program=pair.program,author_subtype_A=pair.author_subtype_A,author_subtype_B=pair.author_subtype_B,family_A=pair.family_A,family_B=pair.family_B,reference=ref.reference,original_paired_cases=original_n,B_captured_paired_cases=n,excluded_for_B_capture=original_n-n,status=base,extra_response_estimable=False)
            row.update({key:np.nan for key in METRICS});row.update(median_target_support_A=float(np.median(arrays['expressed_target_genes'][ia,j])) if n else np.nan,median_target_support_B=float(np.median(arrays['expressed_target_genes'][ib,j])) if n else np.nan,median_target_support_third=float(np.median(barrays['expressed_target_genes'][mi,it,j])) if n else np.nan)
            if q26 is not None:
                if not all(np.isfinite(y).all() for y in [fa,fb,fc]):row['status']='NONFINITE_FUNCTION_ORIGINAL_B_CAPTURE_CASES_NOT_REMOVED'
                elif not all(np.any(y>0) for y in [fgA,fgB,fgC]):row['status']='BACKGROUND_ONLY_ALL_ZERO_FOREGROUND_IN_ONE_OR_MORE_LINEAGES'
                elif any(np.ptp(y)<=1e-12 for y in [fa,fb,fc]):row['status']='CONSTANT_FUNCTION_IN_ONE_OR_MORE_LINEAGES'
                else:
                    rc=standardize(rankdata(fc));ra=standardize(rankdata(fa));Zc=np.column_stack([Z,rc]);D=np.column_stack([Zc,ra])
                    if not identifiable(D,29):row['status']='NONIDENTIFIABLE_COMPLETE26_PLUS_THIRD_FUNCTION_DESIGN'
                    else:
                        q3=np.linalg.qr(Zc,mode='reduced')[0];main=correlation(fa,fb,q3)
                        if not np.isfinite(main):row['status']='NO_TARGET_FUNCTION_VARIATION_AFTER_COMPLETE_CONTROLS'
                        else:
                            row['function_AB_rho20_same_subset']=correlation(fa,fb,q20);row['function_AB_rho26_same_subset']=correlation(fa,fb,q26);row['function_AB_rho26_plus_third']=main
                            row['function_A_third_rho26']=correlation(fa,fc,q26);row['function_B_third_rho26']=correlation(fb,fc,q26);row['programme_AB_rho26']=correlation(x,v,q26);row['programme_AB_rho26_plus_third']=correlation(x,v,q3);row['programme_A_third_rho26']=correlation(x,fc,q26);row['programme_B_third_rho26']=correlation(v,fc,q26)
                            row['target_mean_AB_rho26']=correlation(fgA,fgB,q26);row['target_mean_AB_rho26_plus_third']=correlation(fgA,fgB,q3);row['function_change_after_third']=main-row['function_AB_rho26_same_subset']
                            Zcr=np.column_stack([Zc,response])
                            if identifiable(np.column_stack([Zcr,ra]),30):row['extra_response_estimable']=True;row['function_AB_rho26_plus_third_response']=correlation(fa,fb,np.linalg.qr(Zcr,mode='reduced')[0])
            rows.append(row)
        (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),description_rows=len(rows),pairs=len(coverage)))+'\n')
    table=pd.DataFrame(rows);capt=pd.DataFrame(coverage);assert len(table)==28600 and len(capt)==572 and not table.duplicated(['record','reference']).any();table.to_csv(out/'private_all28600_third_lineage_descriptions.tsv',sep='\t',index=False,na_rep='NA');capt.to_csv(out/'private_all572_explicit_B_capture_opportunities.tsv',sep='\t',index=False,na_rep='NA')
    terminal=dict(status='COMPLETE_ALL28600_THIRD_LINEAGE_DESCRIPTIONS_AND_ALL572_EXPLICIT_CAPTURE_OPPORTUNITIES_AUDIT_PENDING_NO_NEW_P',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),description_rows=len(table),pair_slots=len(capt),status_counts=table.status.value_counts().to_dict(),files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_P_BH=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
