"""Complete private third-lineage summaries/maps after full independent OLS acceptance."""
import hashlib,json,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from report_NSCLC_functional_coordination_private_v2 import heatmaps


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);source=r/'NSCLC_third_lineage_background_descriptions_20261016_v1';proof=source/'private_independent_third_lineage_audit.json';deadline=time.monotonic()+14400
    receipt=json.loads((r/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_third_lineage_background_v1/launch_receipt.json').read_text())
    while not proof.exists():
        p=Path('/proc')/str(receipt['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split();assert s[19]==str(receipt['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(5)
    assert json.loads(proof.read_text())['status'].startswith('PASS_ALL28600_THIRD_LINEAGE') and sha(source/'private_terminal.json')==json.loads(proof.read_text())['producer_terminal_sha256']
    out=r/'NSCLC_third_lineage_complete_review_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'private_display_plan.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),source_audit_sha256=sha(proof),layout_helper_sha256=sha(Path(__file__).with_name('report_NSCLC_functional_coordination_private_v2.py')),scope='All572x50 descriptions, all8 original two-end parents and all50 functions, no effect-selected display or new inference. Same-case medians share subjects, not independent replications.'),indent=2)+'\n')
    t=pd.read_csv(source/'private_all28600_third_lineage_descriptions.tsv',sep='\t');capture=pd.read_csv(source/'private_all572_explicit_B_capture_opportunities.tsv',sep='\t');rows=[]
    fields=['function_AB_rho20_same_subset','function_AB_rho26_same_subset','function_AB_rho26_plus_third','function_AB_rho26_plus_third_response','function_A_third_rho26','function_B_third_rho26','programme_AB_rho26','programme_AB_rho26_plus_third','programme_A_third_rho26','programme_B_third_rho26','target_mean_AB_rho26','target_mean_AB_rho26_plus_third','function_change_after_third']
    for (context,program,reference),f in t.groupby(['context','program','reference'],sort=True):
        v=f[f.function_AB_rho26_plus_third.notna()];row=dict(context=context,program=program,parent=f'C{context}_{program}',reference=reference,planned_subtype_pairs=len(f),describable_subtype_pairs=len(v),minimum_B_subset_cases=int(v.B_captured_paired_cases.min()) if len(v) else -1,maximum_B_subset_cases=int(v.B_captured_paired_cases.max()) if len(v) else -1,same_cases_shared_subtypes_not_independent_replication=True)
        for field in fields:row['median_'+field]=v[field].median()
        rows.append(row)
    a=pd.DataFrame(rows);assert len(a)==400;a.to_csv(out/'private_all400_parent_function_background_summaries.tsv',sep='\t',index=False,na_rep='NA')
    c=[]
    for (context,program),f in capture.groupby(['context','program'],sort=True):
        c.append(dict(context=context,program=program,planned_subtype_pairs=len(f),B_subset_identifiable_pairs=int(f.base_status.eq('DESCRIBABLE_NO_NEW_P').sum()),pairs_with_changed_case_set=int(f.excluded_for_B_capture.gt(0).sum()),median_original_cases=f.original_paired_cases.median(),median_B_subset_cases=f.B_captured_paired_cases.median(),minimum_original_cases=int(f.original_paired_cases.min()),maximum_original_cases=int(f.original_paired_cases.max()),minimum_B_subset_cases=int(f.B_captured_paired_cases.min()),maximum_B_subset_cases=int(f.B_captured_paired_cases.max())))
    pd.DataFrame(c).to_csv(out/'private_all8_parent_capture_summaries.tsv',sep='\t',index=False)
    parents=sorted(a.parent.unique());refs=sorted(a.reference.unique());before=a.pivot(index='parent',columns='reference',values='median_function_AB_rho26_same_subset').reindex(index=parents,columns=refs);after=a.pivot(index='parent',columns='reference',values='median_function_AB_rho26_plus_third').reindex(index=parents,columns=refs)
    bounds=heatmaps(out/'private_all_parent_functions_before_after_B_background.png',[before.to_numpy(),after.to_numpy()],parents,refs,['A. Same B-captured cases: all 26 clinical and three-lineage composition/quality controls','B. Identical cases and controls + ranked B-lineage function'],
                    'All50 functions, fixed learning-gene exclusions. Shared-case subtype medians are not independent replications.\nB cells can participate biologically. Conditioning is not causal mediation or proof of artifact. No new P/q.')
    np.savez_compressed(out/'private_complete_display_arrays.npz',parents=np.asarray(parents,str),references=np.asarray(refs,str),before_B_function=before.to_numpy(),after_B_function=after.to_numpy())
    # Closed-key reread validates all summary numbers, apart from already independent science audits.
    disk=pd.read_csv(out/'private_all400_parent_function_background_summaries.tsv',sep='\t')
    for row in disk.itertuples():
        f=t[t.context.eq(row.context)&t.program.eq(row.program)&t.reference.eq(row.reference)];v=f[f.function_AB_rho26_plus_third.notna()];assert row.planned_subtype_pairs==len(f) and row.describable_subtype_pairs==len(v)
        assert np.allclose([getattr(row,'median_'+field) for field in fields],[v[field].median() for field in fields],atol=1e-12,rtol=1e-10,equal_nan=True)
    text=['# Complete third-lineage functional-background evidence (private)','',
          'All572 original subtype pairs x50 functions remain. Every before/after comparison uses the identical explicitly B-captured case subset. Full independent OLS/case/control/denominator audit precedes this report.',
          'Medians share author cases and are not independent replication or meta-analysis. No new P/BH/CI, direct ligand-receptor, causal mediation, protein activity or baseline prediction claim.',
          'B-lineage function can represent real shared ecology, technical/sample environment or a mediator. Attenuation is not automatic artifact; retention is not proof of pair specificity. B subtype mixture and incompletely reviewed patient identity remain limits.','',
          '## Complete artifacts','- private_all400_parent_function_background_summaries.tsv','- private_all8_parent_capture_summaries.tsv','- private_all_parent_functions_before_after_B_background.png','- All28600 original descriptions and572 capture rows remain in their immutable source directory.','',
          '## Interpretation','Preserve prior clinical/LR negative families. Connect only identical context/MCP/endpoint evidence. Assess third-lineage-shared functional programmes and remaining pair coordination before planning independent source/clinical/spatial validation; do not append same-source cutoffs to force significance.']
    (out/'PRIVATE_COMPLETE_THIRD_LINEAGE_REVIEW.md').write_text('\n'.join(text)+'\n')
    terminal=dict(status='COMPLETE_ALL400_PARENT_FUNCTION_AND_ALL8_CAPTURE_SUMMARIES_DATA_AND_RENDERED_BOUNDS_PASS_VISUAL_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),figure=bounds,files={p.name:sha(p) for p in out.iterdir() if p.is_file()},new_P_BH=False,actual_visual_review=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
