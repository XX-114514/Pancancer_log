"""Independent full evidence joins, scalar classes and forest coverage audit."""
import hashlib,json,math,os,time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from PIL import Image,ImageStat

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8388608),b''):h.update(b)
    return h.hexdigest()

def same(a,b):assert np.allclose(a,b,atol=1e-12,rtol=1e-12,equal_nan=True),(a,b)

def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);O=R/'NSCLC_automatic_ligand_response_review_20261016_v1'
    launch=json.loads((R/'HNSCC_receiver_followup_job_control_20261016_v2/review_NSCLC_automatic_ligand_response_v1/launch_receipt.json').read_text());deadline=time.monotonic()+110000
    while not (O/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists();s=p.read_text().rsplit(')',1)[1].split()
        assert s[19]==str(launch['proc_start_ticks']) and s[0]!='Z' and time.monotonic()<deadline;time.sleep(20)
    output=O/'private_independent_audit.json';assert not output.exists();term=json.loads((O/'private_terminal.json').read_text())
    for name,h in term['files'].items():assert sha(O/name)==h
    I=R/'NSCLC_automatic_ligand_response_results_20261016_v1';ia=I/'private_independent_audit.json';assert sha(ia)==term['result_audit_sha256']
    original=pd.read_csv(I/'private_all_unique_ligand_response_results.tsv',sep='\t').set_index('record')
    data=pd.read_csv(O/'private_all44528_annotated_evidence.tsv.gz',sep='\t').set_index('record');assert original.index.tolist()==data.index.tolist()
    for col in original.columns:
        if pd.api.types.is_numeric_dtype(original[col]):same(data[col],original[col])
        else:assert data[col].fillna('').tolist()==original[col].fillna('').tolist()
    C=R/'NSCLC_frozen_signature_clinical_results_20261016_v1/private_all_planned_signature_clinical_results.tsv';assert sha(C)==term['clinical_result_sha256']
    clinical=pd.read_csv(C,sep='\t').set_index(['context','program','type','sub_cell_type'])
    Q=R/'NSCLC_two_end_coordination_results_20261016_v1/private_all572_two_end_results.tsv';assert sha(Q)==term['coordination_result_sha256'];coord=pd.read_csv(Q,sep='\t').set_index('record')
    V=Path(json.loads((R/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    P=V/'original_CM_programme_pair_associations_20261015_v2/all_context_coverage_reasons.tsv';assert sha(P)==term['training_context_sha256'];contexts=pd.read_csv(P,sep='\t').set_index('context')
    categories={};candidate_ids=[];bh_ids=[]
    for row in data.reset_index().itertuples(index=False):
        if not math.isfinite(float(row.sequential_wild_p)):category='NOT_ESTIMABLE'
        elif row.BH_q_complete_unique_family>=.05:category='INSUFFICIENT_COMPLETE_FAMILY_EVIDENCE'
        else:
            bh_ids.append(row.record);interval=math.isfinite(row.CI_low) and math.isfinite(row.CI_high) and (row.CI_low>0 or row.CI_high<0)
            stable=row.LOO_estimable>=.95*row.cases and row.LOO_same_sign>=.95*row.LOO_estimable
            category='SOURCE_LIMITED_STATISTICAL_CANDIDATE' if interval and stable else 'BH_SUPPORTED_ROBUSTNESS_INCOMPLETE'
        assert row.evidence_class==category;categories[category]=categories.get(category,0)+1
        if category=='SOURCE_LIMITED_STATISTICAL_CANDIDATE':candidate_ids.append(row.record)
        context=contexts.loc[row.context];assert row.training_gse==context.gse_id and row.training_stratum==context.stratum
        expected=';'.join('2026年10月'+cm+'（v20261010）' for cm in context.CMs.split(';') if cm);assert row.original_CM_provenance==expected
        rt=row.receiver_type;st='typeB' if rt=='typeA' else 'typeA'
        assert row.original_receiver_endpoint==(context.endpoint_A if rt=='typeA' else context.endpoint_B)
        assert row.original_sender_endpoint==(context.endpoint_B if st=='typeB' else context.endpoint_A)
        receiver=clinical.loc[(row.context,row.program,rt,row.receiver_author_subtype)];sender=clinical.loc[(row.context,row.program,st,row.sender_author_subtype)]
        same([row.full_receiver_signature_clinical_beta,row.full_receiver_signature_clinical_q],[receiver.beta_MPR_minus_nonMPR,receiver.BH_q_full_planned_family])
        same([row.full_sender_signature_clinical_beta,row.full_sender_signature_clinical_q],[sender.beta_MPR_minus_nonMPR,sender.BH_q_full_planned_family])
        q=coord.loc[row.parent_pair];same([row.original_two_end_coordination_rho,row.original_two_end_coordination_q],[q.partial_rank_rho,q.BH_q_full572])
        assert row.descriptive_combination_not_joint_FDR and not row.native_author_subtype_to_old_node_mapping_validated
    supported=pd.read_csv(O/'private_all_source_limited_statistical_candidates.tsv',sep='\t');bh=pd.read_csv(O/'private_all_primary_BH_supported_rows.tsv',sep='\t')
    assert supported.record.tolist()==candidate_ids and bh.record.tolist()==bh_ids
    summary=pd.read_csv(O/'private_all16_context_direction_summaries.tsv',sep='\t');assert len(summary)==16
    for r in summary.itertuples():
        f=data[data.context.eq(r.context)&data.program.eq(r.program)&data.receiver_type.eq(r.receiver_type)]
        assert r.planned==len(f) and r.evaluated==f.sequential_wild_p.notna().sum()
        assert r.primary_BH_rows==f.BH_q_complete_unique_family.lt(.05).sum() and r.BY_rows==f.BY_q_dependence_sensitivity.lt(.05).sum()
        assert r.robust_rows==f.evidence_class.eq('SOURCE_LIMITED_STATISTICAL_CANDIDATE').sum()
    expected=data[data.sequential_wild_p.notna()].reset_index().sort_values(['context','program','receiver_type','record'])
    forest=pd.read_csv(O/'private_all_estimable_forest_data.tsv',sep='\t');assert forest.record.tolist()==expected.record.tolist()
    for col in ['partial_rank_rho','CI_low','CI_high','BH_q_complete_unique_family','BY_q_dependence_sensitivity','raw_spearman','sender_programme_adjusted_rho','pathology_adjusted_rho','cases']:
        same(forest[col],expected[col])
    pages=json.loads((O/'private_all_forest_page_records.json').read_text());assert [r for p in pages for r in p['records']]==expected.record.tolist()
    assert len(pages)==math.ceil(len(expected)/40) and len({r for p in pages for r in p['records']})==len(expected)
    for name in ['private_coverage_and_adjustment.png',*[p['file'] for p in pages]]:
        with Image.open(O/name) as im:
            assert im.width>=1000 and im.height>=800 and ImageStat.Stat(im.convert('L')).var[0]>10
    result=dict(status='PASS_ALL44528_VALUES_PROVENANCE_CLINICAL_COORDINATION_JOINS_SCALAR_CLASSES_ALL_SUMMARIES_AND_FULL_FOREST_DATA_COVERAGE',
        finished_UTC=datetime.now(timezone.utc).isoformat(),audit_code_sha256=sha(__file__),producer_terminal_sha256=sha(O/'private_terminal.json'),
        all_evidence_rows=len(data),all_estimated_forest_rows=len(expected),forest_pages=len(pages),evidence_class_counts=categories,
        all_pages_individually_visually_reviewed=False,new_tests_or_P_or_BH=False,receptor_specific_or_causal_axis_validated=False)
    output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':main()
