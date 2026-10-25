"""Add reviewed public design evidence without promoting local assay/clinical truth.

Preserves source denominators and original registry. No individual records,
counts or model inputs are read or changed by this metadata overlay.
"""
import argparse,datetime,json,hashlib,os,subprocess
from pathlib import Path


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(root):
    cfg=json.loads((root/'config/design_extension_spec_v1.json').read_text())
    for b in cfg['bindings'].values():
        if sha(b['path'])!=b['sha256']:raise ValueError('Source design evidence binding changed')
    roles=Path(cfg['source_roles']);terminal=json.loads((roles/'provenance/terminal_v1.json').read_text());file=roles/'outputs/source_roles_and_coverage_v1.json'
    if sha(file)!=terminal['output_sha256']:raise ValueError('Accepted source-role output changed')
    original=json.loads(file.read_text());evidence=json.loads((root/'config/public_design_extension_v1.json').read_text());lookup={r['source']:r for r in evidence['rows']};sources={r['source'] for r in original['sources']}
    if len(lookup)!=len(evidence['rows']) or set(lookup)-sources:raise ValueError('Public design/source coverage differs')
    records=[]
    for r in original['sources']:
        out={'source':r['source'],'n_archive_rows':r['n_archive_cells'],'n_expression_units_not_patients':r['n_archive_expression_unit_keys'],'native_count_origin_preserved':r['native_count_origin'],
             'original_capture_and_assay_decision':r['design_decision'],'original_analysis_roles':r['analysis_roles'],'original_cross_study_alias_status':r['cross_study_alias_status'],
             'current_conservative_analysis_roles':dict(r['analysis_roles']),'new_patient_stage_endpoint_qualification':False,'native_slot_full_length_or_UMI_provenance_not_promoted':True,'original_archive_model_counts_or_labels_changed':False}
        if r['source'] in lookup:
            e=lookup[r['source']];out['matched_public_design_evidence']=e;out['exact_local_asset_version_or_library_join_accepted']=False
            if r['source']=='GSE115978':
                out['current_conservative_analysis_roles']['whole_tissue_composition']='do_not_use_unweighted_merged_CD45_sorted_strata_as_native_tissue_proportions';out['current_conservative_analysis_roles']['within_selected_lineage_composition']='requires_sorted_fraction_and_full_length_assay_scope';out['capture_view_restriction_evidence']='official_study_CD45_FACS_design_not_a_new_patient_or_asset_acceptance'
            elif r['source']=='GSE132465':out['author_upstream_identity_discordance_filter_not_reconstructable_from_retained_archive']=True
            elif r['source']=='GSE183904':out['anatomical_mixture_requires_record_level_primary_gastric_vs_peritoneal_join']=True
        else:out['matched_public_design_evidence']=None
        records.append(out)
    if {r['source'] for r in records}!=sources or sum(r['n_archive_rows'] for r in records)!=sum(r['n_archive_cells'] for r in original['sources']):raise ValueError('Design overlay changed archive denominator')
    value={'schema':'source_design_qualification_extension_v1','sources':records,'n_sources':len(records),'n_sources_with_new_public_design_review':len(lookup),'source_and_patient_independence_not_newly_accepted':True,'composition_clinical_analysis_run':False,'native_expression_read':False,'frozen_benchmark_and_model_inputs_unchanged':True}
    output=root/'outputs/source_design_qualification_extension_v1.json'
    with output.open('x') as f:json.dump(value,f,indent=2,ensure_ascii=False)
    readback=json.loads(output.read_text())
    if readback!=value:raise ValueError('Design overlay readback mismatch')
    terminal={'status':'SOURCE_DESIGN_REVIEW_OVERLAY_COMPLETE_NO_PATIENT_ASSET_OR_CLINICAL_PROMOTION','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(records),'n_reviewed_public_design_sources':len(lookup),'n_archive_rows_preserved':sum(r['n_archive_rows'] for r in records),'output_sha256':sha(output),'original_source_role_output_sha256':sha(file),'readback_verified':True,'new_model_or_annotation':False,'native_expression_read':False}
    with (root/'provenance/terminal_v1.json').open('x') as f:json.dump(terminal,f,indent=2)
    print(json.dumps(terminal))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);main(Path(p.parse_args().output))
