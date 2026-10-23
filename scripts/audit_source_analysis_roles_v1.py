"""Source-level role/QC/identity coverage without exporting individual rows.

Existing reviewed unit identities are preserved. A native identifier is not
automatically a reviewed patient, and a tissue label is not a capture census.
"""
import argparse
import collections
import datetime
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
import pandas as pd

from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection

STATUS_FIELDS=['assay','suspension_type','capture_scope','analysis_role','prior_QC_status','QC_prior',
 'patient_identity_status','patient_identity_verified','clinical_validation_eligible','clinical_endpoint_definition_status',
 'source_alias_acceptance_status','preintegration_cell_eligible','preintegration_release_eligible','downstream_release_eligible',
 'residual_doublet_flag','ambient_status','ambient_high_loss_review','analysis_correction_applied',
 'annotation_exclude_from_primary','CM_composition_eligible','ex_vivo_stimulation']
AUTHOR_FIELDS=['source_author_label','author_cell_type','cell_type','ev_author_label']


def usable(values):
    s=np.asarray(values).astype(str)
    lower=np.char.lower(s)
    return ~np.isin(lower,['','nan','none','na','n/a','unknown','unreviewed','not_available','not_applicable','<na>'])


def unit_identity_coverage(units,patient,verified):
    unit_code,unit_values=pd.factorize(np.asarray(units).astype(str),sort=True)
    patient=np.asarray(patient).astype(str);valid=usable(patient)
    verify=np.isin(np.char.lower(np.asarray(verified).astype(str)),['true','1','verified','reviewed_identity'])
    if len(patient)!=len(unit_code) or len(verify)!=len(unit_code):raise ValueError('Unit/patient axes differ')
    n=len(unit_values);total=np.bincount(unit_code,minlength=n)
    identified=np.bincount(unit_code[valid],minlength=n)
    reviewed=np.bincount(unit_code[valid&verify],minlength=n)
    pc,pv=pd.factorize(patient[valid],sort=True)
    pairs=np.unique(np.column_stack([unit_code[valid],pc]),axis=0)
    distinct=np.bincount(pairs[:,0],minlength=n) if len(pairs) else np.zeros(n,int)
    consistent=(distinct==1)&(identified==total)
    reviewed_unit=consistent&(reviewed==total)
    return {'n_expression_units_not_patients':n,'n_units_no_usable_patient_key':int((identified==0).sum()),
            'n_units_multiple_patient_keys':int((distinct>1).sum()),'n_units_partial_patient_key_coverage':int(((identified>0)&(identified<total)).sum()),
            'n_units_complete_single_patient_key_unreviewed_or_reviewed':int(consistent.sum()),
            'n_units_complete_single_patient_key_and_all_rows_verified':int(reviewed_unit.sum()),
            'n_distinct_observed_patient_keys_not_verified_patient_count':len(pv),
            'native_key_presence_is_not_identity_review':True,'patient_phase_endpoint_qualification_still_separate':True}


def source_role(design):
    kind=design['capture_class']
    if kind in ['immune_enriched','T_enriched_nuclei']:
        return {'whole_tissue_composition':'not_eligible_capture_enriched',
                'within_selected_lineage_composition':'candidate_requires_unit_material_and_capture_comparability',
                'expression_reference':'candidate_with_capture_and_assay_scope'}
    if kind=='cell_product':
        return {'whole_tissue_composition':'not_eligible_cell_product','within_selected_lineage_composition':'product_specific_view',
                'expression_reference':'product_specific_reference'}
    if kind=='selected_nonmalignant_subset':
        return {'whole_tissue_composition':'not_eligible_selected_subset','within_selected_lineage_composition':'candidate_with_selection_scope',
                'expression_reference':'candidate_with_selection_scope'}
    if kind=='whole_tissue_sampling_candidate':
        return {'whole_tissue_composition':'candidate_capture_not_certified_by_tissue_name',
                'within_selected_lineage_composition':'candidate_requires_unit_review','expression_reference':'candidate_native_QC_and_assay_review'}
    if kind=='unknown':
        return {'whole_tissue_composition':'unqualified_capture_unknown','within_selected_lineage_composition':'exploratory_source_conditioned_only',
                'expression_reference':'candidate_native_QC_and_assay_review'}
    raise ValueError('Unreviewed capture class')


def main(root):
    started=time.monotonic();cfg=json.loads((root/'config/source_role_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path'])!=binding['sha256']:raise ValueError('Source-role frozen input changed')
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_source_analysis_roles_v1.py')],capture_output=True,text=True,timeout=120)
    (root/'provenance/test_stdout.txt').write_text(test.stdout);(root/'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode!=0 or 'Ran 4 tests' not in test.stderr or 'skipped' in test.stderr.lower():raise ValueError('Source-role counterexamples failed')
    reader=NativeArchiveCollection(cfg['archive_manifest'])
    bridge=json.loads(Path(cfg['source_bridge']).read_text());catalog={r['dataset']:r for r in bridge['source_registry']}
    if set(catalog)!=set(reader.partitions):raise ValueError('Source role coverage differs from selected archive')
    records=[];legacy_reviewed=collections.Counter()
    for ds,part in reader.archive_partitions.items():
        original=catalog[ds];units=reader.frozen_expression_units(ds)
        n=len(reader.rows(ds));nunits=len(np.unique(units))
        if n!=part['n_archive_rows'] or nunits!=part['n_expression_unit_keys']:raise ValueError('Archive source denominator changed')
        decision=cfg['source_decisions'].get(ds,{'capture_class':'unknown','assay_branch':'unreviewed','basis':'frozen_registry_capture_unreviewed'})
        row={'source':ds,'provider':original['provider'],'n_archive_cells':n,'n_frozen_context_cells':part['n_frozen_context_rows'],
             'n_archive_expression_unit_keys':nunits,'n_frozen_context_expression_units':original['n_expression_units'],
             'design_decision':decision,'analysis_roles':source_role(decision),'catalog_design':original['catalog_design'],
             'native_count_origin':part['count_origin'],'native_cell_key_contract':part['cell_key_contract'],
             'archive_membership_is_not_model_or_clinical_eligibility':True,
             'new_clinical_view_status':'requires_qualified_patient_phase_endpoint_join_prior_accepted_analyses_preserved',
             'cross_study_alias_status':'unresolved_unless_explicit_overlap_hold',
             'known_possible_overlap_group':cfg['possible_overlap_groups'].get(ds),
             'possible_overlap_group_is_not_proven_duplicate':True}
        if reader.is_legacy(ds):
            meta=original['frozen_expression_unit_metadata'];row['frozen_expression_unit_metadata']=meta
            if any(sum(v.values())!=original['n_expression_units'] for v in meta.values()):raise ValueError('Legacy source metadata denominator differs')
            legacy_reviewed.update(meta['patient_review_status'])
            row['patient_identity_evidence']='frozen_review_status_preserved_not_overwritten_by_catalog_defaults'
            row['archive_units_without_frozen_context_metadata']=nunits-original['n_expression_units']
            row['native_QC_metadata_not_reaudited_here']=True
        else:
            names=set(reader.field_names(ds));row['native_status_cell_counts']={};row['author_label_coverage']={}
            for field in STATUS_FIELDS:
                if field not in names:
                    row['native_status_cell_counts'][field]={'status':'field_absent'};continue
                values=reader.field(ds,field).astype(str)
                counts={str(k):int(v) for k,v in collections.Counter(values).items()}
                if sum(counts.values())!=n:raise ValueError('Native status denominator differs')
                row['native_status_cell_counts'][field]={'status':'observed_field_not_automatically_accepted','counts':counts}
                reader._fields.clear()
            for field in AUTHOR_FIELDS:
                if field not in names:continue
                values=reader.field(ds,field).astype(str);ok=usable(values)
                row['author_label_coverage'][field]={'n_nonmissing_cells':int(ok.sum()),'n_distinct_nonmissing_labels':len(np.unique(values[ok])),
                    'not_independent_truth':True,'label_mapping_and_original_paper_provenance_still_required':True}
                reader._fields.clear()
            identity_fields=[p for p in ['patient_id','patient_uid','donor_id'] if p in names]
            row['native_patient_key_coverage']={}
            verified=reader.field(ds,'patient_identity_verified') if 'patient_identity_verified' in names else np.full(n,'unknown')
            reader._fields.clear()
            for field in identity_fields:
                values=reader.field(ds,field)
                row['native_patient_key_coverage'][field]=unit_identity_coverage(units,values,verified)
                reader._fields.clear()
            row['frozen_collection_design']=original['frozen_collection_design']
            row['patient_identity_evidence']='native_flags_and_key_consistency_audited_no_new_patient_identity_claim'
        records.append(row)
        save(root/'checkpoints'/(ds+'.json'),{'source':ds,'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    if sum(r['n_frozen_context_expression_units'] for r in records)!=sum(r['n_expression_units'] for r in catalog.values()):raise ValueError('Context unit scope differs')
    payload={'sources':records,'selected_source_scope_only':True,'held_whole_sources_not_silently_reintroduced':True,
             'legacy_review_status_units':dict(legacy_reviewed),'individual_rows_exported':False,
             'source_keys_are_not_independent_studies':True,'source_key_uniqueness_does_not_close_alias_review':True,
             'model_training_does_not_use_clinical_outcomes':True,'new_composition_or_clinical_analysis_run':False}
    out=root/'outputs/source_roles_and_coverage_v1.json';save(out,payload)
    size=sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    if size>cfg['output_cap_bytes']:raise ValueError('Source role output cap exceeded')
    save(root/'provenance/terminal_v1.json',{'status':'SOURCE_ROLE_AND_NATIVE_METADATA_COVERAGE_COMPLETE_NOT_NEW_CLINICAL_QUALIFICATION',
         'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources':len(records),'n_tests_actual_no_skips':4,
         'output_sha256':sha256_file(out),'wall_seconds':time.monotonic()-started,
         'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'written_bytes_before_terminal':size,
         'individual_rows_exported':False,'native_expression_read':False,'new_clinical_analysis':False})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root/'provenance/failure_v1.json',{'error_class':type(error).__name__,'error':str(error)})
        raise
