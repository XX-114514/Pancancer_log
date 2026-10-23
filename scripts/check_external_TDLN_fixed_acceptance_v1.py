"""Artificial aggregate fixtures only; no native or clinical records."""
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import math
import tempfile
import accept_external_TDLN_fixed_spatial_v1 as consumer


def artificial_rows():
    rows = []
    for endpoint in consumer.ENDS:
        for estimand in consumer.DIRECT:
            absent = endpoint == consumer.ENDS[2]
            zero = estimand == 'cycling_minus_noncycling'
            negative = endpoint == consumer.ENDS[1]
            effect = None if absent else (0 if zero else (-1 if negative else 1))
            row = dict(endpoint=endpoint, estimand=estimand, cases=0 if absent else 7,
                       eligible_images=0 if absent else 14,
                       positive_cases=0 if absent or zero or negative else 7,
                       negative_cases=0 if absent or zero or not negative else 7,
                       zero_cases=7 if zero and not absent else 0,
                       median_effect=effect, bootstrap_lower=effect, bootstrap_upper=effect,
                       exact_sign_p=None if absent else (1 if zero else .015625),
                       primary_p=None if absent else (1 if zero else .015625),
                       status='STRUCTURAL_MARKER_ABSENCE' if absent else 'ESTIMABLE')
            rows.append(row)
    for endpoint in consumer.ENDS:
        absent = endpoint == consumer.ENDS[2]
        negative = endpoint == consumer.ENDS[1]
        effect = None if absent else (-.5 if negative else .5)
        rows.append(dict(endpoint=endpoint, estimand=consumer.LOCAL, grid=4,
                         cases=0 if absent else 7, eligible_images=0 if absent else 14,
                         positive_cases=0 if absent or negative else 7,
                         negative_cases=7 if negative else 0, zero_cases=0,
                         median_effect=effect, bootstrap_lower=effect, bootstrap_upper=effect,
                         permutations=999, lower_tail_count=0 if absent or negative else 999,
                         upper_tail_count=999 if negative else 0, source_null_varies=not absent,
                         permutation_p=None if absent else .002,
                         primary_p=None if absent else .002,
                         status='STRUCTURAL_MARKER_ABSENCE' if absent else 'ESTIMABLE'))
    # Hand-computed oracle: two .002, six .015625, two 1, five missing slots.
    harmonic = 3.3182289932289937
    for i, row in enumerate(rows):
        q = None if row['primary_p'] is None else (1 if row['primary_p'] == 1 else
                                                  (.015 if i >= 12 else .029296875))
        row['BH'] = q
        row['BY'] = None if q is None else min(1, q*harmonic)
    return rows


def artificial_summary():
    unit = dict(registered_case_ROI_partition_gate_passed=True,
                literal_case_labels_not_accepted_independent_patients=True,
                alias_checks=[dict(full_primary_coverage_and_same_partition=True)]*2,
                duplicate_case_ROI_object_rows=0, invalid_object_rows=0,
                missing_rows_by_field=dict(case=0, roi=0, image=0, sample=0))
    states = []
    for state, markers in consumer.STATE_MARKERS.items():
        absent = state == consumer.ENDS[2]
        states.append(dict(state=state, physical_markers=list(markers),
                           fixed_relative_cutoffs=consumer.STATE_CUTOFFS[state],
                           marker_set_available=not absent,
                           measurement_status='STRUCTURAL_MARKER_ABSENCE' if absent else
                                              'RELATIVE_STATE_MEASUREMENT_AVAILABLE',
                           cell_rows=None if absent else 70))
    phenotype = dict(parent_selection_excludes_target_state_markers=True,
                     parent_agreement_is_case_balanced_across_images=True,
                     missing_or_unmeasured_states_retained_as_NA=True,
                     state_records=states,resolved_parent_group_counts=dict(CD4=1,CD8a=1))
    spatial = dict(planned_primary_tests=15,structural_missing_endpoints_are_NA_not_negative=True,
                   literal_case_keys_not_accepted_independent_patients=True,
                   primary_rows=artificial_rows(),sensitivities=[
                       dict(endpoint=e,grid=g,permutation_p=None) for e in consumer.ENDS for g in (2,8)])
    audit = dict(status='PASS_FULL_EXTERNAL_FIXED_SPATIAL_NUMERIC_AUDIT',maximum_numeric_error=0,
                 maximum_exact_distance_error=0,direct_tables_all_recalculated=True,
                 all_registered_case_ROI_aggregation_bootstrap_tails_and_BH_recalculated=True)
    return dict(status='COMPLETE_EXTERNAL_FIXED_STATE_SPATIAL_EXTENSION_AUDITED',unit=unit,
                clinical_columns_not_selected=True,author_dyad_or_distance_columns_not_selected=True,
                independent_rank_masks_exactly_match=True,phenotype=phenotype,assay_sensitivity={},
                registered_case_keys_not_independent_patients_or_clinical_validation=True,
                relative_construct_not_exact_RNA_CM_projection_or_direct_LR_mechanism=True,
                spatial=spatial,audit=audit)


def rejected(call):
    try:
        call()
    except (AssertionError, ValueError):
        return True
    raise AssertionError('ARTIFICIAL_COUNTEREXAMPLE_WAS_NOT_REJECTED')


def check():
    analysis = dict(physical_panel=sorted(set(m for k,v in consumer.STATE_MARKERS.items()
                                            if k != consumer.ENDS[2] for m in v)),
                    spatial=dict(minimum_cases=5,primary_grid=4,permutations=999))
    cfg = analysis['spatial'];available = dict(zip(consumer.ENDS,[True,True,False]))
    summary = artificial_summary()
    good = consumer.validate_summary(summary,analysis)
    assert good['Python_aggregate_audit']['estimable_primary_slots'] == 10
    assert all(abs(a-b)<1e-14 for a,b in zip(
        consumer.adjusted_family([.001,.01,.04]+[None]*12)[:3],[.015,.075,.2]))
    assert consumer.direction_p(32,0,0,5) == 2**-31
    assert consumer.direction_p(0,0,7,5) == 1 and consumer.direction_p(1,0,0,5) is None
    changed = deepcopy(summary);changed['spatial']['primary_rows'].pop()
    rejected(lambda:consumer.validate_summary(changed,analysis))
    changed = deepcopy(summary);changed['spatial']['primary_rows'][0]['BH'] /= 2
    rejected(lambda:consumer.validate_summary(changed,analysis))
    changed = deepcopy(summary);changed['spatial']['primary_rows'][8]['median_effect'] = 0
    rejected(lambda:consumer.validate_summary(changed,analysis))
    changed = deepcopy(summary);changed['phenotype']['state_records'][0]['fixed_relative_cutoffs'] = [.5,.5]
    rejected(lambda:consumer.validate_summary(changed,analysis))
    changed = deepcopy(summary);changed['unit']['patient_ids'] = ['SYNTHETIC_CANARY']
    rejected(lambda:consumer.validate_summary(changed,analysis))
    changed = deepcopy(summary);changed['phenotype']['state_records'][0]['cell_rows'] = [1,2]
    rejected(lambda:consumer.validate_summary(changed,analysis))
    unitfail = dict(status='EXTERNAL_UNIT_GATE_NOT_PASSED_NO_SPATIAL_EFFECTS',
                    clinical_columns_not_selected=True,unit=dict(
                        registered_case_ROI_partition_gate_passed=False,failure_reasons=['ARTIFICIAL_ALIAS_CONFLICT']))
    assert 'UNIT_GATE_FAILURE' in consumer.validate_summary(unitfail,analysis)['status']
    parentfail = deepcopy(summary);parentfail['status'] = 'EXTERNAL_PARENT_GATE_NOT_PASSED_NO_SPATIAL_EFFECTS'
    parentfail.pop('spatial');parentfail['phenotype']['resolved_parent_group_counts']['CD4'] = 0
    assert 'PARENT_GATE_FAILURE' in consumer.validate_summary(parentfail,analysis)['status']
    # Preserve the small artificial provenance tree; no native fixture files.
    task = Path(tempfile.mkdtemp(prefix='pancancer_artificial_fixed_acceptance_'))
    code = task/'generic_code';code.mkdir();(code/'a.R').write_text('# ARTIFICIAL ONLY\n')
    def write(name,value):
        (task/name).write_text(json.dumps(value))
    assert consumer.validate_completed(task,code) is None
    plan = dict(frozen_UTC='2026-01-01T00:00:00+00:00',expected_file_bytes=8,
                expected_file_MD5='a'*32,expected_file_SHA256='b'*64,range_chunk_bytes=2,
                maximum_expanded_bytes=32,analysis=analysis,R_components=['a.R'],
                code_SHA256={'a.R':consumer.digest(code/'a.R')},R_code_SHA256=consumer.digest(code/'a.R'))
    proof = dict(complete_ranges=4,exact_ranges_and_lengths=True,no_gaps_or_overlaps=True,
                 whole_official_MD5_verified=True,no_range_payload_files=True,
                 no_alternate_endpoint_or_authentication_protocol_hash_retries=True)
    started = dict(started_UTC='2026-01-01T00:01:00+00:00')
    identity = dict(bytes=8,MD5='a'*32,SHA256='b'*64,range_assembly_proof=proof,
                    observed_UTC='2026-01-01T00:02:00+00:00')
    source = dict(kind='COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED',started_UTC=started['started_UTC'],
                  source_R_exit_code=0,source_payload_RAM_released=True,no_dynamic_followups_or_native_cache=True,
                  no_native_individual_clinical_coordinate_state_or_graph_files=True,
                  no_individual_values_or_identifiers_transferred=True,input_bytes=8,input_MD5='a'*32,
                  input_SHA256='b'*64,R_code_SHA256=plan['R_code_SHA256'],range_assembly_proof=proof,
                  pipe_technical_proof=dict(expanded_stream_bytes=16,no_native_or_expanded_payload_files=True),
                  completed_UTC='2026-01-01T00:03:00+00:00',summary=summary)
    write('private_frozen_fixed_native_workflow_v1.json',plan)
    write('private_complete_input_identity_v1.json',identity)
    write('private_actual_source_start_v1.json',started)
    write('private_complete_source_fixed_workflow_v1.json',source)
    exited = dict(status='SOURCE_EXIT0_FIXED_WORKFLOW_AND_RAM_RELEASE_COMPLETE',actual_exit_code=0,
                  frozen_contract_SHA256=consumer.digest(task/'private_frozen_fixed_native_workflow_v1.json'),
                  complete_source_SHA256=consumer.digest(task/'private_complete_source_fixed_workflow_v1.json'),
                  observed_UTC='2026-01-01T00:04:00+00:00')
    write('private_actual_exit_v1.json',exited)
    assert consumer.validate_completed(task,code)['status'] == 'ACCEPTED_COMPLETED_SOURCE_WORKFLOW_BINDING'
    plan['input_transport'] = 'bounded_local_opaque_ranges_over_existing_SSH'
    write('private_frozen_fixed_native_workflow_v1.json',plan)
    source['range_assembly_proof'].update(ordered_stream_no_local_file_sized_buffer=True,
        opaque_compressed_relay_over_existing_SSH=True,
        whole_bytes_independently_hashed_on_source_before_R=True,
        range_header_proof_from_frozen_local_streamer=True)
    identity['range_assembly_proof'] = source['range_assembly_proof']
    write('private_complete_input_identity_v1.json',identity)
    write('private_complete_source_fixed_workflow_v1.json',source)
    exited['frozen_contract_SHA256'] = consumer.digest(task/'private_frozen_fixed_native_workflow_v1.json')
    exited['complete_source_SHA256'] = consumer.digest(task/'private_complete_source_fixed_workflow_v1.json')
    write('private_actual_exit_v1.json',exited)
    rejected(lambda:consumer.validate_completed(task,code))
    local = dict(kind='COMPLETE_ORDERED_PUBLIC_INPUT_PROOF',MD5='a'*32,SHA256='b'*64,
                 range_assembly_proof=source['range_assembly_proof'])
    write('private_complete_local_opaque_stream_v1.json',local)
    assert consumer.validate_completed(task,code)['local_complete_stream_proof_SHA256'] is not None
    local['SHA256'] = 'c'*64;write('private_complete_local_opaque_stream_v1.json',local)
    rejected(lambda:consumer.validate_completed(task,code))
    local['SHA256'] = 'b'*64;write('private_complete_local_opaque_stream_v1.json',local)
    changed = deepcopy(exited);changed['actual_exit_code'] = 1
    write('private_actual_exit_v1.json',changed)
    rejected(lambda:consumer.validate_completed(task,code))
    write('private_actual_exit_v1.json',exited)
    (code/'a.R').write_text('# ARTIFICIAL ALTERED\n')
    rejected(lambda:consumer.validate_completed(task,code))
    write('private_actual_failure_exit_v1.json',dict(status='ARTIFICIAL_FAILURE'))
    rejected(lambda:consumer.validate_completed(task,code))
    return dict(status='PASS_ARTIFICIAL_FIXED_SPATIAL_ACCEPTANCE_AND_COUNTEREXAMPLES',
                source_exit_input_and_code_binding_checked=True,
                known_direction_and_15_slot_BH_BY_oracles_checked=True,
                missing_not_zero_and_changed_cutoff_and_individual_payload_rejected=True,
                unit_parent_failure_and_pending_not_misreported_as_spatial_evidence=True,
                no_real_native_clinical_or_individual_data_read=True,
                observed_UTC=datetime.now(timezone.utc).isoformat())


if __name__ == '__main__':
    print(json.dumps(check()))
