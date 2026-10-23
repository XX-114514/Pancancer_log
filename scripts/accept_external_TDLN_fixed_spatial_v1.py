"""Accept a completed source-aggregate workflow, including honest gate failures.

No native payload, individual records or arrays are read. Binding checks cover
producer proofs, not an independent local rehash of the remote native payload.
Exact direction tails and the full 15-slot BH/BY family are checked in Python.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import math
import os


ENDS = ('E0_PD1_TIM3', 'E1_CD8_MKI67', 'E2_CCR7_CD127')
DIRECT = ('noncycling', 'cycling', 'combined_MH', 'cycling_minus_noncycling')
LOCAL = 'local_composition_degree_and_cycle_conditional_null'
STATE_MARKERS = {
    'R0_FOXP3_CTLA4': ('FoxP3', 'CTLA-4'),
    'R1_FOXP3_CTLA4_MKI67': ('FoxP3', 'CTLA-4', 'Ki-67'),
    'R2_general_CD4_MKI67': ('Ki-67',),
    'E0_PD1_TIM3': ('PD-1', 'TIM-3'),
    'E1_CD8_MKI67': ('Ki-67',),
    'E2_CCR7_CD127': ('CCR7', 'CD127'),
}
STATE_CUTOFFS = {
    'R0_FOXP3_CTLA4': [.75, .60], 'R1_FOXP3_CTLA4_MKI67': [.75, .60, .75],
    'R2_general_CD4_MKI67': [.75], 'E0_PD1_TIM3': [.70, .70],
    'E1_CD8_MKI67': [.75], 'E2_CCR7_CD127': [.70, .70],
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def integer(value):
    return type(value) is int and value >= 0


def same(a, b, tolerance=1e-10):
    if a is None or b is None:
        return a is None and b is None
    return finite(a) and finite(b) and abs(a-b) <= tolerance


def adjusted_family(values, by=False):
    """Missing positions stay in the denominator and are displayed missing."""
    m = len(values)
    ranked = sorted(range(m), key=lambda i: values[i] if values[i] is not None else 1)
    harmonic = sum(1/j for j in range(1, m+1)) if by else 1
    running, out = 1.0, [None]*m
    for rank in range(m, 0, -1):
        i = ranked[rank-1]
        p = 1 if values[i] is None else values[i]
        running = min(running, p*m*harmonic/rank)
        if values[i] is not None:
            out[i] = running
    return out


def direction_p(pos, neg, zero, minimum):
    n, d = pos+neg+zero, pos+neg
    if n < minimum:
        return None
    if d == 0:
        return 1.0
    return min(1.0, 2*sum(math.comb(d, j) for j in range(min(pos, neg)+1))/2**d)


def validate_rows(rows, spatial_config, available):
    expected = [(e, d) for e in ENDS for d in DIRECT]+[(e, LOCAL) for e in ENDS]
    assert len(rows) == 15
    assert [(v['endpoint'], v['estimand']) for v in rows] == expected
    probabilities = []
    for index, row in enumerate(rows):
        counts = [row[k] for k in ('cases', 'positive_cases', 'negative_cases', 'zero_cases')]
        assert all(integer(n) for n in counts)
        n, pos, neg, zero = counts
        assert pos+neg+zero == n
        assert integer(row['eligible_images']) and n <= row['eligible_images']
        if n == 0:
            assert row['median_effect'] is None
            assert row['bootstrap_lower'] is None and row['bootstrap_upper'] is None
        else:
            assert finite(row['median_effect'])
            assert finite(row['bootstrap_lower']) and finite(row['bootstrap_upper'])
            assert row['bootstrap_lower'] <= row['bootstrap_upper']
        if index < 12:
            p = direction_p(pos, neg, zero, spatial_config['minimum_cases'])
            assert same(p, row['exact_sign_p'])
        else:
            assert row['grid'] == spatial_config['primary_grid']
            assert row['permutations'] == spatial_config['permutations']
            low, high = row['lower_tail_count'], row['upper_tail_count']
            assert integer(low) and integer(high)
            assert low <= row['permutations'] and high <= row['permutations']
            infer = n >= spatial_config['minimum_cases'] and row['source_null_varies'] is True
            p = min(1.0, 2*(min(low, high)+1)/(row['permutations']+1)) if infer else None
            assert same(p, row['permutation_p'])
        assert same(p, row['primary_p'])
        if not available[row['endpoint']]:
            assert row['status'] == 'STRUCTURAL_MARKER_ABSENCE'
            assert n == 0 and row['eligible_images'] == 0
            assert row['median_effect'] is None and p is None
        probabilities.append(p)
    for method in ('BH', 'BY'):
        expected_q = adjusted_family(probabilities, by=method == 'BY')
        assert all(same(q, row[method]) for q, row in zip(expected_q, rows)), 'FAMILY_ADJUSTMENT_DIFFERENCE'
    return dict(status='PASS_ALL_FIXED_ROWS_TAILS_AND_FULL_FAMILY_ADJUSTMENT',
                fixed_primary_slots=len(rows), estimable_primary_slots=sum(v is not None for v in probabilities),
                structurally_unmeasured_slots=sum(not available[r['endpoint']] for r in rows),
                individual_bootstrap_and_geometry_audit_remains_source_bound=True)


def validate_summary(summary, analysis):
    reject_individual_payload(summary)
    assert summary['clinical_columns_not_selected'] is True
    unit = summary['unit']
    if summary['status'] == 'EXTERNAL_UNIT_GATE_NOT_PASSED_NO_SPATIAL_EFFECTS':
        assert unit['registered_case_ROI_partition_gate_passed'] is False
        assert unit['failure_reasons'] and 'spatial' not in summary and 'phenotype' not in summary
        return dict(status='COMPLETED_WITH_UNIT_GATE_FAILURE_NO_SPATIAL_EVIDENCE', unit=unit)
    assert unit['registered_case_ROI_partition_gate_passed'] is True
    assert unit['literal_case_labels_not_accepted_independent_patients'] is True
    assert all(v['full_primary_coverage_and_same_partition'] is True for v in unit['alias_checks'])
    assert unit['duplicate_case_ROI_object_rows'] == unit['invalid_object_rows'] == 0
    assert summary['independent_rank_masks_exactly_match'] is True
    assert summary['author_dyad_or_distance_columns_not_selected'] is True
    phenotype = summary['phenotype']
    assert phenotype['parent_selection_excludes_target_state_markers'] is True
    assert phenotype['parent_agreement_is_case_balanced_across_images'] is True
    assert phenotype['missing_or_unmeasured_states_retained_as_NA'] is True
    states = {v['state']: v for v in phenotype['state_records']}
    assert set(states) == set(STATE_MARKERS)
    available = {}
    for name, markers in STATE_MARKERS.items():
        row = states[name]
        assert row['physical_markers'] == list(markers)
        assert row['fixed_relative_cutoffs'] == STATE_CUTOFFS[name]
        measured = set(markers) <= set(analysis['physical_panel'])
        assert row['marker_set_available'] is measured
        if not measured:
            assert row['measurement_status'] == 'STRUCTURAL_MARKER_ABSENCE'
            assert row['cell_rows'] is None
        if row['measurement_status'] != 'RELATIVE_STATE_MEASUREMENT_AVAILABLE':
            assert row['cell_rows'] is None
        available[name] = measured
    if summary['status'] == 'EXTERNAL_PARENT_GATE_NOT_PASSED_NO_SPATIAL_EFFECTS':
        counts = phenotype['resolved_parent_group_counts']
        assert counts['CD4'] == 0 or counts['CD8a'] == 0
        assert 'spatial' not in summary
        return dict(status='COMPLETED_WITH_PARENT_GATE_FAILURE_NO_SPATIAL_EVIDENCE', unit=unit,
                    phenotype=phenotype, assay_sensitivity=summary['assay_sensitivity'])
    assert summary['status'] == 'COMPLETE_EXTERNAL_FIXED_STATE_SPATIAL_EXTENSION_AUDITED'
    assert summary['registered_case_keys_not_independent_patients_or_clinical_validation'] is True
    assert summary['relative_construct_not_exact_RNA_CM_projection_or_direct_LR_mechanism'] is True
    spatial, audit = summary['spatial'], summary['audit']
    assert spatial['planned_primary_tests'] == 15
    assert spatial['structural_missing_endpoints_are_NA_not_negative'] is True
    assert spatial['literal_case_keys_not_accepted_independent_patients'] is True
    assert audit['status'] == 'PASS_FULL_EXTERNAL_FIXED_SPATIAL_NUMERIC_AUDIT'
    assert finite(audit['maximum_numeric_error']) and audit['maximum_numeric_error'] < 1e-8
    assert finite(audit['maximum_exact_distance_error']) and audit['maximum_exact_distance_error'] < 1e-7
    assert audit['direct_tables_all_recalculated'] is True
    assert audit['all_registered_case_ROI_aggregation_bootstrap_tails_and_BH_recalculated'] is True
    numeric = validate_rows(spatial['primary_rows'], analysis['spatial'], available)
    sensitivity = spatial['sensitivities']
    assert len(sensitivity) == 6
    assert {(r['endpoint'], r['grid']) for r in sensitivity} == {(e, g) for e in ENDS for g in (2, 8)}
    assert all(r['permutation_p'] is None and 'primary_p' not in r for r in sensitivity)
    return dict(status='ACCEPTED_BOUND_SOURCE_STATE_SPATIAL_AGGREGATES',
                unit=unit, phenotype=phenotype, assay_sensitivity=summary['assay_sensitivity'],
                spatial=spatial, source_numeric_audit=audit, Python_aggregate_audit=numeric,
                not_exact_RNA_CM_projection_or_direct_communication_or_clinical_validation=True)


def reject_individual_payload(value):
    forbidden = {'internal', 'case_values', 'case_matrix', 'patient_ids', 'sample_ids', 'cell_ids',
                 'Patient.ID', 'ROI.ID', 'Image.name', 'sample_id', 'ObjectNumber',
                 'pg_clusters_corrected_40', 'coordinates', 'neighbors', 'bootstrap_indices',
                 'Group', 'Group.location', 'RFS.time'}
    if isinstance(value, dict):
        assert not forbidden.intersection(value), 'UNEXPECTED_INDIVIDUAL_OR_CLINICAL_PAYLOAD'
        for key in ('fine', 'image', 'case'):
            if key in value:
                assert value[key] is None or integer(value[key]), 'UNEXPECTED_IDENTIFIER_VECTOR'
        if 'cell_rows' in value:
            assert value['cell_rows'] is None or integer(value['cell_rows'])
        for child in value.values():
            reject_individual_payload(child)
    elif isinstance(value, list):
        assert len(value) <= 128, 'UNEXPECTED_LONG_ARRAY'
        for child in value:
            reject_individual_payload(child)


def validate_completed(task, code):
    for name in ('private_actual_failure_exit_v1.json', 'private_source_not_accepted_v1.json',
                 'private_native_gate_failure_v1.json', 'private_local_relay_failure_v1.json',
                 'private_local_observation_failure_v1.json'):
        if (task/name).exists():
            raise ValueError('TERMINAL_SOURCE_FAILURE_REQUIRES_REVIEW')
    names = ['private_frozen_fixed_native_workflow_v1.json', 'private_complete_input_identity_v1.json',
             'private_complete_source_fixed_workflow_v1.json', 'private_actual_exit_v1.json',
             'private_actual_source_start_v1.json']
    if not all((task/n).is_file() for n in names):
        return None
    plan, identity, source, exited, started = (json.loads((task/n).read_text()) for n in names)
    assert exited['status'] == 'SOURCE_EXIT0_FIXED_WORKFLOW_AND_RAM_RELEASE_COMPLETE'
    assert type(exited['actual_exit_code']) is int and exited['actual_exit_code'] == 0
    assert exited['frozen_contract_SHA256'] == digest(task/names[0])
    assert exited['complete_source_SHA256'] == digest(task/names[2])
    assert source['kind'] == 'COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED'
    assert source['started_UTC'] == started['started_UTC']
    assert type(source['source_R_exit_code']) is int and source['source_R_exit_code'] == 0
    for flag in ('source_payload_RAM_released', 'no_dynamic_followups_or_native_cache',
                 'no_native_individual_clinical_coordinate_state_or_graph_files',
                 'no_individual_values_or_identifiers_transferred'):
        assert source[flag] is True
    for source_key, identity_key, plan_key in (
        ('input_bytes', 'bytes', 'expected_file_bytes'), ('input_MD5', 'MD5', 'expected_file_MD5'),
        ('input_SHA256', 'SHA256', 'expected_file_SHA256')):
        assert source[source_key] == identity[identity_key] == plan[plan_key]
    for name, sha in plan['code_SHA256'].items():
        assert Path(name).name == name and name not in ('.', '..')
        assert digest(code/name) == sha, 'FROZEN_CODE_BYTES_CHANGED'
    joined = '\n'.join((code/n).read_text() for n in plan['R_components'])
    assert hashlib.sha256(joined.encode()).hexdigest() == source['R_code_SHA256'] == plan['R_code_SHA256']
    assert identity['range_assembly_proof'] == source['range_assembly_proof']
    proof = source['range_assembly_proof']
    assert proof['complete_ranges'] == math.ceil(plan['expected_file_bytes']/plan['range_chunk_bytes'])
    for key in ('exact_ranges_and_lengths', 'no_gaps_or_overlaps', 'whole_official_MD5_verified',
                'no_range_payload_files', 'no_alternate_endpoint_or_authentication_protocol_hash_retries'):
        assert proof[key] is True
    relay_digest = None
    if plan.get('input_transport') == 'bounded_local_opaque_ranges_over_existing_SSH':
        relay_path = task/'private_complete_local_opaque_stream_v1.json'
        assert relay_path.is_file(), 'LOCAL_COMPLETE_STREAM_PROOF_REQUIRED'
        local = json.loads(relay_path.read_text())
        assert local['kind'] == 'COMPLETE_ORDERED_PUBLIC_INPUT_PROOF'
        assert local['MD5'] == source['input_MD5'] and local['SHA256'] == source['input_SHA256']
        assert all(proof.get(k) == v for k,v in local['range_assembly_proof'].items())
        for key in ('ordered_stream_no_local_file_sized_buffer',
                    'opaque_compressed_relay_over_existing_SSH',
                    'whole_bytes_independently_hashed_on_source_before_R',
                    'range_header_proof_from_frozen_local_streamer'):
            assert proof[key] is True
        relay_digest = digest(relay_path)
    pipe = source['pipe_technical_proof']
    assert integer(pipe['expanded_stream_bytes']) and 0 < pipe['expanded_stream_bytes'] <= plan['maximum_expanded_bytes']
    assert pipe['no_native_or_expanded_payload_files'] is True
    times = [plan['frozen_UTC'], started['started_UTC'], identity['observed_UTC'],
             source['completed_UTC'], exited['observed_UTC']]
    parsed = [datetime.fromisoformat(t) for t in times]
    assert all(t.tzinfo is not None for t in parsed) and parsed == sorted(parsed)
    result = validate_summary(source['summary'], plan['analysis'])
    return dict(status='ACCEPTED_COMPLETED_SOURCE_WORKFLOW_BINDING',
                observed_UTC=datetime.now(timezone.utc).isoformat(), result=result,
                complete_source_SHA256=digest(task/names[2]), frozen_contract_SHA256=digest(task/names[0]),
                consumer_code_SHA256=digest(Path(__file__)),
                local_complete_stream_proof_SHA256=relay_digest,
                input_producer_proofs_bound_not_independent_local_payload_rehash=True)


def main():
    task = Path(os.environ['PANCANCER_EXTERNAL_FIXED_TASK'])
    result = validate_completed(task, Path(__file__).parent)
    if result is None:
        print('{"status":"PENDING_REQUIRED_COMPLETE_INPUT_SOURCE_AND_EXIT_NO_OUTPUT_CREATED"}')
        raise SystemExit(3)
    target = task/'private_accepted_fixed_spatial_workflow_v1.json'
    with target.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps(dict(status=result['status'], private_output_SHA256=digest(target))))


if __name__ == '__main__':
    main()
