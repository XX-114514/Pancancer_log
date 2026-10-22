"""Bind a completed RAM qualification to its input, code and actual exit.

Only source aggregates are consumed. Numeric marker coverage is a measurement
candidate, not a cell state, independently reviewed patient, or clinical result.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import math
import os
import re


MARKER_SETS = {
    'paired_regulatory_PD1_TIM3_with_proliferation_context':
        ['CD3', 'CD4', 'CD8', 'FOXP3', 'CTLA4', 'PD1', 'TIM3', 'MKI67'],
    'paired_regulatory_CCR7_IL7R_with_proliferation_context':
        ['CD3', 'CD4', 'CD8', 'FOXP3', 'CTLA4', 'CCR7', 'IL7R', 'MKI67'],
    'regulatory_TBX21_marker_dimension': ['CD3', 'CD4', 'FOXP3', 'CTLA4', 'TBX21'],
}
CONTEXT_DIMENSIONS = ['CD86', 'PDL1', 'CD28', 'CD11c', 'CD68', 'CD163', 'LAMP3',
                      'GZMB', 'GZMK', 'TCF7', 'IFNG', 'CD45', 'CD19', 'CD20', 'PLA2G2D']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def positive_integer(value):
    return type(value) is int and value > 0


def marker_status(row, cells):
    if row is None or type(row.get('matching_rows')) is not int:
        return 'UNREVIEWED_MATCH_NOT_ASSUMED_ZERO'
    if row['matching_rows'] == 0:
        return 'NO_DECLARED_MATCH_ON_THIS_AXIS'
    if row['matching_rows'] != 1:
        return 'AMBIGUOUS_DECLARED_MATCH'
    if row.get('numeric_RAM_assay_checked') is not True:
        return 'NUMERIC_ASSAY_NOT_QUALIFIED'
    keys = ['finite_cells', 'nonfinite_cells', 'nonzero_finite_cells']
    if not all(type(row.get(k)) is int for k in keys):
        return 'NUMERIC_SUMMARY_UNREVIEWED'
    finite, nonfinite, nonzero = (row[k] for k in keys)
    if finite < 0 or nonfinite < 0 or finite+nonfinite != cells or not 0 <= nonzero <= finite:
        raise ValueError('INCONSISTENT_SOURCE_MARKER_COUNTS')
    if finite == 0:
        return 'NO_FINITE_MEASUREMENTS'
    if nonzero == 0:
        return 'ALL_FINITE_VALUES_ZERO'
    if row.get('finite_values_vary') is not True:
        return 'FINITE_VARIATION_NOT_QUALIFIED'
    if finite != cells or nonfinite != 0:
        return 'PARTIAL_FINITE_MEASUREMENT_COMMON_CELLS_UNQUALIFIED'
    return 'UNIQUE_FINITE_VARIABLE_NUMERIC_ROW_CANDIDATE'


def summarize_measurement(summary):
    assert summary['status'] == 'EXTERNAL_TDLN_STRUCTURE_MARKER_NUMERIC_QUALIFICATION_ONLY'
    assert summary['no_states_thresholds_graphs_effects_or_clinical_values_computed'] is True
    structure = summary['structural']; cells = structure['source_cells']
    assert positive_integer(cells)
    assays = {v['position']: v for v in summary['schema']['assays']}
    axes = []
    for profile in summary['predeclared_marker_numeric_profiles']:
        position = profile['assay_position']; assay = assays[position]
        dims = assay['dimensions']
        assert len(dims) == 2 and positive_integer(dims[0]) and dims[0] <= 256 and dims[1] == cells
        rows = profile['markers']
        assert len({v['marker'] for v in rows}) == len(rows), 'DUPLICATE_MARKER_SUMMARY_LABEL'
        by_name = {v['marker']: v for v in rows}
        statuses = {m: marker_status(by_name.get(m), cells)
                    for m in sorted(set(CONTEXT_DIMENSIONS).union(*map(set, MARKER_SETS.values())))}
        if assay.get('numeric_in_RAM') is not True:
            assert not any(v == 'UNIQUE_FINITE_VARIABLE_NUMERIC_ROW_CANDIDATE' for v in statuses.values())
        axes.append(dict(assay_position=position, axis=profile['axis'], markers=statuses))
    sets = {}
    for name, markers in MARKER_SETS.items():
        complete = [dict(assay_position=v['assay_position'], axis=v['axis']) for v in axes
                    if all(v['markers'][m] == 'UNIQUE_FINITE_VARIABLE_NUMERIC_ROW_CANDIDATE' for m in markers)]
        partial = [dict(assay_position=v['assay_position'], axis=v['axis']) for v in axes
                   if all(v['markers'][m] in ['UNIQUE_FINITE_VARIABLE_NUMERIC_ROW_CANDIDATE',
                                             'PARTIAL_FINITE_MEASUREMENT_COMMON_CELLS_UNQUALIFIED'] for m in markers)
                   and any(v['markers'][m] == 'PARTIAL_FINITE_MEASUREMENT_COMMON_CELLS_UNQUALIFIED' for m in markers)]
        sets[name] = dict(status='COMPLETE_NUMERIC_MARKER_CANDIDATE_ON_ONE_AXIS' if complete
                         else ('POTENTIAL_FINITE_SUBSET_COMMON_CELLS_UNQUALIFIED' if partial
                               else 'NO_COMPLETE_QUALIFIED_AXIS'), required_markers=markers,
                         complete_axis_candidates=complete,
                         potential_finite_subset_axes=partial,
                         partial_finite_measurement_does_not_discard_dataset=True,
                         no_cross_assay_or_axis_union=True,
                         duplicate_axes_are_not_independent_evidence=True)
    return dict(status='SOURCE_MEASUREMENT_CANDIDATES_ONLY_NOT_STATES_OR_CLINICAL_VALIDATION',
                marker_sets=sets, axis_qualification=axes,
                actual_cells_not_exported=True, lineage_or_function_accepted=False,
                coordinate_scale_or_pathology_compartment_accepted=False,
                patient_independence_or_clinical_endpoint_accepted=False,
                no_missing_marker_filled_as_negative=True)


def validate_completed(task, code):
    if (task/'private_actual_failure_exit_v1.json').exists() or (task/'private_source_not_accepted_v1.json').exists():
        raise ValueError('SOURCE_TERMINAL_FAILURE_REQUIRES_REVIEW')
    names = ['private_frozen_native_qualification_v1.json', 'private_complete_input_identity_v1.json',
             'private_complete_source_qualification_v1.json', 'private_actual_exit_v1.json',
             'private_actual_source_start_v1.json']
    if not all((task/n).is_file() for n in names):
        return None
    plan, identity, source, exited, started = (json.loads((task/n).read_text()) for n in names)
    assert exited['status'] == 'SOURCE_EXIT0_QUALIFICATION_AND_RAM_RELEASE_COMPLETE'
    assert type(exited['actual_exit_code']) is int and exited['actual_exit_code'] == 0
    assert exited['frozen_contract_SHA256'] == digest(task/names[0])
    assert exited['complete_source_SHA256'] == digest(task/names[2])
    assert source['kind'] == 'COMPLETE_EXTERNAL_TDLN_QUALIFICATION_RAM_RELEASED'
    assert source['started_UTC'] == started['started_UTC']
    assert type(source['source_R_exit_code']) is int and source['source_R_exit_code'] == 0
    assert source['source_payload_RAM_released'] is True
    assert source['no_dynamic_followups_or_native_cache'] is True
    assert source['no_native_individual_clinical_coordinate_state_or_graph_files'] is True
    assert source['no_individual_values_or_identifiers_transferred'] is True
    assert source['input_bytes'] == identity['bytes'] == plan['expected_file_bytes']
    assert source['input_MD5'] == identity['MD5'] == plan['expected_file_MD5']
    assert source['input_SHA256'] == identity['SHA256']
    assert re.fullmatch('[0-9a-f]{64}', source['input_SHA256'])
    for name, sha in plan['code_SHA256'].items():
        assert '/' not in name and name not in ['.', '..']
        assert digest(code/name) == sha, 'FROZEN_CODE_BYTES_CHANGED'
    R = '\n'.join((code/n).read_text() for n in plan['R_components'])
    assert hashlib.sha256(R.encode()).hexdigest() == source['R_code_SHA256'] == plan['R_code_SHA256']
    for proof in [identity['range_assembly_proof'], source['range_assembly_proof']]:
        assert proof['complete_ranges'] == math.ceil(plan['expected_file_bytes']/plan['range_chunk_bytes'])
        for k in ['exact_ranges_and_lengths', 'no_gaps_or_overlaps',
                  'whole_official_MD5_verified', 'no_range_payload_files',
                  'no_alternate_endpoint_or_authentication_protocol_hash_retries']:
            assert proof[k] is True
    assert identity['range_assembly_proof'] == source['range_assembly_proof']
    pipe = source['pipe_technical_proof']
    assert positive_integer(pipe['expanded_stream_bytes'])
    assert pipe['expanded_stream_bytes'] <= plan['maximum_expanded_bytes']
    assert pipe['no_native_or_expanded_payload_files'] is True
    times = [started['started_UTC'], identity['observed_UTC'], source['completed_UTC']]
    assert list(map(datetime.fromisoformat, times)) == sorted(map(datetime.fromisoformat, times))
    result = summarize_measurement(source['summary'])
    return dict(status='ACCEPTED_BOUND_SOURCE_QUALIFICATION_NO_NEW_EFFECTS',
                observed_UTC=datetime.now(timezone.utc).isoformat(), measurement=result,
                complete_source_SHA256=digest(task/names[2]), frozen_contract_SHA256=digest(task/names[0]),
                consumer_code_SHA256=digest(Path(__file__)),
                input_producer_proofs_bound_not_an_independent_native_payload_rehash=True)


def main():
    task = Path(os.environ['PANCANCER_EXTERNAL_TDLN_TASK'])
    result = validate_completed(task, Path(__file__).parent)
    if result is None:
        print('{"status":"PENDING_REQUIRED_COMPLETE_INPUT_SOURCE_AND_EXIT_NO_OUTPUT_CREATED"}')
        raise SystemExit(3)
    target = task/'private_accepted_measurement_qualification_v1.json'
    assert not target.exists()
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=result['status'], private_output_SHA256=digest(target))))


if __name__ == '__main__':
    main()
