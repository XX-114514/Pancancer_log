"""Counterexamples for premature, mismatched and split-axis qualification.

Only tiny artificial control JSON/code fixtures are files; no native data.
Fixtures are retained. No deletion or real input/worker restart is performed.
"""
from pathlib import Path
from copy import deepcopy
import hashlib
import json
import tempfile
import summarize_external_TDLN_qualification_v1 as consumer


def summary():
    markers = sorted(set(consumer.CONTEXT_DIMENSIONS).union(*map(set, consumer.MARKER_SETS.values())))
    rows = [dict(marker=m, matching_rows=1, numeric_RAM_assay_checked=True,
                 finite_cells=5, nonfinite_cells=0, nonzero_finite_cells=3, finite_values_vary=True) for m in markers]
    return dict(status='EXTERNAL_TDLN_STRUCTURE_MARKER_NUMERIC_QUALIFICATION_ONLY',
                no_states_thresholds_graphs_effects_or_clinical_values_computed=True,
                structural=dict(source_cells=5), schema=dict(assays=[dict(position=1, dimensions=[30, 5], numeric_in_RAM=True)]),
                predeclared_marker_numeric_profiles=[dict(assay_position=1, axis='SYNTHETIC_AXIS', markers=rows)])


def complete(value):
    return consumer.summarize_measurement(value)['marker_sets']['paired_regulatory_PD1_TIM3_with_proliferation_context']['status'] == 'COMPLETE_NUMERIC_MARKER_CANDIDATE_ON_ONE_AXIS'


def write_fixture(root, bad=None):
    root.mkdir(); code = root/'generic_code';code.mkdir()
    (code/'part.R').write_text('# Artificial generic R component; no native data\n')
    sha = consumer.digest(code/'part.R')
    plan = dict(expected_file_bytes=1024, expected_file_MD5='a'*32, range_chunk_bytes=1024,
                maximum_expanded_bytes=4096, code_SHA256={'part.R': sha}, R_components=['part.R'], R_code_SHA256=sha)
    proof = dict(complete_ranges=1, exact_ranges_and_lengths=True, no_gaps_or_overlaps=True,
                 whole_official_MD5_verified=True, no_range_payload_files=True,
                 no_alternate_endpoint_or_authentication_protocol_hash_retries=True)
    start = dict(started_UTC='2026-01-01T00:00:00+00:00')
    identity = dict(bytes=1024, MD5='a'*32, SHA256='b'*64, observed_UTC='2026-01-01T00:00:01+00:00', range_assembly_proof=proof)
    source = dict(kind='COMPLETE_EXTERNAL_TDLN_QUALIFICATION_RAM_RELEASED', started_UTC=start['started_UTC'],
                  completed_UTC='2026-01-01T00:00:02+00:00', source_R_exit_code=0, source_payload_RAM_released=True,
                  no_dynamic_followups_or_native_cache=True, no_native_individual_clinical_coordinate_state_or_graph_files=True,
                  no_individual_values_or_identifiers_transferred=True, input_bytes=1024, input_MD5='a'*32,
                  input_SHA256='b'*64, R_code_SHA256=sha, range_assembly_proof=proof, summary=summary(),
                  pipe_technical_proof=dict(expanded_stream_bytes=2048, no_native_or_expanded_payload_files=True))
    if bad == 'MD5':identity['MD5']='c'*32
    for name, value in [('private_frozen_native_qualification_v1.json', plan), ('private_actual_source_start_v1.json', start),
                        ('private_complete_input_identity_v1.json', identity), ('private_complete_source_qualification_v1.json', source)]:
        (root/name).write_text(json.dumps(value)+'\n')
    exited = dict(status='SOURCE_EXIT0_QUALIFICATION_AND_RAM_RELEASE_COMPLETE', actual_exit_code=0,
                  frozen_contract_SHA256=consumer.digest(root/'private_frozen_native_qualification_v1.json'),
                  complete_source_SHA256=consumer.digest(root/'private_complete_source_qualification_v1.json'))
    if bad == 'EXIT_BINDING':exited['complete_source_SHA256']='d'*64
    if bad == 'CODE':(code/'part.R').write_text('# Changed artificial component\n')
    (root/'private_actual_exit_v1.json').write_text(json.dumps(exited)+'\n')
    return code


def main():
    base = summary(); assert complete(base)
    for mode in ['missing', 'zero', 'ambiguous', 'nonfinite', 'unreviewed']:
        bad = deepcopy(base); rows=bad['predeclared_marker_numeric_profiles'][0]['markers']
        row = next(v for v in rows if v['marker']=='PD1')
        if mode == 'missing':rows.remove(row)
        if mode == 'zero':row.update(nonzero_finite_cells=0, finite_values_vary=False)
        if mode == 'ambiguous':row['matching_rows']=2
        if mode == 'nonfinite':row.update(finite_cells=4, nonfinite_cells=1)
        if mode == 'unreviewed':row['numeric_RAM_assay_checked']='true'
        assert not complete(bad), mode
        if mode == 'nonfinite':
            decision=consumer.summarize_measurement(bad)['marker_sets']['paired_regulatory_PD1_TIM3_with_proliferation_context']
            assert decision['status']=='POTENTIAL_FINITE_SUBSET_COMMON_CELLS_UNQUALIFIED'
            assert len(decision['potential_finite_subset_axes'])==1
    split = deepcopy(base); profile=split['predeclared_marker_numeric_profiles'][0]
    other=deepcopy(profile);other['axis']='DIFFERENT_SYNTHETIC_AXIS'
    for row in profile['markers']:
        if row['marker'] in ['PD1','TIM3']:row['matching_rows']=0
    for row in other['markers']:
        if row['marker'] not in ['PD1','TIM3']:row['matching_rows']=0
    split['predeclared_marker_numeric_profiles'].append(other)
    assert not complete(split), 'CROSS_AXIS_UNION_MUST_NOT_PASS'
    root=Path(tempfile.mkdtemp(prefix='pancancer_TDLN_consumer_synthetic_'))
    assert consumer.validate_completed(root, root) is None
    failed=root/'terminal_failure';failed.mkdir()
    (failed/'private_actual_failure_exit_v1.json').write_text('{"actual_exit_code":2}\n')
    try:consumer.validate_completed(failed,root)
    except ValueError:pass
    else:raise AssertionError('TERMINAL_FAILURE_MUST_NOT_BE_LABELLED_PENDING')
    code=write_fixture(root/'complete')
    assert consumer.validate_completed(root/'complete',code)['status']=='ACCEPTED_BOUND_SOURCE_QUALIFICATION_NO_NEW_EFFECTS'
    for label in ['MD5','EXIT_BINDING','CODE']:
        case=root/label;code=write_fixture(case,label)
        try:consumer.validate_completed(case,code)
        except AssertionError:pass
        else:raise AssertionError('MISMATCH_ACCEPTED_'+label)
    print('{"status":"PASS_MEASUREMENT_AND_INPUT_BINDING_COUNTEREXAMPLES_NO_NATIVE_DATA"}')


if __name__ == '__main__':
    main()
