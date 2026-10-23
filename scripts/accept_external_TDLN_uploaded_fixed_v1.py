"""Bind a completed read-only user-upload workflow to the original fixed R.

Existing numeric, missingness and individual-payload checks are reused unchanged.
The user-provided input file is retained; no extra native cache is created.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
from accept_external_TDLN_fixed_spatial_v1 import digest, integer, validate_summary


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
    for flag in ('source_payload_RAM_released', 'no_dynamic_followups_or_additional_native_cache', 'existing_user_upload_read_only_and_retained',
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
    assert plan['input_transport'] == 'preexisting_user_SFTP_file_read_only_RAM'
    assert identity['uploaded_file_read_proof'] == source['uploaded_file_read_proof']
    proof = source['uploaded_file_read_proof']
    assert proof['status'] == 'COMPLETE_USER_UPLOAD_READ_ONLY_RAM_INPUT_VERIFIED'
    assert proof['bytes_read'] == plan['expected_file_bytes']
    for flag in ('exact_length_before_and_after', 'owned_regular_file',
                 'same_descriptor_and_path_inode', 'unchanged_size_mtime_ctime_and_mode',
                 'whole_official_MD5_verified', 'whole_frozen_SHA256_verified',
                 'second_descriptor_pass_whole_hash_verified',
                 'no_network_requests', 'no_input_write_rename_or_permission_change',
                 'no_additional_native_or_expanded_payload_files', 'existing_user_file_retained'):
        assert proof[flag] is True
    prior_path = Path(plan['prior_qualification_evidence_path'])
    assert digest(prior_path) == plan['prior_qualification_evidence_SHA256']
    prior = json.loads(prior_path.read_text())
    assert prior['record_metadata_SHA256'] == source['record_metadata_SHA256'] == plan['prior_record_metadata_SHA256']
    for k in ('input_MD5', 'input_SHA256', 'input_bytes'):
        assert prior[k] == source[k]
    supervisor_path = task/'private_source_supervisor_actual_exit_v1.json'
    assert digest(supervisor_path) == exited['source_supervisor_exit_receipt_SHA256']
    supervisor = json.loads(supervisor_path.read_text())
    assert supervisor['actual_worker_exit_code'] == supervisor['effective_exit_code'] == 0
    assert supervisor['complete_source_frames'] == 1 and supervisor['failure'] is None
    assert supervisor['source_plan_SHA256'] == digest(task/names[0])
    assert exited['exit_code_scope'] == 'ACTUAL_SOURCE_WORKER_WAIT_RESULT_FROM_SUPERVISOR'
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
                source_actual_worker_exit_receipt_SHA256=digest(supervisor_path),
                input_producer_proofs_bound_not_independent_local_payload_rehash=True)


def main():
    task = Path(os.environ['PANCANCER_EXTERNAL_FIXED_TASK'])
    result = validate_completed(task, Path(__file__).parent)
    if result is None:
        print('{"status":"PENDING_COMPLETE_UPLOAD_INPUT_SOURCE_AND_ACTUAL_EXIT_NO_OUTPUT_CREATED"}')
        raise SystemExit(3)
    path = task/'private_accepted_uploaded_fixed_spatial_workflow_v1.json'
    with path.open('x') as stream:
        json.dump(result,stream,indent=2)
        stream.write('\n')
    print(json.dumps(dict(status=result['status'],private_output_SHA256=digest(path))))


if __name__=='__main__':
    main()
