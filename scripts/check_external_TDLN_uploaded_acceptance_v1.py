"""Artificial completed-upload receipts: full binding and scientific guards."""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
import tempfile
import accept_external_TDLN_uploaded_fixed_v1 as consumer
from check_external_TDLN_fixed_acceptance_v1 import artificial_summary, rejected
from accept_external_TDLN_fixed_spatial_v1 import STATE_MARKERS, ENDS


def check():
    root=Path(tempfile.mkdtemp(prefix='ARTIFICIAL_uploaded_receipt_checks_'))
    code=root/'code';code.mkdir()
    (code/'ARTIFICIAL_fixed.R').write_text('# ARTIFICIAL R IDENTITY ONLY\n')
    rsha=consumer.digest(code/'ARTIFICIAL_fixed.R')
    input_values=dict(input_bytes=64,input_MD5='a'*32,input_SHA256='b'*64)
    prior=dict(input_values,record_metadata_SHA256='c'*64)
    prior_path=root/'ARTIFICIAL_prior.json';prior_path.write_text(json.dumps(prior))
    analysis=dict(physical_panel=sorted(set(m for k,v in STATE_MARKERS.items() if k!=ENDS[2] for m in v)),
                  spatial=dict(minimum_cases=5,primary_grid=4,permutations=999))
    plan=dict(input_transport='preexisting_user_SFTP_file_read_only_RAM',
        expected_file_bytes=64,expected_file_MD5='a'*32,expected_file_SHA256='b'*64,
        code_SHA256={'ARTIFICIAL_fixed.R':rsha},R_components=['ARTIFICIAL_fixed.R'],R_code_SHA256=rsha,
        prior_qualification_evidence_path=str(prior_path),prior_qualification_evidence_SHA256=consumer.digest(prior_path),
        prior_record_metadata_SHA256='c'*64,maximum_expanded_bytes=10000,analysis=analysis,
        frozen_UTC='2026-01-01T00:00:00+00:00')
    proof=dict(status='COMPLETE_USER_UPLOAD_READ_ONLY_RAM_INPUT_VERIFIED',bytes_read=64,
        **{k:True for k in ('exact_length_before_and_after','owned_regular_file','same_descriptor_and_path_inode',
            'unchanged_size_mtime_ctime_and_mode','whole_official_MD5_verified','whole_frozen_SHA256_verified',
            'second_descriptor_pass_whole_hash_verified',
            'no_network_requests','no_input_write_rename_or_permission_change',
            'no_additional_native_or_expanded_payload_files','existing_user_file_retained')})
    started=dict(started_UTC='2026-01-01T00:01:00+00:00')
    identity=dict(bytes=64,MD5='a'*32,SHA256='b'*64,uploaded_file_read_proof=proof,
                  observed_UTC='2026-01-01T00:02:00+00:00')
    source=dict(input_values,kind='COMPLETE_EXTERNAL_FIXED_WORKFLOW_RAM_RELEASED',
        started_UTC=started['started_UTC'],completed_UTC='2026-01-01T00:03:00+00:00',source_R_exit_code=0,
        R_code_SHA256=rsha,record_metadata_SHA256='c'*64,summary=artificial_summary(),
        uploaded_file_read_proof=proof,pipe_technical_proof=dict(expanded_stream_bytes=128,no_native_or_expanded_payload_files=True),
        **{k:True for k in ('source_payload_RAM_released','no_dynamic_followups_or_additional_native_cache',
            'existing_user_upload_read_only_and_retained','no_native_individual_clinical_coordinate_state_or_graph_files',
            'no_individual_values_or_identifiers_transferred')})
    checks=[]
    def fixture(name, mutate=None):
        task=root/name;task.mkdir();p,i,s,t=deepcopy((plan,identity,source,started))
        supervisor=dict(actual_worker_exit_code=0,effective_exit_code=0,complete_source_frames=1,failure=None)
        if mutate:mutate(p,i,s,supervisor)
        def put(name,value):
            (task/name).write_text(json.dumps(value))
        put('private_frozen_fixed_native_workflow_v1.json',p)
        put('private_complete_input_identity_v1.json',i)
        put('private_complete_source_fixed_workflow_v1.json',s)
        put('private_actual_source_start_v1.json',t)
        supervisor['source_plan_SHA256']=consumer.digest(task/'private_frozen_fixed_native_workflow_v1.json')
        put('private_source_supervisor_actual_exit_v1.json',supervisor)
        put('private_actual_exit_v1.json',dict(status='SOURCE_EXIT0_FIXED_WORKFLOW_AND_RAM_RELEASE_COMPLETE',actual_exit_code=0,
            exit_code_scope='ACTUAL_SOURCE_WORKER_WAIT_RESULT_FROM_SUPERVISOR',
            frozen_contract_SHA256=consumer.digest(task/'private_frozen_fixed_native_workflow_v1.json'),
            complete_source_SHA256=consumer.digest(task/'private_complete_source_fixed_workflow_v1.json'),
            source_supervisor_exit_receipt_SHA256=consumer.digest(task/'private_source_supervisor_actual_exit_v1.json'),
            observed_UTC='2026-01-01T00:04:00+00:00'))
        return task
    task=fixture('ARTIFICIAL_POSITIVE_NEGATIVE_AND_UNMEASURED')
    result=consumer.validate_completed(task,code)
    assert result['result']['Python_aggregate_audit']['estimable_primary_slots']==10
    checks.append('COMPLETE_BINDING_WITH_POSITIVE_NEGATIVE_AND_STRUCTURAL_NA')
    cases={
        'BAD_INPUT_PROOF':lambda p,i,s,x:s['uploaded_file_read_proof'].update(whole_frozen_SHA256_verified=False),
        'WRONG_INPUT_SHA':lambda p,i,s,x:s.update(input_SHA256='d'*64),
        'NO_ACTUAL_WORKER_ZERO':lambda p,i,s,x:x.update(actual_worker_exit_code=7),
        'INVALID_BH':lambda p,i,s,x:s['summary']['spatial']['primary_rows'][0].update(BH=0),
        'INDIVIDUAL_CANARY':lambda p,i,s,x:s['summary'].update(patient_ids=['ARTIFICIAL_CANARY']),
        'WRONG_PRIOR_METADATA':lambda p,i,s,x:s.update(record_metadata_SHA256='d'*64),
        'WRONG_R':lambda p,i,s,x:s.update(R_code_SHA256='e'*64),
        'WRONG_TRANSPORT':lambda p,i,s,x:p.update(input_transport='unverified_file'),
    }
    for name,mutate in cases.items():
        bad=fixture('ARTIFICIAL_'+name,mutate)
        rejected(lambda:consumer.validate_completed(bad,code));checks.append(name+'_REJECTED')
    pending=root/'ARTIFICIAL_PENDING';pending.mkdir()
    assert consumer.validate_completed(pending,code) is None
    checks.append('PENDING_HAS_NO_ACCEPTANCE')
    print(json.dumps(dict(status='PASS_ARTIFICIAL_UPLOADED_CONSUMER',checks=checks,
                         no_real_native_or_individual_data_read=True)))


if __name__=='__main__':
    check()
