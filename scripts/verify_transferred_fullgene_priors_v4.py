"""Verify a fresh offline recovery; preserve v1/v2/v3 failure history.

The target is transferred before running this script into a new v4 folder.
LR/TF are copied from verified v2 inputs. All official MD5/local SHA checks
precede appending the previously absent target to frozen biological inputs.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import shutil
import subprocess
import sys
import traceback


def now():
    return datetime.now(timezone.utc).isoformat()


def read(path):
    return json.loads(Path(path).read_text())


def digest(path, kind='sha256'):
    result = hashlib.new(kind)
    with Path(path).open('rb') as handle:
        for part in iter(lambda: handle.read(8388608), b''):
            result.update(part)
    return result.hexdigest()


def write(path, value):
    temporary = Path(str(path) + '.atomic_tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    os.replace(temporary, path)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    annotation = Path(os.environ['PANCANCER_ANNOTATION_ROOT'])
    output = root / 'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1'
    code = root / 'fullgene_receiver_target_code_v1'
    destination = output / 'verified_prior_inputs_v4'
    contract = read(code / 'fullgene_offline_resource_recovery_v4_spec.json')
    plan = read(output / 'predeclared_plan.json')
    assert {x['key']: (x['size'], x['checksum']) for x in plan['prior']['files']} == {
        x['file']: (x['bytes'], x['official_md5']) for x in contract['files']}
    old = read(root / 'family_LR_fullgene_assets_v3.execution.json')
    assert old['status'] == 'TERMINAL_FAILED_ARIA2_RESOURCE_FETCH'
    for key in ['pid', 'child_pid']:
        assert not Path('/proc/' + str(old[key])).exists()
    manifest = root / 'family_LR_fullgene_assets_v4.execution.json'
    assert not manifest.exists() and destination.is_dir()
    assert not (destination / 'download_receipts.json').exists()
    target = 'ligand_target_matrix_nsga2r_final.rds'
    assert (destination / target).is_file()
    assert not (output / 'verified_prior_inputs_v2' / target).exists()
    register = root / 'CM_program_screen_evidence_register_v1.json'
    checkpoint = root / 'cross_source_cell_identity_programme_localization_v1/current_execution_checkpoint.json'
    control = annotation / 'control/handoffs/20261011_GSE243013_preflight_v1.json'
    current = read(register)
    assert current['schema'] == 'native_LR_program_evidence_register_v16' and current['active_jobs'] == []
    assert digest(current['latest_verified_handoff']) == plan['parent_handoff_sha256']
    for label, path in [('register', register), ('root', control), ('checkpoint', checkpoint)]:
        snapshot = output / 'execution_state_snapshots' / ('assets_v4_before_' + label + '.json')
        assert not snapshot.exists()
        shutil.copy2(path, snapshot)
    execution = dict(status='RUNNING', started_UTC=now(), pid=os.getpid(), child_pid=None,
                     role='offline_prior_verification_not_scientific_model', script_sha256=digest(__file__),
                     contract_sha256=digest(code / 'fullgene_offline_resource_recovery_v4_spec.json'))

    def state(active, status):
        entry = read(register)
        entry['active_jobs'] = active
        write(register, entry)
        handoff = read(control)
        handoff.update(execution_status=status, active_jobs=active, last_verified_utc=now())
        write(control, handoff)
        write(checkpoint, dict(verified_utc=now(), latest_complete_handoff=entry['latest_verified_handoff'],
                               execution_status=status, active_jobs=active,
                               verified_execution={str(manifest): execution},
                               scientific_release=False, axis_validation_complete=False))

    write(manifest, execution)
    state([dict(phase='fullgene_assets_v4', pid=os.getpid(), role='offline_resource_verifier',
                status_manifest=str(manifest), observed_UTC=now())], 'RUNNING_FULLGENE_OFFLINE_RESOURCE_VERIFICATION')
    result = 1
    try:
        subprocess.run(['python3', str(annotation / 'control/refresh_state.py')], check=True)
        records = []
        for item in contract['files']:
            started = now()
            path = destination / item['file']
            if item['file'] != target:
                assert not path.exists()
                source = output / 'verified_prior_inputs_v2' / item['file']
                assert source.stat().st_size == item['bytes']
                assert digest(source, 'md5') == item['official_md5'].split(':')[1]
                assert digest(source) == item['sha256']
                shutil.copy2(source, path)
            assert path.stat().st_size == item['bytes']
            assert digest(path, 'md5') == item['official_md5'].split(':')[1]
            assert digest(path) == item['sha256']
            records.append(dict(item, observed_md5='md5:' + digest(path, 'md5'), started_UTC=started,
                                finished_UTC=now(), method='offline_transferred_verified_local_target' if item['file'] == target else 'copied_official_MD5_verified_v2'))
        destination_v2 = output / 'verified_prior_inputs_v2' / target
        with (destination / target).open('rb') as source, destination_v2.open('xb') as handle:
            shutil.copyfileobj(source, handle, length=8388608)
        assert digest(destination_v2) == next(x['sha256'] for x in records if x['file'] == target)
        write(destination / 'download_receipts.json', dict(status='ALL3_OFFICIAL_MD5_VERIFIED',
              method='offline_transfer_verification_no_public_redownload', records=records))
        result = 0
    except Exception as error:
        execution['error_type'] = type(error).__name__
        traceback.print_exc()
    finally:
        execution.update(status='TERMINAL_SUCCESS_RESOURCE_ONLY' if result == 0 else 'TERMINAL_FAILED_OFFLINE_RESOURCE_VERIFICATION',
                         finished_UTC=now(), exit_code=result)
        write(manifest, execution)
        state([], execution['status'])
        try:
            subprocess.run(['python3', str(annotation / 'control/refresh_state.py')], check=True)
        except Exception as error:
            execution.update(status='TERMINAL_STATE_REFRESH_FAILED', exit_code=1,
                             resource_exit_code=result, refresh_error_type=type(error).__name__)
            write(manifest, execution)
            state([], execution['status'])
            result = 1
    print(json.dumps(execution), flush=True)
    raise SystemExit(result)


if __name__ == '__main__':
    main()
