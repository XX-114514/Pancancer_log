"""Restartable finite recovery for failed public native inputs only.

Run the same command again after interruption. Previously accepted inputs are
rehashed and reused. Source/archive paths and donor mappings remain private.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import os
import time
import traceback
import urllib.error
import uuid
import anndata as ad
from resumable_public_zip_v1 import DiskRange, extract_verified, sha256, atomic_json, exclusive_lock, IntegrityError
from native_coordinate_contract_v1 import native_micron_coordinates


def recover(root, hours=18, rounds=3):
    root = Path(root).resolve()
    config = json.loads((root / 'config/cohort_PRIVATE.json').read_text())
    contract = json.loads((root / 'config/recovery_download_PRIVATE.json').read_text())
    out = root / 'download_recovery_v1'
    out.mkdir(exist_ok=True)
    deadline = time.monotonic() + hours * 3600
    events = out / 'events_PRIVATE.jsonl'

    def event(**data):
        with events.open('a') as handle:
            handle.write(json.dumps({'UTC': datetime.now(timezone.utc).isoformat(), **data}) + '\n')

    def coordinate_check(path):
        data = ad.read_h5ad(path, backed='r')
        try:
            return native_micron_coordinates(data)[1]
        finally:
            data.file.close()

    def accepted(case):
        folder = root / 'inputs_direct_v1' / case['unit']
        complete, ready = folder / 'complete_PRIVATE.json', folder / 'ready_PRIVATE.json'
        if not complete.exists() or not ready.exists():
            return False
        item = json.loads(complete.read_text())['members'][0]
        expected = case['receipt']['members'][0]
        path = folder / item['local_name']
        if (item['member'] != case['expected_member'] or path.stat().st_size != expected['bytes']
                or sha256(path) != expected['sha256']):
            raise IntegrityError('Previously accepted source no longer matches reference')
        coordinate_check(path)
        return True

    with exclusive_lock(out / 'recovery_writer.lock'):
        mapping = root / 'inputs_direct_v1/mapping'
        mappings = json.loads((mapping / 'complete_PRIVATE.json').read_text())
        expected_mapping = {v['member']: v for v in config['mapping_receipt']['members']}
        if len(mappings['members']) != 18 or len(expected_mapping) != 18:
            raise IntegrityError('Expected all18 accepted mapping files before recovery')
        for item in mappings['members']:
            path = mapping / item['local_name']
            expected = expected_mapping[item['member']]
            if path.stat().st_size != expected['bytes'] or sha256(path) != expected['sha256']:
                raise IntegrityError('Previously accepted mapping changed')
        event(stage='mapping', status='REVALIDATED_ALL18_REUSED_NO_REDOWNLOAD')
        pending = []
        for case in config['cases']:
            if accepted(case):
                event(stage=case['unit'], status='REVALIDATED_NATIVE_REUSED_NO_REDOWNLOAD')
            else:
                pending.append(case)
        event(stage='all6', status='RECOVERY_STARTED', failed_or_incomplete=len(pending),
              block_bytes=contract['block_bytes'], rounds=rounds, finite_hours=hours)
        blocked = set()
        for round_number in range(1, rounds + 1):
            for case in pending:
                unit = case['unit']
                if unit in blocked or time.monotonic() >= deadline:
                    continue
                try:
                    if accepted(case):
                        continue
                    event(stage=unit, status='RESUME_OR_START_FAILED_NATIVE', round=round_number)
                    expected = contract['members'][unit]
                    reader = DiskRange(contract['url'], out / 'archive_block_cache',
                                       contract['archive_bytes'], block_bytes=contract['block_bytes'],
                                       attempts=6, timeout=60, max_new_bytes=8 * 1024**3,
                                       deadline_seconds=max(1, deadline - time.monotonic()),
                                       event=lambda **data: event(stage=unit, **data))
                    receipt = extract_verified(reader, case['expected_member'], expected, out / unit)
                    item = receipt['members'][0]
                    path = out / unit / item['local_name']
                    coordinates = coordinate_check(path)
                    folder = root / 'inputs_direct_v1' / unit
                    target = folder / 'native_PRIVATE.h5ad'
                    if target.exists():
                        if target.stat().st_size != expected['bytes'] or sha256(target) != expected['sha256']:
                            raise IntegrityError('Existing target changed; retain for investigation')
                    else:
                        os.link(path, target)
                    complete = folder / 'complete_PRIVATE.json'
                    if not complete.exists():
                        atomic_json(complete, case['receipt'])
                    else:
                        old = json.loads(complete.read_text())['members'][0]
                        if old['member'] != item['member'] or old['sha256'] != expected['sha256']:
                            raise IntegrityError('Existing source receipt mismatch')
                    atomic_json(folder / 'ready_PRIVATE.json', {
                        'UTC': datetime.now(timezone.utc).isoformat(),
                        'server_CRC_size_and_independent_SHA_match': True,
                        'coordinate_contract': coordinates, 'download_recovery_version': 1})
                    event(stage=unit, status='RECOVERED_CRC_SIZE_SHA_COORDINATES_ACCEPTED_READY',
                          network_attempted_bytes=reader.network_bytes,
                          cached_block_reads=reader.reused_blocks)
                except urllib.error.HTTPError as error:
                    if error.code in (401, 403, 404, 410):
                        blocked.add(unit)
                    event(stage=unit, status='HTTP_FAILURE_RETAINED', code=error.code,
                          needs_source_or_access_review=unit in blocked, round=round_number)
                except Exception as error:
                    path = out / (unit + '_' + uuid.uuid4().hex + '_failure_PRIVATE.txt')
                    path.write_text(traceback.format_exc())
                    event(stage=unit, status='FAILED_THIS_ROUND_CACHE_RETAINED',
                          exception_type=type(error).__name__, round=round_number)
            pending = [case for case in config['cases'] if not accepted(case)]
            if not pending or time.monotonic() >= deadline:
                break
        terminal = {'UTC': datetime.now(timezone.utc).isoformat(),
                    'requested_native_sources': len(config['cases']),
                    'accepted_native_sources': len(config['cases']) - len(pending),
                    'still_failed_or_incomplete': [case['unit'] for case in pending],
                    'mapping_files_revalidated_and_reused': 18,
                    'status': 'ALL_ACCEPTED' if not pending else 'INCOMPLETE_RESUMABLE_REPEAT_SAME_COMMAND',
                    'old_failed_files_retained': True, 'local_deletions': 0}
        atomic_json(out / 'latest_terminal_PRIVATE.json', terminal)
        event(stage='all6', **terminal)
        return terminal


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run_root', type=Path)
    parser.add_argument('--hours', type=float, default=18)
    parser.add_argument('--rounds', type=int, default=3)
    args = parser.parse_args()
    if not 0 < args.hours <= 48 or not 1 <= args.rounds <= 5:
        parser.error('Finite bounds required: hours in(0,48], rounds in[1,5]')
    result = recover(args.run_root, args.hours, args.rounds)
    print(json.dumps(result), flush=True)
    raise SystemExit(0 if result['status'] == 'ALL_ACCEPTED' else 2)
