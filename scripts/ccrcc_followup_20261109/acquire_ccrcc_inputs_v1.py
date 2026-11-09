"""Minimal official ccRCC source acquisition on server; preserve old partials."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import sys
import time
import traceback
import uuid
import zipfile
from resumable_public_zip_v1 import DiskRange, atomic_json, exclusive_lock, sha256

root = Path(sys.argv[1]).resolve()
contract = json.loads((root / 'config/input_contract_PRIVATE.json').read_text())
out = root / 'acquisition_v1'
out.mkdir(exist_ok=True)
inputs = root / 'inputs'
inputs.mkdir(exist_ok=True)
deadline = time.monotonic() + 12 * 3600


def event(**data):
    with (out / 'events_PRIVATE.jsonl').open('a') as f:
        f.write(json.dumps({'UTC': datetime.now(timezone.utc).isoformat(), **data}) + '\n')


with exclusive_lock(out / 'writer.lock'):
    accepted = {}
    for round_number in range(1, 4):
        for item in contract['members']:
            name = item['member']
            if name in accepted or time.monotonic() >= deadline:
                continue
            folder = out / item['alias']
            folder.mkdir(exist_ok=True)
            receipt_path = folder / 'complete_PRIVATE.json'
            try:
                if receipt_path.exists():
                    receipt = json.loads(receipt_path.read_text())
                    path = folder / receipt['local_name']
                    assert path.stat().st_size == item['bytes']
                    assert sha256(path) == receipt['sha256']
                    if item.get('sha256'):
                        assert receipt['sha256'] == item['sha256']
                    accepted[name] = receipt
                else:
                    event(stage=item['alias'], status='START_OR_RESUME_OFFICIAL_MEMBER', round=round_number)
                    reader = DiskRange(contract['url'], out / 'archive_cache', contract['archive_bytes'],
                                       block_bytes=1024**2, attempts=6, timeout=60,
                                       max_new_bytes=2 * 1024**3,
                                       deadline_seconds=max(1, deadline - time.monotonic()),
                                       event=lambda **data: event(stage=item['alias'], **data))
                    attempt = folder / ('attempt_' + uuid.uuid4().hex + '.partial')
                    digest, size = hashlib.sha256(), 0
                    with zipfile.ZipFile(reader) as archive:
                        info = archive.getinfo(name)
                        assert info.file_size == item['bytes'] and info.CRC == item['CRC32']
                        with archive.open(info) as source, attempt.open('xb') as target:
                            for block in iter(lambda: source.read(1024**2), b''):
                                target.write(block)
                                digest.update(block)
                                size += len(block)
                            target.flush()
                            os.fsync(target.fileno())
                    assert size == item['bytes']
                    if item.get('sha256'):
                        assert digest.hexdigest() == item['sha256']
                    receipt = {'member': name, 'local_name': attempt.name, 'bytes': size,
                               'sha256': digest.hexdigest(), 'CRC32': info.CRC,
                               'member_CRC_validated': True,
                               'independent_previous_SHA_match': bool(item.get('sha256'))}
                    atomic_json(receipt_path, receipt)
                    accepted[name] = receipt
                target = inputs / item['local_name']
                if not target.exists():
                    os.link(folder / accepted[name]['local_name'], target)
                assert target.stat().st_size == item['bytes']
                assert sha256(target) == accepted[name]['sha256']
                event(stage=item['alias'], status='CRC_SIZE_SHA_ACCEPTED', bytes=item['bytes'])
                core = contract['members'][:6]
                if all(v['member'] in accepted for v in core) and not (inputs / 'complete_PRIVATE.json').exists():
                    records = [{**accepted[v['member']], 'local_name': v['local_name']} for v in core]
                    atomic_json(inputs / 'complete_PRIVATE.json', {
                        'members': records, 'whole_archive_MD5_verified': False,
                        'all6_independent_previous_SHA_match': True})
                    event(stage='core6', status='ALL6_INPUTS_VERIFIED_READY_FOR_BASELINE_REPLAY')
                if name == 'cell_boundaries.parquet':
                    atomic_json(inputs / 'boundary_complete_PRIVATE.json',
                                {**accepted[name], 'local_name': item['local_name']})
            except Exception as error:
                (folder / ('failure_' + uuid.uuid4().hex + '_PRIVATE.txt')).write_text(traceback.format_exc())
                event(stage=item['alias'], status='FAILED_THIS_ROUND_PARTIAL_AND_CACHE_RETAINED',
                      exception_type=type(error).__name__)
        if len(accepted) == len(contract['members']) or time.monotonic() >= deadline:
            break
    atomic_json(out / 'terminal_PRIVATE.json', {
        'UTC': datetime.now(timezone.utc).isoformat(), 'requested_members': len(contract['members']),
        'accepted_members': len(accepted),
        'status': 'ALL7_ACCEPTED' if len(accepted) == len(contract['members']) else 'INCOMPLETE_RESUMABLE',
        'old_source_and_upload_preserved': True, 'local_deletions': 0})
