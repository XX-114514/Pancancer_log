"""Bounded four-connection prefetch of one public member; shared verified cache.

The sequential acquisition and all older transfers remain untouched. Cache locks
prevent competing writers; repeated BusyError is a finite wait, not a forced unlock.
"""
from pathlib import Path
import json
import os
import sys
import time
import uuid
import zipfile
from concurrent.futures import ThreadPoolExecutor
from resumable_public_zip_v1 import DiskRange, BusyError, sha256, atomic_json, exclusive_lock

r = Path(sys.argv[1]).resolve()
contract = json.loads((r / 'config/input_contract_PRIVATE.json').read_text())
item = next(v for v in contract['members'] if v['member'] == 'cell_boundaries.parquet')
out = r / 'boundary_acquisition_v1'
out.mkdir(exist_ok=True)
cache = r / 'acquisition_v1/archive_cache'
deadline = time.monotonic() + 4 * 3600


def reader():
    for attempt in range(180):
        try:
            return DiskRange(contract['url'], cache, contract['archive_bytes'],
                             block_bytes=1024**2, attempts=6, timeout=60,
                             max_new_bytes=512 * 1024**2,
                             deadline_seconds=max(1, deadline - time.monotonic()))
        except BusyError:
            if time.monotonic() >= deadline:
                raise
            time.sleep(1)
    raise TimeoutError('Cache identity lock wait expired')


def fetch(index):
    for attempt in range(180):
        if time.monotonic() >= deadline:
            raise TimeoutError('Finite acquisition deadline reached')
        try:
            source = reader()
            block = source._block(index)
            return len(block)
        except BusyError:
            time.sleep(1)
    raise TimeoutError('Cache block lock wait expired')


with exclusive_lock(out / 'writer.lock'):
    source = reader()
    with zipfile.ZipFile(source) as archive:
        info = archive.getinfo(item['member'])
        assert info.file_size == item['bytes'] and info.CRC == item['CRC32']
        start = info.header_offset // (1024**2)
        # Includes the complete ZIP local header, even with a maximal extra field.
        end = (info.header_offset + info.compress_size + 30 + 65535 + 65535) // (1024**2)
    with ThreadPoolExecutor(max_workers=4) as pool:
        accepted_bytes = sum(pool.map(fetch, range(start, end + 1)))
    attempt = out / (uuid.uuid4().hex + '.partial')
    with zipfile.ZipFile(reader()) as archive:
        with archive.open(item['member']) as source, attempt.open('xb') as target:
            for block in iter(lambda: source.read(1024**2), b''):
                target.write(block)
            target.flush(); os.fsync(target.fileno())
    assert attempt.stat().st_size == item['bytes']
    digest = sha256(attempt)
    destination = r / 'inputs' / item['local_name']
    try:
        os.link(attempt, destination)
    except FileExistsError:
        assert sha256(destination) == digest
    receipt = {**item, 'sha256': digest, 'member_CRC_validated': True,
               'independent_previous_SHA_match': False, 'cache_prefetched_bytes': accepted_bytes}
    atomic_json(r / 'inputs/boundary_complete_PRIVATE.json', receipt)
    atomic_json(out / 'complete_PRIVATE.json', receipt)
    print(json.dumps({'boundary_accepted': True, 'bytes': item['bytes'], 'CRC_validated': True,
                      'source_SHA_computed': True, 'max_connections': 4}), flush=True)
