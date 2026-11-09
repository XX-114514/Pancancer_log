"""Three-connection core-member cache prefetch; no input publication or deletion."""
from pathlib import Path
import json
import sys
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from resumable_public_zip_v1 import DiskRange, BusyError, exclusive_lock, atomic_json

r = Path(sys.argv[1]).resolve()
c = json.loads((r / 'config/input_contract_PRIVATE.json').read_text())
out = r / 'core_prefetch_v1'
out.mkdir(exist_ok=True)
deadline = time.monotonic() + 4 * 3600


def reader():
    return DiskRange(c['url'], r / 'acquisition_v1/archive_cache', c['archive_bytes'],
                     block_bytes=1024**2, attempts=6, timeout=60,
                     max_new_bytes=256 * 1024**2,
                     deadline_seconds=max(1, deadline - time.monotonic()))


def fetch(i):
    for attempt in range(180):
        if time.monotonic() >= deadline:
            raise TimeoutError('Core cache prefetch deadline reached')
        try:
            return len(reader()._block(i))
        except BusyError:
            time.sleep(1)
    raise TimeoutError('Held cache block preserved; finite wait expired')


with exclusive_lock(out / 'writer.lock'):
    blocks = set()
    with zipfile.ZipFile(reader()) as z:
        for item in c['members'][:6]:
            info = z.getinfo(item['member'])
            assert info.CRC == item['CRC32'] and info.file_size == item['bytes']
            a = info.header_offset // 1024**2
            b = (info.header_offset + info.compress_size + 30 + 2 * 65535) // 1024**2
            blocks.update(range(a, b + 1))
    with ThreadPoolExecutor(max_workers=3) as pool:
        total = sum(pool.map(fetch, sorted(blocks)))
    atomic_json(out / 'complete.json', {'cached_blocks': len(blocks), 'cached_bytes': total,
                                        'not_member_acceptance': True, 'max_connections': 3})
