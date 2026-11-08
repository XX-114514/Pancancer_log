"""Fault-injected download tests using synthetic bytes, no patient data/network."""
from pathlib import Path
import hashlib
import io
import json
import re
import sys
import urllib.error
import zipfile
import numpy as np
from resumable_public_zip_v1 import DiskRange, extract_verified, IntegrityError, exclusive_lock, BusyError

root = Path(sys.argv[1])
root.mkdir(parents=True, exist_ok=False)
payload = np.random.default_rng(44).integers(0, 256, 12000, dtype=np.uint8).tobytes()
buffer = io.BytesIO()
with zipfile.ZipFile(buffer, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
    archive.writestr('native.h5ad', payload)
raw = buffer.getvalue()
with zipfile.ZipFile(io.BytesIO(raw)) as archive:
    info = archive.getinfo('native.h5ad')
expected = {'bytes': len(payload), 'sha256': hashlib.sha256(payload).hexdigest(), 'CRC32': info.CRC}


class Response:
    def __init__(self, start, end, action):
        self.status = 200 if action == 'ignore_range' else 206
        total = len(raw) + 1 if action == 'changed_total' else len(raw)
        self.headers = {'Content-Range': f'bytes {start}-{end}/{total}'}
        self.data = raw[start:end + 1]
        if action == 'short':
            self.data = self.data[:-1]

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return False

    def read(self, size):
        return self.data[:size]


class Opener:
    def __init__(self, actions=()):
        self.actions = list(actions)
        self.calls = []

    def open(self, request, timeout):
        start, end = map(int, re.fullmatch(r'bytes=(\d+)-(\d+)', request.headers['Range']).groups())
        self.calls.append((start, end))
        action = self.actions.pop(0) if self.actions else 'ok'
        if action in ('429', '403', '503'):
            raise urllib.error.HTTPError(request.full_url, int(action), 'injected', {'Retry-After': '1'}, None)
        if action == 'timeout':
            raise TimeoutError('injected interruption')
        return Response(start, end, action)


def reader(folder, opener, attempts=2):
    return DiskRange('https://example.test/archive.zip', root / folder, len(raw),
                     block_bytes=1024, attempts=attempts, opener=opener,
                     sleeper=lambda seconds: None, max_new_bytes=1000000)


# Truncation and 429 must retry the same exact block, not accept partial bytes.
op = Opener(['short', '429', 'ok'])
r = reader('retry', op, 3)
assert r.read(300) == raw[:300] and len(op.calls) == 3
assert op.calls[0] == op.calls[1] == op.calls[2]

# Forced interruption after one accepted block: new reader reuses the first block.
op = Opener(['ok', 'timeout', 'timeout'])
r = reader('resume', op)
try:
    r.read(2200)
except RuntimeError:
    pass
else:
    raise AssertionError('Interrupted transfer reported success')
op2 = Opener()
again = reader('resume', op2)
assert again.read(2200) == raw[:2200]
assert all(start >= 1024 for start, end in op2.calls)

# Corrupt cache is retained and fetched again, rather than trusted by length.
meta = json.loads((root / 'resume/block_00000000.json').read_text())
corrupt = root / 'resume' / meta['file']
corrupt.write_bytes(b'X' * 1024)
op3 = Opener()
again = reader('resume', op3)
assert again.read(100) == raw[:100] and corrupt.exists() and len(op3.calls) == 1

# Full member acceptance: ZIP CRC, frozen size and independent SHA; idempotent.
op = Opener()
r = reader('member_cache', op)
receipt = extract_verified(r, 'native.h5ad', expected, root / 'member')
assert (root / 'member' / receipt['members'][0]['local_name']).read_bytes() == payload
network_calls = len(op.calls)
assert extract_verified(r, 'native.h5ad', expected, root / 'member') == receipt
assert len(op.calls) == network_calls
try:
    extract_verified(reader('bad_sha_cache', Opener()), 'native.h5ad',
                     {**expected, 'sha256': '0' * 64}, root / 'bad_sha')
except IntegrityError:
    pass
else:
    raise AssertionError('Wrong full-member SHA accepted')
assert not (root / 'bad_sha/complete_PRIVATE.json').exists()
assert list((root / 'bad_sha').glob('attempt_*'))

# A full ZIP-member extraction interrupted in its compressed stream can restart
# from the disk cache; the failed uncompressed attempt remains as evidence.
class FailAtSecondBlock(Opener):
    def open(self, request, timeout):
        start = int(request.headers['Range'].split('=')[1].split('-')[0])
        if start == 1024:
            self.calls.append((start, start + 1023))
            raise TimeoutError('interrupted compressed stream')
        return super().open(request, timeout)


try:
    extract_verified(reader('zip_restart_cache', FailAtSecondBlock()),
                     'native.h5ad', expected, root / 'zip_restart_member')
except RuntimeError:
    pass
else:
    raise AssertionError('Interrupted member extraction accepted')
assert not (root / 'zip_restart_member/complete_PRIVATE.json').exists()
old_attempts = list((root / 'zip_restart_member').glob('attempt_*'))
assert old_attempts
resume_op = Opener()
extract_verified(reader('zip_restart_cache', resume_op), 'native.h5ad', expected,
                 root / 'zip_restart_member')
assert all(path.exists() for path in old_attempts)
assert not any(start == 0 for start, end in resume_op.calls)

# Range semantics and archive identity cannot silently change on resume.
for action in ['ignore_range', 'changed_total']:
    try:
        reader(action, Opener([action, action])).read(100)
    except RuntimeError:
        pass
    else:
        raise AssertionError('Bad range semantics accepted')
op = Opener(['403', 'ok'])
try:
    reader('forbidden', op).read(100)
except urllib.error.HTTPError as error:
    assert error.code == 403 and len(op.calls) == 1
else:
    raise AssertionError('403 bypassed')
with exclusive_lock(root / 'busy.lock'):
    try:
        with exclusive_lock(root / 'busy.lock'):
            pass
    except BusyError:
        pass
    else:
        raise AssertionError('Concurrent writer accepted')
try:
    DiskRange('https://example.test/other.zip', root / 'resume', len(raw), block_bytes=1024)
except IntegrityError:
    pass
else:
    raise AssertionError('Different archive reused old cache')
print('PASS_SYNTHETIC_ONLY: truncated/429 retries, crash-resume blocks, corrupt-cache refetch, CRC+size+SHA, idempotence, wrong-range/403/source-change/concurrent-writer rejection')
