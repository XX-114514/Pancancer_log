"""Persistent verified HTTP-range cache and restartable ZIP member extraction.

No source deletion. Cached blocks and failed extraction attempts survive a crash.
Whole-member CRC/size/SHA are required; partial files never become ready inputs.
"""
from pathlib import Path
from contextlib import contextmanager
from datetime import datetime, timezone
import fcntl
import hashlib
import http.client
import io
import json
import os
import random
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
import zipfile


class IntegrityError(RuntimeError):
    pass


class BusyError(RuntimeError):
    pass


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path, value):
    """Publish only complete metadata; preserve the former version separately."""
    path = Path(path)
    temp = path.with_name(path.name + '.' + uuid.uuid4().hex + '.writing')
    with temp.open('x') as handle:
        json.dump(value, handle, indent=2)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    if path.exists():
        os.link(path, path.with_name(path.name + '.' + uuid.uuid4().hex + '.history'))
    os.replace(temp, path)


@contextmanager
def exclusive_lock(path):
    with Path(path).open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise BusyError('Another writer owns this download') from error
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


class DiskRange(io.RawIOBase):
    """Virtual archive; every accepted block persists with an independent SHA."""
    def __init__(self, url, cache, archive_bytes, block_bytes=1024**2,
                 attempts=6, timeout=60, max_new_bytes=4 * 1024**3,
                 deadline_seconds=12 * 3600, opener=None, sleeper=time.sleep,
                 event=None):
        if urllib.parse.urlparse(url).scheme not in ('http', 'https'):
            raise ValueError('HTTP(S) public source required')
        if archive_bytes <= 0 or block_bytes <= 0 or attempts < 1:
            raise ValueError('Positive size, block size and attempts required')
        self.url, self.size = url, int(archive_bytes)
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.block, self.attempts, self.timeout = block_bytes, attempts, timeout
        self.network_bytes, self.reused_blocks = 0, 0
        self.max_new_bytes = max_new_bytes
        self.deadline = time.monotonic() + deadline_seconds
        self.opener = opener or urllib.request.build_opener(urllib.request.ProxyHandler({}))
        self.sleeper, self.event = sleeper, event or (lambda **data: None)
        self.pos = 0
        self.memory = {}
        self.identity = {'url_sha256': hashlib.sha256(url.encode()).hexdigest(),
                         'archive_bytes': self.size, 'block_bytes': self.block}
        identity_path = self.cache / 'identity.json'
        with exclusive_lock(self.cache / 'identity.lock'):
            if identity_path.exists():
                if json.loads(identity_path.read_text()) != self.identity:
                    raise IntegrityError('Cache belongs to a different archive or block layout')
            else:
                atomic_json(identity_path, self.identity)

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        value = offset if whence == 0 else self.pos + offset if whence == 1 else self.size + offset
        if not 0 <= value <= self.size:
            raise ValueError('Seek outside known archive')
        self.pos = value
        return value

    def _download(self, start, size):
        last = None
        for attempt in range(self.attempts):
            if time.monotonic() >= self.deadline:
                raise TimeoutError('Finite download deadline reached')
            if self.network_bytes + size > self.max_new_bytes:
                raise RuntimeError('Finite attempted range byte budget reached')
            self.network_bytes += size  # Count failed attempts conservatively too.
            wait = min(30., 2. ** attempt + random.random())
            try:
                request = urllib.request.Request(self.url, headers={
                    'Range': f'bytes={start}-{start + size - 1}',
                    'Accept-Encoding': 'identity', 'User-Agent': 'PanCancer-public-research/2.0'})
                with self.opener.open(request, timeout=self.timeout) as response:
                    content_range = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)',
                                                response.headers.get('Content-Range', ''))
                    if response.status != 206 or not content_range:
                        raise IntegrityError('Server did not honor exact HTTP range')
                    if tuple(map(int, content_range.groups())) != (start, start + size - 1, self.size):
                        raise IntegrityError('Range boundaries or archive total changed')
                    if response.headers.get('Content-Encoding', 'identity') != 'identity':
                        raise IntegrityError('Unexpected HTTP content encoding')
                    body = response.read(size + 1)
                    if len(body) != size:
                        raise IntegrityError('Truncated or oversized range body')
                return body
            except urllib.error.HTTPError as error:
                if error.code in (401, 403, 404, 410):
                    raise  # Access/permanent failures need review, not endless retries.
                if error.code not in (408, 425, 429, 500, 502, 503, 504):
                    raise
                retry_after = error.headers.get('Retry-After', '') if error.headers else ''
                if retry_after.isdigit():
                    wait = min(30., max(wait, float(retry_after)))
                last = error
            except (OSError, EOFError, http.client.HTTPException, IntegrityError) as error:
                last = error
            self.event(status='RANGE_RETRY', offset=start, requested_bytes=size,
                       attempt=attempt + 1, exception_type=type(last).__name__)
            if attempt + 1 < self.attempts:
                remaining = self.deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError('Finite download deadline reached') from last
                self.sleeper(min(wait, remaining))
        raise RuntimeError('Range retry budget exhausted') from last

    def _block(self, index):
        if index in self.memory:
            return self.memory[index]
        start = index * self.block
        size = min(self.block, self.size - start)
        metadata = self.cache / f'block_{index:08d}.json'
        with exclusive_lock(self.cache / f'block_{index:08d}.lock'):
            if metadata.exists():
                record = json.loads(metadata.read_text())
                filename = Path(record['file'])
                if filename.name != str(filename):
                    raise IntegrityError('Invalid cached block filename')
                path = self.cache / filename
                if path.exists() and path.stat().st_size == size and sha256(path) == record['sha256']:
                    self.reused_blocks += 1
                    data = path.read_bytes()
                    if len(self.memory) >= 4:
                        self.memory.pop(next(iter(self.memory)))
                    self.memory[index] = data
                    return data
                self.event(status='CORRUPT_CACHE_BLOCK_RETAINED_REFETCH', offset=start)
            data = self._download(start, size)
            path = self.cache / f'block_{index:08d}_{uuid.uuid4().hex}.bin'
            with path.open('xb') as handle:
                handle.write(data)
                handle.flush()
                os.fsync(handle.fileno())
            atomic_json(metadata, {'file': path.name, 'bytes': size,
                                   'sha256': hashlib.sha256(data).hexdigest()})
            if len(self.memory) >= 4:
                self.memory.pop(next(iter(self.memory)))
            self.memory[index] = data
            return data

    def read(self, size=-1):
        size = self.size - self.pos if size < 0 else min(size, self.size - self.pos)
        result = []
        while size > 0:
            data = self._block(self.pos // self.block)
            offset = self.pos % self.block
            take = min(size, len(data) - offset)
            if take <= 0:
                raise IntegrityError('Zero-progress archive read')
            result.append(data[offset:offset + take])
            self.pos += take
            size -= take
        return b''.join(result)


def extract_verified(reader, member, expected, out):
    """Resume compressed blocks; restart inflation safely into a fresh attempt.

    expected requires independently known uncompressed bytes, CRC32 and SHA256.
    A prior accepted file is fully rechecked before reuse, without HTTP access.
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with exclusive_lock(out / 'writer.lock'):
        complete = out / 'complete_PRIVATE.json'
        if complete.exists():
            receipt = json.loads(complete.read_text())
            item = receipt['members'][0]
            path = out / item['local_name']
            if (item['member'] == member and path.stat().st_size == expected['bytes']
                    and sha256(path) == expected['sha256']
                    and item['CRC32'] == expected['CRC32']):
                return receipt
            raise IntegrityError('Accepted output changed; preserve and investigate')
        attempt = out / ('attempt_' + uuid.uuid4().hex + '_PRIVATE.h5ad')
        digest, size = hashlib.sha256(), 0
        try:
            with zipfile.ZipFile(reader) as archive:
                info = archive.getinfo(member)
                if info.file_size != expected['bytes'] or info.CRC != expected['CRC32']:
                    raise IntegrityError('ZIP member identity differs from frozen inventory')
                with archive.open(info) as source, attempt.open('xb') as target:
                    for data in iter(lambda: source.read(1024**2), b''):
                        target.write(data)
                        digest.update(data)
                        size += len(data)
                    target.flush()
                    os.fsync(target.fileno())
            if size != expected['bytes'] or digest.hexdigest() != expected['sha256']:
                raise IntegrityError('Extracted member size/SHA differs from independent reference')
            item = {'member': member, 'local_name': attempt.name, 'bytes': size,
                    'sha256': digest.hexdigest(), 'CRC32': expected['CRC32'],
                    'member_CRC_validated': True}
            receipt = {'UTC': datetime.now(timezone.utc).isoformat(), 'members': [item],
                       'compressed_block_resume': True, 'whole_archive_MD5_verified': False,
                       'network_attempted_range_bytes_this_reader': reader.network_bytes,
                       'cache_blocks_reused_this_reader': reader.reused_blocks}
            atomic_json(complete, receipt)
            return receipt
        except Exception as error:
            atomic_json(out / ('failure_' + uuid.uuid4().hex + '_PRIVATE.json'), {
                'UTC': datetime.now(timezone.utc).isoformat(),
                'exception_type': type(error).__name__, 'partial_file_retained': attempt.name,
                'next_run': 'repeat same command; accepted compressed blocks resume'})
            raise
