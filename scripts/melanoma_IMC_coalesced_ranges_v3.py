"""Same strict HTTPS transport, with bounded aligned RAM read caching.

Tiny ZIP footer/header reads are served from 64KiB physical ranges. All wire
requests still use the original status/identity/TLS/deadline/byte guards.
No compressed bytes are written to disk. clear_cache forces a wire replay.
"""
from collections import OrderedDict
from pathlib import Path
import hashlib
import importlib.util

path = Path(__file__).with_name('melanoma_IMC_verified_dns_transport_v2.py')
assert hashlib.sha256(path.read_bytes()).hexdigest() == '12ffe4ac9c1c37efefbe42e8153928473da3cee8997a2086118d46ba81d7c25c'
spec = importlib.util.spec_from_file_location('guarded_coalesced_base', path)
_base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_base)


class RangeFile(_base.RangeFile):
    def __init__(self, source, url, size, checkpoint=None):
        super().__init__(source, url, size, checkpoint)
        self.block_bytes = 64*1024
        self.cache_blocks = 8
        self.cache = OrderedDict()

    def clear_cache(self):
        self.cache.clear()

    def _fetch_block(self, start):
        count = min(self.block_bytes, self.size-start)
        position = self.position
        self.position = start
        try:
            # Includes overfetch/retries in the base wire budget. The logical
            # file pointer is restored even on a rejected or timed-out read.
            return _base.RangeFile._read_chunk(self, count)
        finally:
            self.position = position

    def read(self, n=-1):
        if n < 0:
            n = self.size-self.position
        n = min(n, self.size-self.position)
        assert 0 <= n <= 16*1024**2
        parts = []
        while n:
            # One complete tail block contains EOCD, ZIP64 probe and directory.
            start = max(0, self.size-self.block_bytes) if self.position >= self.size-self.block_bytes else self.position//self.block_bytes*self.block_bytes
            if start not in self.cache:
                block = self._fetch_block(start)
                assert len(block) == min(self.block_bytes, self.size-start)
                self.cache[start] = block
                if len(self.cache) > self.cache_blocks:
                    self.cache.popitem(last=False)
            self.cache.move_to_end(start)
            block = self.cache[start]
            offset = self.position-start
            count = min(n, len(block)-offset)
            assert count > 0
            parts.append(block[offset:offset+count])
            self.position += count
            n -= count
        return b''.join(parts)
