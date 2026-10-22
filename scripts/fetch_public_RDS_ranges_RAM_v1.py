"""Bounded exact HTTP ranges, whole-object MD5, and process RAM only."""
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.request import Request, build_opener, ProxyHandler
from urllib.parse import urlparse
import hashlib
import io
import threading
import time


class RangeInputError(RuntimeError):
    def __init__(self, code_name):
        self.code_name = code_name
        super().__init__(code_name)


def fetch_ranges(url, *, expected_bytes, expected_MD5, chunk_bytes, workers,
                 request_timeout, maximum_range_seconds, progress=None,
                 opener_factory=None):
    assert urlparse(url).scheme == 'https' and urlparse(url).hostname == 'zenodo.org'
    assert 0 < expected_bytes <= 1024**3 and 1024 <= chunk_bytes <= 16*1024**2
    assert 1 <= workers <= 4 and 0 < request_timeout <= 60
    ranges = [(start, min(expected_bytes-1, start+chunk_bytes-1))
              for start in range(0, expected_bytes, chunk_bytes)]
    assert len(ranges) <= 1024
    factory = opener_factory or (lambda: build_opener(ProxyHandler({})))
    stop = threading.Event()

    def read_range(start, end):
        if stop.is_set():
            raise RangeInputError('CANCELLED')
        began = time.monotonic()
        with factory().open(Request(url, headers={'Range': f'bytes={start}-{end}'}),
                            timeout=request_timeout) as response:
            if response.status != 206 or urlparse(response.url).hostname != 'zenodo.org':
                raise RangeInputError('NOT_EXACT_OFFICIAL_206')
            if response.headers.get('Content-Range') != f'bytes {start}-{end}/{expected_bytes}':
                raise RangeInputError('CONTENT_RANGE_MISMATCH')
            length = end-start+1
            declared = response.headers.get('Content-Length')
            if declared is not None and int(declared) != length:
                raise RangeInputError('CONTENT_LENGTH_MISMATCH')
            buffer = io.BytesIO()
            while True:
                if stop.is_set():
                    raise RangeInputError('CANCELLED')
                if time.monotonic()-began > maximum_range_seconds:
                    raise RangeInputError('RANGE_ELAPSED_BUDGET')
                chunk = response.read1(min(65536, length-buffer.tell()+1))
                if not chunk:
                    break
                buffer.write(chunk)
                if buffer.tell() > length:
                    raise RangeInputError('RANGE_LENGTH_EXCESS')
            if buffer.tell() != length:
                raise RangeInputError('RANGE_LENGTH_SHORT')
            return start, end, buffer.getvalue()

    output = io.BytesIO()
    output.seek(expected_bytes-1)
    output.write(b'\0')
    completed, count = set(), 0
    pool = ThreadPoolExecutor(max_workers=workers)
    futures = []
    try:
        futures = [pool.submit(read_range, start, end) for start, end in ranges]
        for future in as_completed(futures):
            start, end, value = future.result()
            assert (start, end) not in completed and len(value) == end-start+1
            output.seek(start)
            output.write(value)
            count += len(value)
            completed.add((start, end))
            futures.remove(future)
            del value
            if progress:
                progress(count)
        assert completed == set(ranges) and count == expected_bytes
        view = output.getbuffer()
        md5, sha = hashlib.md5(view), hashlib.sha256(view)
        view.release()
        if md5.hexdigest() != expected_MD5:
            raise RangeInputError('WHOLE_OBJECT_MD5_MISMATCH')
        output.seek(expected_bytes)
        return output, md5, sha, dict(complete_ranges=len(ranges), workers=workers,
            exact_ranges_and_lengths=True, no_gaps_or_overlaps=True,
            whole_official_MD5_verified=True, no_range_payload_files=True,
            no_retries_or_alternate_endpoints=True)
    except BaseException:
        output.close()
        raise
    finally:
        stop.set()
        for future in futures:
            future.cancel()
        pool.shutdown(wait=True, cancel_futures=True)
