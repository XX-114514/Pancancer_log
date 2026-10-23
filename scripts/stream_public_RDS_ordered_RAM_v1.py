"""Bounded ordered public HTTPS ranges to a binary sink; no file-sized buffer.

Continuation rules follow fetch_public_RDS_ranges_resume_RAM_v1. This is a new
streaming protocol, not a replacement or monkeypatch of that frozen fetcher.
The receiver must independently verify whole bytes/MD5/SHA before parsing.
"""
from concurrent.futures import ThreadPoolExecutor
from http.client import IncompleteRead
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from urllib.request import Request, build_opener, ProxyHandler
import hashlib
import io
import ssl
import threading
import time


class StreamInputError(RuntimeError):
    def __init__(self, code_name):
        super().__init__(code_name)
        self.code_name = code_name


def stream_ranges(url, sink, *, expected_bytes, expected_MD5, expected_SHA256,
                  chunk_bytes=1048576, workers=4, request_timeout=45,
                  maximum_range_seconds=480, maximum_attempts=4,
                  maximum_attempt_seconds=90, retry_backoff_seconds=1,
                  progress=None, opener_factory=None):
    assert urlparse(url).scheme == 'https' and urlparse(url).hostname == 'zenodo.org'
    assert type(expected_bytes) is int and 0 < expected_bytes <= 2*1024**3
    assert 1024 <= chunk_bytes <= 16*1024**2 and 1 <= workers <= 4
    assert 1 <= maximum_attempts <= 4 and 0 < request_timeout <= 60
    assert 0 < maximum_attempt_seconds <= 120 and 0 <= retry_backoff_seconds <= 2
    spans = [(s, min(expected_bytes-1, s+chunk_bytes-1)) for s in range(0, expected_bytes, chunk_bytes)]
    assert len(spans) <= 2048
    factory = opener_factory or (lambda: build_opener(ProxyHandler({})))
    abort = threading.Event();lock = threading.Lock();root_errors = []

    def read_span(start, end):
        prefix = io.BytesIO();began = time.monotonic();etag = None
        continuation_bytes = 0
        for attempt in range(1, maximum_attempts+1):
            if abort.is_set():
                raise StreamInputError('CANCELLED')
            if time.monotonic()-began > maximum_range_seconds:
                raise StreamInputError('RANGE_ELAPSED_BUDGET')
            offset = start+prefix.tell();remaining = end-offset+1
            if attempt > 1:
                continuation_bytes += prefix.tell()
            headers = {'Range':f'bytes={offset}-{end}'}
            if etag is not None:
                headers['If-Match'] = etag
            qualified = False
            deadline = min(began+maximum_range_seconds,time.monotonic()+maximum_attempt_seconds)
            try:
                with factory().open(Request(url,headers=headers),timeout=request_timeout) as response:
                    assert_url = urlparse(response.url)
                    if response.status != 206 or assert_url.scheme != 'https' or assert_url.hostname != 'zenodo.org':
                        raise StreamInputError('NOT_EXACT_OFFICIAL_206')
                    if response.headers.get('Content-Range') != f'bytes {offset}-{end}/{expected_bytes}':
                        raise StreamInputError('CONTENT_RANGE_MISMATCH')
                    declared = response.headers.get('Content-Length')
                    if declared is not None and int(declared) != remaining:
                        raise StreamInputError('CONTENT_LENGTH_MISMATCH')
                    observed = response.headers.get('ETag')
                    if etag is not None and observed is not None and observed != etag:
                        raise StreamInputError('OBJECT_ETAG_CHANGED')
                    if observed is not None:
                        if len(observed)>256 or '\r' in observed or '\n' in observed:
                            raise StreamInputError('INVALID_ETAG_HEADER')
                        if etag is None and observed.startswith('"') and observed.endswith('"'):
                            etag = observed
                    qualified = True
                    while prefix.tell() < end-start+1:
                        if abort.is_set():
                            raise StreamInputError('CANCELLED')
                        if time.monotonic()>deadline:
                            raise StreamInputError('RANGE_ATTEMPT_ELAPSED_BUDGET')
                        block = response.read1(min(65536,end-start+2-prefix.tell()))
                        if not block:
                            break
                        if prefix.tell()+len(block)>end-start+1:
                            raise StreamInputError('RANGE_LENGTH_EXCESS')
                        prefix.write(block)
                    if time.monotonic()-began>maximum_range_seconds:
                        raise StreamInputError('RANGE_ELAPSED_BUDGET')
                    if prefix.tell()==end-start+1:
                        if declared is None and response.read1(1):
                            raise StreamInputError('RANGE_LENGTH_EXCESS')
                        return prefix.getvalue(),attempt,continuation_bytes
                    raise StreamInputError('RANGE_LENGTH_SHORT')
            except HTTPError:
                raise
            except StreamInputError as error:
                if error.code_name not in ('RANGE_LENGTH_SHORT','RANGE_ATTEMPT_ELAPSED_BUDGET'):
                    raise
                if attempt==maximum_attempts:
                    raise
            except IncompleteRead as error:
                if qualified and error.partial:
                    if prefix.tell()+len(error.partial)>end-start+1:
                        raise StreamInputError('RANGE_LENGTH_EXCESS') from None
                    prefix.write(error.partial)
                if attempt==maximum_attempts and prefix.tell()<end-start+1:
                    raise
            except URLError as error:
                if not isinstance(error.reason,(TimeoutError,ConnectionError,ssl.SSLEOFError)):
                    raise
                if attempt==maximum_attempts and prefix.tell()<end-start+1:
                    raise
            except (TimeoutError,ConnectionError,ssl.SSLEOFError):
                if attempt==maximum_attempts and prefix.tell()<end-start+1:
                    raise
            if prefix.tell()==end-start+1:
                return prefix.getvalue(),attempt,continuation_bytes
            if abort.wait(retry_backoff_seconds*2**(attempt-1)):
                raise StreamInputError('CANCELLED')
        raise AssertionError('UNREACHABLE_RANGE_LOOP')

    def guarded_read_span(start,end):
        try:
            return read_span(start,end)
        except BaseException as error:
            with lock:
                if not root_errors:
                    root_errors.append(error)
            abort.set()
            raise

    pool = ThreadPoolExecutor(max_workers=workers);pending = {}
    md5, sha = hashlib.md5(),hashlib.sha256()
    count = attempts = continued = preserved = max_attempts = 0
    try:
        for index in range(min(workers,len(spans))):
            pending[index] = pool.submit(guarded_read_span,*spans[index])
        for index in range(len(spans)):
            future = pending.pop(index)
            try:
                block, tries, kept = future.result()
            except BaseException:
                if root_errors:
                    raise root_errors[0] from None
                raise
            assert len(block)==spans[index][1]-spans[index][0]+1
            if sink.write(block) != len(block):
                raise StreamInputError('SINK_SHORT_WRITE')
            md5.update(block);sha.update(block);count += len(block)
            attempts += tries;continued += tries>1;preserved += kept;max_attempts = max(max_attempts,tries)
            del block,future
            next_index = index+workers
            if next_index<len(spans):
                pending[next_index] = pool.submit(guarded_read_span,*spans[next_index])
            if progress:
                progress(count)
        assert count==expected_bytes
        if md5.hexdigest()!=expected_MD5:
            raise StreamInputError('WHOLE_OBJECT_MD5_MISMATCH')
        if sha.hexdigest()!=expected_SHA256:
            raise StreamInputError('WHOLE_OBJECT_SHA256_MISMATCH')
        return dict(complete_ranges=len(spans),workers=workers,exact_ranges_and_lengths=True,
                    no_gaps_or_overlaps=True,whole_official_MD5_verified=True,
                    whole_known_SHA256_verified=True,no_range_payload_files=True,
                    no_alternate_endpoint_or_authentication_protocol_hash_retries=True,
                    bounded_verified_suffix_continuation_only=True,
                    ordered_stream_no_local_file_sized_buffer=True,
                    maximum_pending_ranges=workers,range_chunk_bytes=chunk_bytes,
                    actual_request_attempts=attempts,logical_ranges_requiring_continuation=continued,
                    prefix_bytes_not_redownloaded_across_continuations=preserved,
                    maximum_actual_attempts_per_range=max_attempts),md5.hexdigest(),sha.hexdigest()
    finally:
        abort.set()
        for future in pending.values():
            future.cancel()
        pool.shutdown(wait=True,cancel_futures=True)
