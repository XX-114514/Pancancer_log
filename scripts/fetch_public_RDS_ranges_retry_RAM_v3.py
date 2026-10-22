"""Bounded same-range retries for truncated/transient public input only.

The accepted core still verifies exact ranges, whole coverage and official MD5.
HTTP/authentication/protocol/hash errors are never retried. Nothing is persisted.
"""
from urllib.request import build_opener, ProxyHandler
from urllib.error import HTTPError, URLError
from urllib.parse import urlparse
from http.client import IncompleteRead
import io
import ssl
import threading
import time

_core = {}
exec(RANGE_CORE_SOURCE, _core)
RangeInputError = _core['RangeInputError']


class _CompleteRange(io.BytesIO):
    def __init__(self, value, url, headers):
        super().__init__(value)
        self.url, self.headers, self.status = url, headers, 206


def fetch_ranges(url, *, expected_bytes, expected_MD5, chunk_bytes, workers,
                 request_timeout, maximum_range_seconds, maximum_attempts=4,
                 retry_backoff_seconds=1, progress=None, opener_factory=None):
    assert 1 <= maximum_attempts <= 4 and 0 <= retry_backoff_seconds <= 2
    factory = opener_factory or (lambda: build_opener(ProxyHandler({})))
    abort, lock, attempts = threading.Event(), threading.Lock(), {}
    root_errors = []

    def reject(error):
        with lock:
            if not root_errors:
                root_errors.append(error)
        abort.set()
        raise error

    class CompleteRangeOpener:
        def open(self, request, timeout):
            span = request.get_header('Range')
            start, end = map(int, span.removeprefix('bytes=').split('-'))
            length, deadline = end-start+1, time.monotonic()+maximum_range_seconds
            for attempt in range(1, maximum_attempts+1):
                if abort.is_set():
                    raise RangeInputError('CANCELLED')
                if time.monotonic() > deadline:
                    reject(RangeInputError('RANGE_ELAPSED_BUDGET'))
                with lock:
                    attempts[span] = attempts.get(span, 0)+1
                try:
                    with factory().open(request, timeout=timeout) as response:
                        if response.status != 206 or urlparse(response.url).scheme != 'https' or urlparse(response.url).hostname != 'zenodo.org':
                            raise RangeInputError('NOT_EXACT_OFFICIAL_206')
                        if response.headers.get('Content-Range') != f'bytes {start}-{end}/{expected_bytes}':
                            raise RangeInputError('CONTENT_RANGE_MISMATCH')
                        declared = response.headers.get('Content-Length')
                        if declared is not None and int(declared) != length:
                            raise RangeInputError('CONTENT_LENGTH_MISMATCH')
                        buffer = io.BytesIO()
                        while True:
                            if abort.is_set():
                                raise RangeInputError('CANCELLED')
                            if time.monotonic() > deadline:
                                raise RangeInputError('RANGE_ELAPSED_BUDGET')
                            value = response.read1(min(65536, length-buffer.tell()+1))
                            if not value:
                                break
                            buffer.write(value)
                            if buffer.tell() > length:
                                raise RangeInputError('RANGE_LENGTH_EXCESS')
                        if buffer.tell() != length:
                            error = RangeInputError('RANGE_LENGTH_SHORT')
                            error.range_start, error.range_end = start, end
                            error.expected_range_bytes, error.received_range_bytes = length, buffer.tell()
                            raise error
                        return _CompleteRange(buffer.getvalue(), response.url, response.headers)
                except HTTPError as error:
                    reject(error)
                except RangeInputError as error:
                    if error.code_name != 'RANGE_LENGTH_SHORT' or attempt == maximum_attempts:
                        reject(error)
                except (TimeoutError, ConnectionError, IncompleteRead, ssl.SSLEOFError) as error:
                    if attempt == maximum_attempts:
                        reject(error)
                except URLError as error:
                    transient = isinstance(error.reason, (TimeoutError, ConnectionError, ssl.SSLEOFError))
                    if not transient or attempt == maximum_attempts:
                        reject(error)
                except BaseException as error:
                    reject(error)
                if abort.wait(retry_backoff_seconds*(2**(attempt-1))):
                    raise RangeInputError('CANCELLED')
            raise AssertionError('UNREACHABLE_RETRY_LOOP')

    try:
        result = _core['fetch_ranges'](url, expected_bytes=expected_bytes,
            expected_MD5=expected_MD5, chunk_bytes=chunk_bytes, workers=workers,
            request_timeout=request_timeout, maximum_range_seconds=maximum_range_seconds,
            progress=progress, opener_factory=lambda: CompleteRangeOpener())
    except BaseException as error:
        if root_errors and not isinstance(error, (TimeoutError, KeyboardInterrupt, SystemExit)):
            raise root_errors[0] from None
        raise
    buffer, md5, sha, proof = result
    proof.pop('no_retries_or_alternate_endpoints')
    proof.update(no_alternate_endpoint_or_authentication_protocol_hash_retries=True,
        bounded_same_range_transport_retries_only=True, maximum_attempts_per_range=maximum_attempts,
        actual_request_attempts=sum(attempts.values()), ranges_requiring_retry=sum(v > 1 for v in attempts.values()),
        maximum_actual_attempts_per_range=max(attempts.values(), default=0))
    return buffer, md5, sha, proof
