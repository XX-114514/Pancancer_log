"""Bounded verified suffix continuation after short public HTTP range bodies.

Only RAM prefixes with valid official 206/range/length headers are retained.
The unchanged core checks complete logical ranges, whole coverage and MD5.
HTTP/authentication/protocol/hash errors are terminal; no alternate endpoints.
"""
from urllib.request import Request, build_opener, ProxyHandler
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


def _header(headers, name):
    for key, value in headers.items():
        if key.lower() == name.lower():
            return value
    return None


class _LogicalRange(io.BytesIO):
    def __init__(self, value, url, start, end, total):
        super().__init__(value)
        self.url, self.status = url, 206
        # These describe a reassembled logical range, not a single HTTP response.
        self.headers = {'Content-Range': f'bytes {start}-{end}/{total}',
                        'Content-Length': str(end-start+1)}


def fetch_ranges(url, *, expected_bytes, expected_MD5, chunk_bytes, workers,
                 request_timeout, maximum_range_seconds, maximum_attempts=4,
                 retry_backoff_seconds=1, progress=None, opener_factory=None,
                 maximum_attempt_seconds=90):
    assert 1 <= maximum_attempts <= 4 and 0 <= retry_backoff_seconds <= 2
    assert 0 < maximum_attempt_seconds <= 120
    factory = opener_factory or (lambda: build_opener(ProxyHandler({})))
    abort, lock = threading.Event(), threading.Lock()
    attempts, continued, preserved = {}, set(), 0
    root_errors = []

    def reject(error):
        with lock:
            if not root_errors:
                root_errors.append(error)
        abort.set()
        raise error

    class ResumableRangeOpener:
        def open(self, request, timeout):
            nonlocal preserved
            span = request.get_header('Range')
            start, end = map(int, span.removeprefix('bytes=').split('-'))
            began = time.monotonic();deadline = began+maximum_range_seconds
            length = end-start+1;buffer = io.BytesIO();etag = None
            attempt = 0

            def annotate(error):
                error.range_start, error.range_end = start, end
                error.expected_range_bytes, error.received_range_bytes = length, buffer.tell()
                error.range_attempt_number = attempt
                error.range_elapsed_seconds = max(0, time.monotonic()-began)
                return error

            for attempt in range(1, maximum_attempts+1):
                if abort.is_set():
                    raise RangeInputError('CANCELLED')
                if time.monotonic() > deadline:
                    reject(annotate(RangeInputError('RANGE_ELAPSED_BUDGET')))
                if buffer.tell() == length:
                    return _LogicalRange(buffer.getvalue(), url, start, end, expected_bytes)
                offset = start+buffer.tell();remaining = end-offset+1
                headers = {'Range': f'bytes={offset}-{end}'}
                if etag is not None:
                    headers['If-Match'] = etag
                with lock:
                    attempts[span] = attempts.get(span, 0)+1
                    if offset > start:
                        continued.add(span);preserved += buffer.tell()
                attempt_deadline = min(deadline, time.monotonic()+maximum_attempt_seconds)
                headers_qualified = False
                try:
                    with factory().open(Request(url, headers=headers), timeout=timeout) as response:
                        if response.status != 206 or urlparse(response.url).scheme != 'https' or urlparse(response.url).hostname != 'zenodo.org':
                            raise RangeInputError('NOT_EXACT_OFFICIAL_206')
                        if _header(response.headers, 'Content-Range') != f'bytes {offset}-{end}/{expected_bytes}':
                            raise RangeInputError('CONTENT_RANGE_MISMATCH')
                        declared = _header(response.headers, 'Content-Length')
                        if declared is not None and int(declared) != remaining:
                            raise RangeInputError('CONTENT_LENGTH_MISMATCH')
                        observed_etag = _header(response.headers, 'ETag')
                        if etag is not None and observed_etag is not None and observed_etag != etag:
                            raise RangeInputError('OBJECT_ETAG_CHANGED')
                        if etag is None and observed_etag is not None:
                            if '\r' in observed_etag or '\n' in observed_etag or len(observed_etag) > 256:
                                raise RangeInputError('INVALID_ETAG_HEADER')
                            if observed_etag.startswith('"') and observed_etag.endswith('"'):
                                etag = observed_etag
                        headers_qualified = True
                        while buffer.tell() < length:
                            if abort.is_set():
                                raise RangeInputError('CANCELLED')
                            if time.monotonic() > deadline:
                                raise RangeInputError('RANGE_ELAPSED_BUDGET')
                            if time.monotonic() > attempt_deadline:
                                raise RangeInputError('RANGE_ATTEMPT_ELAPSED_BUDGET')
                            value = response.read1(min(65536, length-buffer.tell()+1))
                            if not value:
                                break
                            if buffer.tell()+len(value) > length:
                                raise RangeInputError('RANGE_LENGTH_EXCESS')
                            buffer.write(value)
                        if time.monotonic() > deadline:
                            raise RangeInputError('RANGE_ELAPSED_BUDGET')
                        if buffer.tell() == length:
                            # Declared framing or EOF must agree with the exact suffix.
                            if declared is None and response.read1(1):
                                raise RangeInputError('RANGE_LENGTH_EXCESS')
                            return _LogicalRange(buffer.getvalue(), url, start, end, expected_bytes)
                        raise RangeInputError('RANGE_LENGTH_SHORT')
                except HTTPError as error:
                    reject(annotate(error))
                except RangeInputError as error:
                    if error.code_name not in ['RANGE_LENGTH_SHORT', 'RANGE_ATTEMPT_ELAPSED_BUDGET'] or attempt == maximum_attempts:
                        reject(annotate(error))
                except IncompleteRead as error:
                    if headers_qualified and error.partial:
                        if buffer.tell()+len(error.partial) > length:
                            reject(annotate(RangeInputError('RANGE_LENGTH_EXCESS')))
                        buffer.write(error.partial)
                    if attempt == maximum_attempts and buffer.tell() < length:
                        reject(annotate(error))
                except (TimeoutError, ConnectionError, ssl.SSLEOFError) as error:
                    if attempt == maximum_attempts and buffer.tell() < length:
                        reject(annotate(error))
                except URLError as error:
                    transient = isinstance(error.reason, (TimeoutError, ConnectionError, ssl.SSLEOFError))
                    if not transient or (attempt == maximum_attempts and buffer.tell() < length):
                        reject(annotate(error))
                except BaseException as error:
                    reject(annotate(error))
                if buffer.tell() == length:
                    if time.monotonic() > deadline:
                        reject(annotate(RangeInputError('RANGE_ELAPSED_BUDGET')))
                    return _LogicalRange(buffer.getvalue(), url, start, end, expected_bytes)
                if abort.wait(retry_backoff_seconds*(2**(attempt-1))):
                    raise RangeInputError('CANCELLED')
            raise AssertionError('UNREACHABLE_CONTINUATION_LOOP')

    try:
        result = _core['fetch_ranges'](url, expected_bytes=expected_bytes,
            expected_MD5=expected_MD5, chunk_bytes=chunk_bytes, workers=workers,
            request_timeout=request_timeout, maximum_range_seconds=maximum_range_seconds,
            progress=progress, opener_factory=lambda: ResumableRangeOpener())
    except BaseException as error:
        if root_errors and not isinstance(error, (TimeoutError, KeyboardInterrupt, SystemExit)):
            raise root_errors[0] from None
        raise
    buffer, md5, sha, proof = result
    proof.pop('no_retries_or_alternate_endpoints')
    proof.update(bounded_verified_suffix_continuation_only=True,
        logical_headers_are_reassembled_not_single_HTTP_response=True,
        all_accepted_prefixes_had_qualified_206_range_and_length_headers=True,
        no_alternate_endpoint_or_authentication_protocol_hash_retries=True,
        maximum_attempts_per_logical_range=maximum_attempts,
        configured_attempt_soft_deadline_seconds=maximum_attempt_seconds,
        logical_range_budget_seconds=maximum_range_seconds,
        per_socket_blocking_grace_seconds=request_timeout,
        actual_request_attempts=sum(attempts.values()),
        logical_ranges_requiring_continuation=len(continued),
        prefix_bytes_not_redownloaded_across_continuations=preserved,
        maximum_actual_attempts_per_range=max(attempts.values(), default=0))
    return buffer, md5, sha, proof
