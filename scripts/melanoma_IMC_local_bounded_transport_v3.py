"""Same-origin default-TLS persistent small ranges for a fixed public object.

Local socket timeout is 35s to match the successful API client; the 45s
per-request alarm, strict default TLS, retries, and byte budgets remain.
The standard HTTPS context factory preserves HTTP/1.1 ALPN negotiation.
No authentication, proxy, host/DNS overrides or full-payload fallback. Every
attempt counts against the requested-byte budget, including failed reads.
"""
import http.client
import signal
import ssl
import time
import urllib.parse


class HTTPRangeStatusError(Exception):
    def __init__(self, status):
        self.status = status
        super().__init__("unexpected_HTTP_range_status")


class TransientHTTPRangeStatusError(HTTPRangeStatusError):
    def __init__(self, status, delay):
        super().__init__(status)
        self.delay = delay


class ResponseLengthError(Exception):
    pass


class RangeFile:
    def __init__(self,source,url,size,checkpoint=None):
        parsed=urllib.parse.urlsplit(url)
        assert parsed.scheme=='https' and parsed.hostname=='zenodo.org'
        assert parsed.port in [None,443] and parsed.username is None and not parsed.fragment
        self.host=parsed.hostname;self.path=parsed.path+('?' + parsed.query if parsed.query else '')
        self.size=size;self.position=0;self.transferred=0;self.requested=0;self.requests=0
        self.identity=None;self.failures=[];self.connection=None;self.checkpoint=checkpoint
        self.started=time.monotonic();self.max_chunk=512*1024
    def tell(self):return self.position
    def seek(self,offset,whence=0):
        self.position=offset if whence==0 else self.position+offset if whence==1 else self.size+offset
        assert 0<=self.position<=self.size
        return self.position
    def seekable(self):return True
    def close(self):
        if self.connection is not None:self.connection.close();self.connection=None
    def progress(self):
        if self.checkpoint:
            self.checkpoint(dict(successful_range_bytes=self.transferred,requested_bytes_including_failed=self.requested,
                requests=self.requests,failure_classes=self.failures,elapsed_seconds=time.monotonic()-self.started,
                no_record_or_payload_values=True))
    def _read_chunk(self,n):
        assert 0<n<=self.max_chunk
        start=self.position;end=start+n-1
        for attempt in range(3):
            assert time.monotonic()-self.started<=75*60
            assert self.requested+n<=64*1024**2
            self.requested+=n;self.requests+=1
            received=None
            old_handler=signal.getsignal(signal.SIGALRM)
            assert signal.getitimer(signal.ITIMER_REAL)==(0.0,0.0),'existing_alarm_not_supported'
            def deadline(signum,frame):raise TimeoutError('bounded_request_deadline')
            signal.signal(signal.SIGALRM,deadline);signal.alarm(45)
            try:
                if self.connection is None:
                    self.connection=http.client.HTTPSConnection(self.host,timeout=35)
                    assert self.connection._context.check_hostname
                    assert self.connection._context.verify_mode==ssl.CERT_REQUIRED
                self.connection.request('GET',self.path,headers={'Range':f'bytes={start}-{end}',
                    'Accept-Encoding':'identity','User-Agent':'Python-urllib/3.10'})
                response=self.connection.getresponse()
                if response.status != 206:
                    # Never consume/log an error body, follow a redirect or use
                    # a status200 full archive. Record only a numeric status.
                    status = response.status
                    retry_after = response.getheader('Retry-After')
                    response.close()
                    delay = attempt+1
                    if retry_after is not None:
                        try:
                            parsed_delay = int(retry_after)
                            if not 0 <= parsed_delay <= 30:
                                raise ValueError()
                            delay = max(delay, parsed_delay)
                        except (ValueError, TypeError):
                            self.failures.append(dict(error_class='HTTPRangeStatusError', HTTP_status=status,
                                retry_deferred=True, requested_bytes=n, received_bytes=0))
                            self.progress()
                            raise HTTPRangeStatusError(status)
                    if status in (429, 500, 502, 503, 504):
                        raise TransientHTTPRangeStatusError(status, delay)
                    self.failures.append(dict(error_class='HTTPRangeStatusError', HTTP_status=status,
                        requested_bytes=n, received_bytes=0))
                    self.progress()
                    raise HTTPRangeStatusError(status)
                assert response.getheader('Content-Range')==f'bytes {start}-{end}/{self.size}'
                identity=(response.getheader('ETag'),response.getheader('Last-Modified'))
                if self.identity is None:self.identity=identity
                else:assert self.identity==identity
                raw=response.read(n+1);received=len(raw)
                if received<n:raise http.client.IncompleteRead(raw,n-received)
                if received>n:raise ResponseLengthError('overlong_range_response')
                self.position+=n;self.transferred+=n;self.progress()
                return raw
            except TransientHTTPRangeStatusError as error:
                self.failures.append(dict(error_class=type(error).__name__,HTTP_status=error.status,
                    requested_bytes=n,received_bytes=0,retry_delay_seconds=error.delay))
                self.close();self.progress()
                if attempt==2:raise
                # No sleep under the per-request alarm; all three attempts
                # and sleeps remain within the total75min/64MiB budgets.
                signal.alarm(0)
                time.sleep(error.delay)
            except ssl.SSLError:
                self.close();raise
            except (TimeoutError,ConnectionResetError,BrokenPipeError,http.client.IncompleteRead,http.client.RemoteDisconnected) as error:
                if isinstance(error,http.client.IncompleteRead):received=len(error.partial)
                self.failures.append(dict(error_class=type(error).__name__,requested_bytes=n,received_bytes=received))
                self.close();self.progress()
                if attempt==2:raise
                time.sleep(attempt+1)
            except Exception:
                self.close();raise
            finally:
                signal.alarm(0);signal.signal(signal.SIGALRM,old_handler)
        raise AssertionError('unreachable')
    def read(self,n=-1):
        if n<0:n=self.size-self.position
        n=min(n,self.size-self.position);assert 0<=n<=16*1024**2
        parts=[]
        while n:
            size=min(n,self.max_chunk);parts.append(self._read_chunk(size));n-=size
        return b''.join(parts)
