"""Same-origin default-TLS persistent small ranges for a fixed public object.

Observed ~17KiB/s needs a 75min technical budget for the selected member.
No authentication, proxy, host/DNS overrides or full-payload fallback. Every
attempt counts against the requested-byte budget, including failed reads.
"""
import http.client
import signal
import ssl
import time
import urllib.parse


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
        self.context=ssl.create_default_context()
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
                    self.connection=http.client.HTTPSConnection(self.host,timeout=15,context=self.context)
                self.connection.request('GET',self.path,headers={'Range':f'bytes={start}-{end}',
                    'Accept-Encoding':'identity','User-Agent':'Python-urllib/3.10'})
                response=self.connection.getresponse()
                assert response.status==206,'no_status_or_redirect_or_full_archive_fallback'
                assert response.getheader('Content-Range')==f'bytes {start}-{end}/{self.size}'
                identity=(response.getheader('ETag'),response.getheader('Last-Modified'))
                if self.identity is None:self.identity=identity
                else:assert self.identity==identity
                raw=response.read(n+1);received=len(raw)
                if received<n:raise http.client.IncompleteRead(raw,n-received)
                if received>n:raise ResponseLengthError('overlong_range_response')
                self.position+=n;self.transferred+=n;self.progress()
                return raw
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
