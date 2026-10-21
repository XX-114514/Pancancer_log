"""RAM-only aligned ranges for one fixed public archive, not a full download.

Explicit larger budget is separate from the closed small-object transports.
Existing loopback proxy is request-scoped. CONNECT/SNI/default TLS retain the
original public host. No redirects, status-200 fallback, disk cache or secrets.
"""
from collections import OrderedDict
import http.client
import signal
import ssl
import time
from urllib.parse import urlsplit


class IntegrityError(Exception):
    pass


class BudgetError(Exception):
    pass


class RangeStatusError(Exception):
    def __init__(self, status):
        self.status = status
        super().__init__('range_status_rejected')


class RangeFile:
    BLOCK = 1024**2
    MAX_REQUESTED = 768*1024**2
    MAX_SECONDS = 7200

    def __init__(self, url, size, port, checkpoint=None, connection_factory=None):
        parsed = urlsplit(url)
        if (parsed.scheme != 'https' or parsed.hostname != 'zenodo.org'
                or parsed.port not in (None,443) or parsed.username is not None
                or parsed.password is not None or parsed.fragment or parsed.query
                or not parsed.path.startswith('/api/records/7324285/files/')):
            raise IntegrityError('source_scope_rejected')
        if not isinstance(port,int) or not 1<=port<=65535 or size<=0:
            raise IntegrityError('transport_parameters_rejected')
        self.path=parsed.path;self.size=size;self.port=port;self.position=0
        self.requested=0;self.transferred=0;self.requests=0;self.failures=[]
        self.started=time.monotonic();self.connection=None;self.identity=None
        self.cache=OrderedDict();self.checkpoint=checkpoint
        self.connection_factory=connection_factory

    def _connect(self):
        if self.connection_factory:return self.connection_factory()
        context=ssl.create_default_context()
        if not context.check_hostname or context.verify_mode!=ssl.CERT_REQUIRED:
            raise IntegrityError('default_TLS_required')
        c=http.client.HTTPSConnection('127.0.0.1',self.port,timeout=35,context=context)
        c.set_tunnel('zenodo.org',443)
        return c

    def close(self):
        if self.connection:self.connection.close();self.connection=None

    def tell(self):return self.position
    def seekable(self):return True
    def seek(self, offset, whence=0):
        if whence not in (0,1,2):raise ValueError('seek_whence')
        p=offset+(0 if whence==0 else self.position if whence==1 else self.size)
        if not 0<=p<=self.size:raise ValueError('seek_bounds')
        self.position=p;return p

    def progress(self):
        if self.checkpoint:self.checkpoint(dict(requested_bytes_including_retries=self.requested,
            successful_range_bytes=self.transferred,requests=self.requests,
            failure_classes=self.failures,elapsed_seconds=time.monotonic()-self.started))

    def _fetch(self, start, n):
        end=start+n-1
        for attempt in range(3):
            if self.requested+n>self.MAX_REQUESTED or time.monotonic()-self.started>self.MAX_SECONDS:
                raise BudgetError('bounded_member_budget_exceeded')
            self.requested+=n;self.requests+=1
            if signal.getitimer(signal.ITIMER_REAL)!=(0.0,0.0):raise IntegrityError('existing_alarm')
            old=signal.getsignal(signal.SIGALRM)
            def deadline(signum,frame):raise TimeoutError('request_deadline')
            signal.signal(signal.SIGALRM,deadline);signal.alarm(45)
            delay=attempt+1;retry=False
            try:
                if self.connection is None:self.connection=self._connect()
                self.connection.request('GET',self.path,headers={
                    'Range':f'bytes={start}-{end}','Accept-Encoding':'identity',
                    'User-Agent':'Python-urllib/3.13'})
                r=self.connection.getresponse()
                if r.status!=206:
                    status=r.status;after=r.getheader('Retry-After');r.close()
                    if status not in (429,500,502,503,504):raise RangeStatusError(status)
                    if after is not None:
                        try:d=int(after)
                        except (ValueError,TypeError):raise RangeStatusError(status)
                        if not 0<=d<=30:raise RangeStatusError(status)
                        delay=max(delay,d)
                    self.failures.append(dict(error_class='TransientHTTPStatus',status=status))
                    retry=True;self.close()
                else:
                    if r.getheader('Content-Range')!=f'bytes {start}-{end}/{self.size}':
                        raise IntegrityError('content_range_mismatch')
                    if r.getheader('Content-Encoding') not in (None,'identity'):
                        raise IntegrityError('content_encoding_mismatch')
                    if r.getheader('Content-Length') is not None and r.getheader('Content-Length')!=str(n):
                        raise IntegrityError('content_length_mismatch')
                    identity=(r.getheader('ETag'),r.getheader('Last-Modified'))
                    # This fixed public endpoint can omit both validators.
                    # Missing is recorded, never invented. The caller MUST
                    # check pinned directories before/after and full member CRC.
                    if self.identity is None:self.identity=identity
                    elif self.identity!=identity:raise IntegrityError('object_identity_changed')
                    raw=r.read(n+1)
                    if len(raw)!=n:raise IntegrityError('range_body_length_mismatch')
                    self.transferred+=n;self.progress();return raw
            except ssl.SSLError:
                self.close();raise
            except (TimeoutError,ConnectionResetError,BrokenPipeError,
                    http.client.IncompleteRead,http.client.RemoteDisconnected) as error:
                self.failures.append(dict(error_class=type(error).__name__))
                retry=True;self.close()
            except Exception:
                self.close();raise
            finally:
                signal.alarm(0);signal.signal(signal.SIGALRM,old)
            self.progress()
            if not retry or attempt==2:raise IntegrityError('range_retry_exhausted')
            time.sleep(delay)
        raise IntegrityError('unreachable')

    def read(self,n=-1):
        if n<0:n=self.size-self.position
        n=min(n,self.size-self.position)
        if not 0<=n<=16*1024**2:raise BudgetError('individual_read_size')
        parts=[]
        while n:
            start=self.position//self.BLOCK*self.BLOCK
            if start not in self.cache:
                self.cache[start]=self._fetch(start,min(self.BLOCK,self.size-start))
                if len(self.cache)>4:self.cache.popitem(last=False)
            self.cache.move_to_end(start);block=self.cache[start]
            offset=self.position-start;take=min(n,len(block)-offset)
            parts.append(block[offset:offset+take]);self.position+=take;n-=take
        return b''.join(parts)
