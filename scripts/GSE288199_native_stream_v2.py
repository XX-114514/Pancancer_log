"""Same reader via a task-owned SSH Unix socket to official HTTPS.

TLS hostname/certificate checks remain on the server. SSH opens the official
TCP connection from the user's machine. No TCP listener or proxy fallback.
"""
import http.client
import os
from pathlib import Path
import socket
import ssl
import urllib.request

import GSE288199_native_stream_v1 as base


class OfficialHTTPSConnection(http.client.HTTPSConnection):
    def connect(self):
        assert self.host=='ftp.ncbi.nlm.nih.gov' and self.port==443
        assert not self._tunnel_host
        path=Path(os.environ['PANCANCER_PUBLIC_HTTPS_RELAY_SOCKET'])
        assert path.parent.parent==Path('/dev/shm') and path.is_socket()
        raw=socket.socket(socket.AF_UNIX,socket.SOCK_STREAM)
        raw.settimeout(self.timeout)
        try:
            raw.connect(str(path))
            self.sock=self._context.wrap_socket(raw,server_hostname=self.host)
        except BaseException:
            raw.close();raise


class OfficialHTTPSHandler(urllib.request.HTTPSHandler):
    def https_open(self,request):
        return self.do_open(OfficialHTTPSConnection,request,context=self._context,
                            check_hostname=self._check_hostname)


def opener():
    context=ssl.create_default_context()
    assert context.check_hostname and context.verify_mode==ssl.CERT_REQUIRED
    return urllib.request.build_opener(urllib.request.ProxyHandler({}),OfficialHTTPSHandler(context=context))


# The origin, parsing, numeric gates, counters, CRC and digests are unchanged.
base.opener=opener
BASE=base.BASE
PREFIX=base.PREFIX
feature_contract=base.feature_contract
census=base.census
header=base.header
chunks=base.chunks
consume=base.consume
read_full=base.read_full
synthetic=base.synthetic
