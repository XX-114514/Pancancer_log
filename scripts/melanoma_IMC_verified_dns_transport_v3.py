"""Strict same-host TLS across current system-DNS candidates only.

Reject a failing address before HTTP; accept only a fully certificate-verified
connection for the unchanged HTTPS hostname. No resolver, proxy, trust, SNI,
URL-host or system settings are changed. HTTP transport failures demote that peer on the next current-DNS resolution.
The base request/range/budget guards remain intact. This is address selection, never use of an invalid certificate.
"""
from collections import Counter
from pathlib import Path
import hashlib
import http.client
import importlib.util
import socket
import ssl
from types import SimpleNamespace

path = Path(__file__).with_name('melanoma_IMC_local_bounded_transport_v3.py')
assert hashlib.sha256(path.read_bytes()).hexdigest() == '1d575fb3c460f0df3511f5fda36f79ef266d633566aefc67d4d4187893ca012d'
spec = importlib.util.spec_from_file_location('guarded_verified_dns_base', path)
_base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_base)
_preferred = None
_diagnostics = Counter()
_failed_http_peers = set()


class VerifiedDNSHTTPSConnection(http.client.HTTPSConnection):
    def connect(self):
        global _preferred
        assert self.host == 'zenodo.org' and self.port == 443
        assert self._tunnel_host is None and self.source_address is None
        assert self._context.check_hostname and self._context.verify_mode == ssl.CERT_REQUIRED
        candidates = []
        for family, kind, protocol, canonical, address in socket.getaddrinfo(self.host, self.port, type=socket.SOCK_STREAM):
            candidate = (family, address)
            if candidate not in candidates:
                candidates.append(candidate)
        candidates.sort(key=lambda value: (value in _failed_http_peers, value != _preferred))
        assert candidates
        last_error = None
        for family, address in candidates[:4]:
            raw = None
            try:
                raw = socket.socket(family, socket.SOCK_STREAM)
                raw.settimeout(6)
                raw.connect(address)
                # SNI/hostname validation remains the original URL hostname.
                secured = self._context.wrap_socket(raw, server_hostname=self.host)
                secured.settimeout(self.timeout)
                self.sock = secured
                self._selected_peer = (family, address)
                _preferred = (family, address)
                _diagnostics['verified_connections'] += 1
                return
            except Exception as error:
                if raw is not None:
                    raw.close()
                # Never swallow the outer hard request alarm.
                if isinstance(error, TimeoutError) and error.args == ('bounded_request_deadline',):
                    raise
                if not isinstance(error, (OSError, ssl.SSLError)):
                    raise
                _diagnostics[type(error).__name__] += 1
                last_error = error
        raise last_error


    def _failed_http(self):
        global _preferred
        peer = getattr(self, '_selected_peer', None)
        if peer is not None:
            _failed_http_peers.add(peer)
            if _preferred == peer: _preferred = None
            _diagnostics['HTTP_transport_failures_demoted'] += 1

    def getresponse(self):
        transient = (TimeoutError, ConnectionResetError, BrokenPipeError,
            http.client.IncompleteRead, http.client.RemoteDisconnected)
        try:
            response = super().getresponse()
        except transient:
            self._failed_http()
            raise
        ordinary_read = response.read
        def checked_read(*args, **kwargs):
            try:
                return ordinary_read(*args, **kwargs)
            except transient:
                self._failed_http()
                raise
        response.read = checked_read
        return response


# Dependency injection is scoped to this own module, not the stdlib or OS.
_base.http = SimpleNamespace(client=SimpleNamespace(
    HTTPSConnection=VerifiedDNSHTTPSConnection,
    IncompleteRead=http.client.IncompleteRead,
    RemoteDisconnected=http.client.RemoteDisconnected))


class RangeFile(_base.RangeFile):
    def tls_diagnostics(self):
        return dict(_diagnostics)
