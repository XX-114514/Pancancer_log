"""Request-only use of a verified existing loopback Clash HTTP listener.

No system proxy/TUN/DNS/routing/trust edits. The CONNECT destination and TLS
hostname remain zenodo.org. Original status/range/identity and wire limits stay.
The caller supplies the port already read from the current local configuration.
"""
from pathlib import Path
import hashlib
import http.client
import importlib.util
import ssl
from types import SimpleNamespace

PROXY_PORT = None
path=Path(__file__).with_name('melanoma_IMC_local_bounded_transport_v3.py')
assert hashlib.sha256(path.read_bytes()).hexdigest()=='1d575fb3c460f0df3511f5fda36f79ef266d633566aefc67d4d4187893ca012d'
spec=importlib.util.spec_from_file_location('strict_proxy_range_base',path)
_base=importlib.util.module_from_spec(spec);spec.loader.exec_module(_base)


class ExistingClashHTTPSConnection(http.client.HTTPSConnection):
    def __init__(self,host,timeout=35):
        assert host=='zenodo.org'
        assert isinstance(PROXY_PORT,int) and 1<=PROXY_PORT<=65535
        context=ssl.create_default_context()
        assert context.check_hostname and context.verify_mode==ssl.CERT_REQUIRED
        super().__init__('127.0.0.1',port=PROXY_PORT,timeout=timeout,context=context)
        self.set_tunnel(host,443)


_base.http=SimpleNamespace(client=SimpleNamespace(
    HTTPSConnection=ExistingClashHTTPSConnection,
    IncompleteRead=http.client.IncompleteRead,
    RemoteDisconnected=http.client.RemoteDisconnected))


class RangeFile(_base.RangeFile):
    def tls_diagnostics(self):
        return dict(existing_loopback_proxy_request_only=True,
            original_destination_SNI_and_default_certificate_checks=True,
            no_global_network_settings_changed=True)
