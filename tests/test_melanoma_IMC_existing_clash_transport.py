"""Request-local proxy configuration cannot change destination or trust rules."""
import importlib.util
from pathlib import Path
import ssl
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('existing_proxy',ROOT/'scripts/melanoma_IMC_existing_clash_transport_v1.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


class ExistingProxy(unittest.TestCase):
    def tearDown(self):m.PROXY_PORT=None
    def test_current_port_only_original_host_connect_and_default_TLS(self):
        m.PROXY_PORT=7897
        c=m.ExistingClashHTTPSConnection('zenodo.org')
        self.assertEqual((c.host,c.port),('127.0.0.1',7897))
        self.assertEqual((c._tunnel_host,c._tunnel_port),('zenodo.org',443))
        self.assertTrue(c._context.check_hostname)
        self.assertEqual(c._context.verify_mode,ssl.CERT_REQUIRED)
        c.close()
    def test_port_must_be_explicit_and_valid(self):
        for port in [None,0,65536,'7897']:
            m.PROXY_PORT=port
            with self.assertRaises(AssertionError):m.ExistingClashHTTPSConnection('zenodo.org')
    def test_other_origins_are_rejected(self):
        m.PROXY_PORT=7897
        with self.assertRaises(AssertionError):m.ExistingClashHTTPSConnection('example.com')
    def test_global_HTTP_client_is_unchanged(self):
        import http.client
        self.assertIsNot(http.client.HTTPSConnection,m.ExistingClashHTTPSConnection)
        self.assertIs(m._base.http.client.HTTPSConnection,m.ExistingClashHTTPSConnection)


if __name__=='__main__':unittest.main()
