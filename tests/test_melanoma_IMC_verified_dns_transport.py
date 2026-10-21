"""Verify same-host DNS selection never accepts an invalid TLS connection."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import socket
import ssl
import unittest
from unittest.mock import patch


class VerifiedAddressSelection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[1]/'scripts'/'melanoma_IMC_verified_dns_transport_v1.py'
        spec = importlib.util.spec_from_file_location('verified_addresses', path)
        cls.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.m)

    def setUp(self):
        self.m._preferred = None
        self.m._diagnostics.clear()
        self.records = [(socket.AF_INET, socket.SOCK_STREAM, 6, '', ('192.0.2.1', 443)),
            (socket.AF_INET, socket.SOCK_STREAM, 6, '', ('192.0.2.2', 443))]
        self.created = []
        class Raw:
            def __init__(instance, *args):
                instance.closed = False
                self.created.append(instance)
            def settimeout(instance, value):pass
            def connect(instance, address):instance.address = address
            def close(instance):instance.closed = True
        self.raw = Raw

    def client(self, wrap):
        result = self.m.VerifiedDNSHTTPSConnection('zenodo.org', timeout=35)
        result._context = SimpleNamespace(check_hostname=True, verify_mode=ssl.CERT_REQUIRED, wrap_socket=wrap)
        return result

    def test_bad_certificate_rejected_before_valid_same_host_candidate(self):
        secured = SimpleNamespace(settimeout=lambda value: None)
        def wrap(raw, server_hostname):
            self.assertEqual(server_hostname, 'zenodo.org')
            if raw.address[0].endswith('.1'):
                raise ssl.SSLCertVerificationError(1, 'synthetic bad hostname')
            return secured
        client = self.client(wrap)
        with patch.object(self.m.socket, 'getaddrinfo', return_value=self.records), patch.object(self.m.socket, 'socket', self.raw):
            client.connect()
        self.assertIs(client.sock, secured)
        self.assertTrue(self.created[0].closed)
        self.assertEqual(self.m._diagnostics['SSLCertVerificationError'], 1)
        self.assertEqual(self.m._diagnostics['verified_connections'], 1)

    def test_all_bad_certificates_abort_no_accepted_socket(self):
        def wrap(raw, server_hostname):raise ssl.SSLCertVerificationError(1, 'synthetic')
        client = self.client(wrap)
        with patch.object(self.m.socket, 'getaddrinfo', return_value=self.records), patch.object(self.m.socket, 'socket', self.raw):
            with self.assertRaises(ssl.SSLCertVerificationError):client.connect()
        self.assertIsNone(client.sock)
        self.assertTrue(all(raw.closed for raw in self.created))

    def test_outer_hard_deadline_is_not_swallowed(self):
        def wrap(raw, server_hostname):raise TimeoutError('bounded_request_deadline')
        client = self.client(wrap)
        with patch.object(self.m.socket, 'getaddrinfo', return_value=self.records), patch.object(self.m.socket, 'socket', self.raw):
            with self.assertRaisesRegex(TimeoutError, 'bounded_request_deadline'):client.connect()
        self.assertEqual(len(self.created), 1)

    def test_no_stale_preferred_address_outside_current_DNS(self):
        self.m._preferred = (socket.AF_INET, ('198.51.100.10', 443))
        client = self.client(lambda raw, server_hostname: SimpleNamespace(settimeout=lambda value: None))
        with patch.object(self.m.socket, 'getaddrinfo', return_value=self.records), patch.object(self.m.socket, 'socket', self.raw):
            client.connect()
        self.assertEqual(self.created[0].address, ('192.0.2.1', 443))

    def test_hostname_verification_cannot_be_disabled(self):
        client = self.client(lambda raw, server_hostname: None)
        client._context.check_hostname = False
        with patch.object(self.m.socket, 'getaddrinfo') as dns:
            with self.assertRaises(AssertionError):client.connect()
        dns.assert_not_called()
        self.assertIsNot(self.m.http.client.HTTPSConnection, self.m.VerifiedDNSHTTPSConnection)


if __name__ == '__main__':
    unittest.main()
