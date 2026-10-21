"""Only transient HTTP transport failure changes current-DNS peer priority."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import socket
import http.client

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('dns_fixtures',ROOT/'tests/test_melanoma_IMC_verified_dns_transport.py')
fixtures=importlib.util.module_from_spec(spec);spec.loader.exec_module(fixtures)


class HTTPPeerFailover(fixtures.VerifiedAddressSelection):
    @classmethod
    def setUpClass(cls):
        spec=importlib.util.spec_from_file_location('http_peer_backend',ROOT/'scripts/melanoma_IMC_verified_dns_transport_v3.py')
        cls.m=importlib.util.module_from_spec(spec);spec.loader.exec_module(cls.m)
    def setUp(self):
        super().setUp();self.m._failed_http_peers.clear()
    def test_body_timeout_demotes_peer_and_preserves_original_exception(self):
        client=self.client(lambda raw,server_hostname:SimpleNamespace(settimeout=lambda value:None))
        client._selected_peer=(socket.AF_INET,self.records[0][4]);self.m._preferred=client._selected_peer
        error=TimeoutError('synthetic_body_timeout')
        def read(*args):raise error
        with patch.object(http.client.HTTPSConnection,'getresponse',return_value=SimpleNamespace(status=206,read=read)):
            response=client.getresponse()
            with self.assertRaises(TimeoutError) as caught:response.read(10)
        self.assertIs(caught.exception,error)
        fresh=self.client(lambda raw,server_hostname:SimpleNamespace(settimeout=lambda value:None))
        with patch.object(self.m.socket,'getaddrinfo',return_value=self.records),patch.object(self.m.socket,'socket',self.raw):fresh.connect()
        self.assertEqual(self.created[0].address,self.records[1][4])
    def test_header_timeout_demotes_only_current_peer(self):
        client=self.client(lambda *args:None);client._selected_peer=(socket.AF_INET,self.records[0][4])
        with patch.object(http.client.HTTPSConnection,'getresponse',side_effect=TimeoutError('synthetic_headers')):
            with self.assertRaises(TimeoutError):client.getresponse()
        self.assertEqual(self.m._failed_http_peers,{client._selected_peer})
    def test_permission_or_redirect_status_does_not_rotate(self):
        for status in [301,403,429]:
            client=self.client(lambda *args:None);client._selected_peer=(socket.AF_INET,self.records[0][4])
            with patch.object(http.client.HTTPSConnection,'getresponse',return_value=SimpleNamespace(status=status,read=lambda *args:b'')):
                self.assertEqual(client.getresponse().status,status)
            self.assertFalse(self.m._failed_http_peers)


if __name__=='__main__':unittest.main()
