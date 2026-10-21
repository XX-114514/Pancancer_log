"""Bounded transient recovery, never permission or integrity bypass."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch


class HTTPStatusRecovery(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        p = Path(__file__).resolve().parents[1]/'scripts/melanoma_IMC_local_bounded_transport_v3.py'
        spec = importlib.util.spec_from_file_location('http_status_test', p)
        cls.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.m)

    def run_fixture(self, codes, retry_after=None, etags=None):
        state = SimpleNamespace(requests=0, error_bodies=0, closed=0, snapshots=[])
        class Connection:
            def __init__(c, *args, **kwargs):
                import ssl
                c._context = SimpleNamespace(check_hostname=True, verify_mode=ssl.CERT_REQUIRED)
            def request(c, *args, **kwargs):
                state.requests += 1
            def close(c):state.closed += 1
            def getresponse(c):
                status = codes[min(state.requests-1, len(codes)-1)]
                class Response:
                    def getheader(r, name):
                        return {'Content-Range':'bytes 0-2/3', 'Retry-After':retry_after,
                            'ETag':etags[state.requests-1] if etags else 'fixed', 'Last-Modified':'fixed'}.get(name)
                    def close(r):pass
                    def read(r, count):
                        if status != 206:state.error_bodies += 1
                        return b'abc'
                r = Response();r.status = status
                return r
        client = self.m.RangeFile(None, 'https://zenodo.org/synthetic', 3,
            checkpoint=lambda v:state.snapshots.append(v))
        return client, state, Connection

    def test_transient_then_exact_range_with_recorded_status(self):
        for code in [429, 500, 502, 503, 504]:
            c, s, factory = self.run_fixture([code,206])
            with patch.object(self.m.http.client, 'HTTPSConnection', factory), patch.object(self.m.time, 'sleep'):
                self.assertEqual(c.read(3), b'abc')
            self.assertEqual((s.requests,c.requested,c.transferred,c.position,s.error_bodies),(2,6,3,3,0))
            self.assertEqual(c.failures[0]['HTTP_status'],code)

    def test_permission_redirect_and_full_archive_are_not_retried(self):
        for code in [200,301,302,401,403,404,416]:
            c,s,f = self.run_fixture([code,206])
            with patch.object(self.m.http.client,'HTTPSConnection',f), patch.object(self.m.time,'sleep') as sleep:
                with self.assertRaises(self.m.HTTPRangeStatusError):c.read(3)
            self.assertEqual(s.requests,1);sleep.assert_not_called()
            self.assertEqual((c.position,s.error_bodies),(0,0))
            self.assertEqual(c.failures[0]['HTTP_status'],code)

    def test_retries_have_three_attempt_limit(self):
        c,s,f=self.run_fixture([503])
        with patch.object(self.m.http.client,'HTTPSConnection',f),patch.object(self.m.time,'sleep'):
            with self.assertRaises(self.m.TransientHTTPRangeStatusError):c.read(3)
        self.assertEqual((s.requests,c.requested,c.transferred,c.position),(3,9,0,0))

    def test_short_retry_after_is_respected_and_long_or_date_is_deferred(self):
        c,s,f=self.run_fixture([429,206],retry_after='12')
        with patch.object(self.m.http.client,'HTTPSConnection',f),patch.object(self.m.time,'sleep') as sleep:
            self.assertEqual(c.read(3),b'abc')
        sleep.assert_called_once_with(12)
        for value in ['31','-1','Wed, 21 Oct 2026 10:00:00 GMT']:
            c,s,f=self.run_fixture([429,206],retry_after=value)
            with patch.object(self.m.http.client,'HTTPSConnection',f),patch.object(self.m.time,'sleep') as sleep:
                with self.assertRaises(self.m.HTTPRangeStatusError):c.read(3)
            self.assertEqual(s.requests,1);sleep.assert_not_called()
            self.assertTrue(c.failures[0]['retry_deferred'])

    def test_content_identity_change_remains_fatal(self):
        c,s,f=self.run_fixture([206],etags=['changed']);c.identity=('original','fixed')
        with patch.object(self.m.http.client,'HTTPSConnection',f),patch.object(self.m.time,'sleep') as sleep:
            with self.assertRaises(AssertionError):c.read(3)
        self.assertEqual(s.requests,1);sleep.assert_not_called()


if __name__ == '__main__':unittest.main()
