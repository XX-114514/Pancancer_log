import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import urllib.request

P=Path(__file__).resolve().parents[1]/'scripts/qualify_CM02_RNA_programme_direct_sources_v3.py'
s=importlib.util.spec_from_file_location('direct_request_scope',P);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


class DirectScope(unittest.TestCase):
    def test_direct_route_and_restoration_after_success_or_error(self):
        original=urllib.request.urlopen;calls=[]
        def direct(url,*a,**kw):calls.append(url);return 'synthetic_result'
        def build(handler):
            self.assertEqual(handler.proxies,{})
            return SimpleNamespace(open=direct)
        with patch.object(urllib.request,'build_opener',side_effect=build):
            self.assertEqual(m.with_direct_requests(lambda:urllib.request.urlopen('https://synthetic.invalid')),'synthetic_result')
            self.assertIs(urllib.request.urlopen,original)
            def fail():
                self.assertIs(urllib.request.urlopen,direct);raise RuntimeError('synthetic')
            with self.assertRaises(RuntimeError):m.with_direct_requests(fail)
        self.assertIs(urllib.request.urlopen,original)
        self.assertEqual(calls,['https://synthetic.invalid'])


if __name__=='__main__':unittest.main()
