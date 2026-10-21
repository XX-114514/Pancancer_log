import hashlib
import importlib.util
import io
from pathlib import Path
import unittest
import urllib.error

P=Path(__file__).resolve().parents[1]/'scripts/qualify_CM02_RNA_programme_cached_sources_v2.py'
s=importlib.util.spec_from_file_location('exact_source_cache',P);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


class Response(io.BytesIO):
    status=200
    def __init__(self,raw,length=None):super().__init__(raw);self.headers={'Content-Length':str(len(raw) if length is None else length)}
    def geturl(self):return 'https://synthetic.invalid/file.gz'


class ExpectedSource(unittest.TestCase):
    def test_full_exact_bytes_accepted_and_truncation_retried_bounded(self):
        raw=bytes([31,139])+b'synthetic';expected=hashlib.sha256(raw).hexdigest();calls=[]
        def opener(*args,**kwargs):
            calls.append(True);return Response(raw,length=len(raw)+1) if len(calls)==1 else Response(raw)
        got,proof=m.acquire(opener,'https://synthetic.invalid/file.gz',expected,sleeper=lambda seconds:None)
        self.assertEqual(got,raw);self.assertEqual(proof['attempts'],2)
        self.assertEqual(proof['transient_failure_classes'],['IncompleteRead'])
    def test_wrong_complete_hash_is_not_retried_or_accepted(self):
        calls=[]
        def opener(*args,**kwargs):calls.append(True);return Response(bytes([31,139])+b'changed')
        with self.assertRaises(AssertionError):m.acquire(opener,'https://synthetic.invalid/file.gz','0'*64,sleeper=lambda s:None)
        self.assertEqual(len(calls),1)
    def test_permission_error_not_retried(self):
        calls=[]
        def opener(*args,**kwargs):
            calls.append(True);raise urllib.error.HTTPError('https://synthetic.invalid',403,'denied',{},None)
        with self.assertRaises(urllib.error.HTTPError):m.acquire(opener,'https://synthetic.invalid/file.gz','0'*64,sleeper=lambda s:None)
        self.assertEqual(len(calls),1)
    def test_absolute_deadline_is_terminal_and_transient_attempts_bounded(self):
        calls=[]
        def deadline(*args,**kwargs):calls.append(True);raise m.SourceDeadlineExceeded()
        with self.assertRaises(m.SourceDeadlineExceeded):m.acquire(deadline,'https://synthetic.invalid/file.gz','0'*64,sleeper=lambda s:None)
        self.assertEqual(len(calls),1)
        calls.clear()
        def failed(*args,**kwargs):calls.append(True);raise urllib.error.URLError('synthetic transient')
        with self.assertRaises(urllib.error.URLError):m.acquire(failed,'https://synthetic.invalid/file.gz','0'*64,sleeper=lambda s:None)
        self.assertEqual(len(calls),3)


if __name__=='__main__':unittest.main()
