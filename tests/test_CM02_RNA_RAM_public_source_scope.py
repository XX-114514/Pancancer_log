import hashlib
import importlib.util
import io
from pathlib import Path
import tarfile
import unittest
from unittest.mock import patch
import urllib.request

P=Path(__file__).resolve().parents[1]/'scripts/qualify_CM02_RNA_programme_RAM_public_sources_v4.py'
s=importlib.util.spec_from_file_location('public_source_RAM',P);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)


def payload(files):
    stream=io.BytesIO()
    with tarfile.open(fileobj=stream,mode='w') as tar:
        for name,body in files:
            info=tarfile.TarInfo(name);info.size=len(body);tar.addfile(info,io.BytesIO(body))
    stream.seek(0);return stream


class PublicSourceRAM(unittest.TestCase):
    def test_two_exact_files_in_memory_no_extraction_and_wrong_sha_rejected(self):
        files=[('a.gz',bytes([31,139])+b'synthetic_A'),('b.gz',bytes([31,139])+b'synthetic_B')]
        expected={name:hashlib.sha256(body).hexdigest() for name,body in files}
        with patch.object(m,'SOURCES',expected):
            cache,proof=m.receive_exact_sources(payload(files))
            self.assertEqual(len(cache),2);self.assertEqual(len(proof),2)
            with self.assertRaises(AssertionError):m.receive_exact_sources(payload([files[0],('b.gz',bytes([31,139])+b'wrong')]))
            with self.assertRaises(AssertionError):m.receive_exact_sources(payload(files+[files[0]]))
    def test_undeclared_file_and_symlink_not_accepted(self):
        with self.assertRaises(AssertionError):m.receive_exact_sources(payload([('../unknown.gz',b'unknown')]))
        body=bytes([31,139])+b'synthetic';expected={'a.gz':hashlib.sha256(body).hexdigest(),'b.gz':hashlib.sha256(body).hexdigest()}
        stream=io.BytesIO()
        with tarfile.open(fileobj=stream,mode='w') as archive:
            link=tarfile.TarInfo('a.gz');link.type=tarfile.SYMTYPE;link.linkname='../private';archive.addfile(link)
            regular=tarfile.TarInfo('b.gz');regular.size=len(body);archive.addfile(regular,io.BytesIO(body))
        stream.seek(0)
        with patch.object(m,'SOURCES',expected),self.assertRaises(AssertionError):m.receive_exact_sources(stream)
    def test_request_scope_and_restoration_on_failure(self):
        original=urllib.request.urlopen
        def work():
            with urllib.request.urlopen('https://synthetic.invalid/a') as response:self.assertEqual(response.read(),b'synthetic')
            with self.assertRaises(AssertionError):urllib.request.urlopen('https://undeclared.invalid')
            raise RuntimeError('synthetic')
        with self.assertRaises(RuntimeError):m.with_cached_sources({'https://synthetic.invalid/a':b'synthetic'},work)
        self.assertIs(urllib.request.urlopen,original)


if __name__=='__main__':unittest.main()
