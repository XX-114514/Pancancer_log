"""RAM-only source qualification: integrity, unsafe members and header privacy."""
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('processed_scope',ROOT/'scripts/qualify_independent_melanoma_processed_tar_v1.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def archive(name='source/Data/ST4_cell_data.txt',link=False):
    stream=io.BytesIO()
    with tarfile.open(fileobj=stream,mode='w') as tar:
        body=b'sample_id\tcell_id\tKi67\tCD4\nSECRET_RECORD\tSECRET_CELL\t2\t3\n'
        entry=tarfile.TarInfo(name);entry.size=len(body)
        if link:entry.type=tarfile.SYMTYPE;entry.linkname='/tmp/target';entry.size=0
        tar.addfile(entry,io.BytesIO(body) if not link else None)
        code=b'stop("AUTHOR_CODE_MUST_NOT_RUN")\n'
        entry=tarfile.TarInfo('source/Fig4_cell_proportion.R');entry.size=len(code);tar.addfile(entry,io.BytesIO(code))
    return gzip.compress(stream.getvalue())


class ProcessedScope(unittest.TestCase):
    def qualify(self,raw):return m.qualify(raw,len(raw),hashlib.md5(raw).hexdigest())
    def test_complete_integrity_and_header_only(self):
        result,tar=self.qualify(archive());tar.close()
        self.assertTrue(result['full_gzip_CRC_checked'])
        self.assertEqual(result['generic_tables'][0]['recognized_generic_fields'],['CD4','Ki67','cell_id','sample_id'])
        encoded=json.dumps(result)
        self.assertNotIn('SECRET_RECORD',encoded);self.assertNotIn('SECRET_CELL',encoded)
        self.assertNotIn('AUTHOR_CODE_MUST_NOT_RUN',encoded)
    def test_reject_absolute_or_parent_paths(self):
        for name in ['/tmp/ST4_cell_data.txt','source/../ST4_cell_data.txt']:
            with self.assertRaises(AssertionError):self.qualify(archive(name))
    def test_reject_symlink_without_extracting(self):
        with self.assertRaises(AssertionError):self.qualify(archive(link=True))
    def test_gzip_crc_corruption_rejected_even_with_matching_outer_digest(self):
        raw=bytearray(archive());raw[-8]^=1
        with self.assertRaises(gzip.BadGzipFile):self.qualify(bytes(raw))
    def test_outer_digest_and_record_array_guards(self):
        raw=archive()
        with self.assertRaises(AssertionError):m.qualify(raw,len(raw),'0'*32)
        with self.assertRaises(AssertionError):m.source_summary_only({'patient_id':'SECRET'})
        with self.assertRaises(AssertionError):m.source_summary_only(list(range(6)))


if __name__=='__main__':unittest.main()
