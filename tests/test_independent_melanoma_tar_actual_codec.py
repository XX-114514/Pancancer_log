"""Extension never substitutes compression signatures or per-table eligibility."""
import importlib.util
from pathlib import Path
import io
import gzip
import hashlib
import tarfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('actual_codec',ROOT/'scripts/qualify_independent_melanoma_processed_tar_v4.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def source(payload):
    out=io.BytesIO()
    with tarfile.open(fileobj=out,mode='w') as tar:
        item=tarfile.TarInfo('Data/ICI_meanIntensity.tsv.gz');item.size=len(payload);tar.addfile(item,io.BytesIO(payload))
    raw=gzip.compress(out.getvalue())
    return raw


class ActualCodec(unittest.TestCase):
    def qualify(self,payload):
        raw=source(payload);result,archive=m.qualify(raw,len(raw),hashlib.md5(raw).hexdigest());archive.close();return result
    def test_plain_table_with_gz_suffix_is_explicit_not_assumed_gzip(self):
        result=self.qualify(b'sample_id\tKi67\nSECRET\t1\n')
        table=result['generic_tables'][0]
        self.assertEqual(table['status'],'TABLE_HEADER_ACCEPTED_NOT_ROWS')
        self.assertEqual(table['detected_codec'],'plain_text_candidate')
        self.assertEqual(table['recognized_generic_fields'],['Ki67','sample_id'])
        self.assertNotIn('SECRET',str(result))
    def test_actual_gzip_magic_is_decompressed(self):
        table=self.qualify(gzip.compress(b'sample_id\tKi67\nSECRET\t1\n'))['generic_tables'][0]
        self.assertEqual(table['status'],'TABLE_HEADER_ACCEPTED_NOT_ROWS');self.assertEqual(table['detected_codec'],'gzip')
    def test_binary_or_truncated_gzip_does_not_accept_table(self):
        for body in [b'\xff\xfeBINARY\n',b'\x1f\x8b']:
            table=self.qualify(body)['generic_tables'][0]
            self.assertEqual(table['status'],'TABLE_HEADER_NOT_ACCEPTED')
            self.assertNotIn('recognized_generic_fields',table)


if __name__=='__main__':unittest.main()
