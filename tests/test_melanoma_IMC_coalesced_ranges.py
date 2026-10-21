"""Independent byte/ZIP/CRC replay for the bounded block cache, no network."""
import importlib.util
import io
from pathlib import Path
import random
import unittest
from unittest.mock import patch
import zipfile


class CoalescedRanges(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[1]/'scripts'/'melanoma_IMC_coalesced_ranges_v1.py'
        spec = importlib.util.spec_from_file_location('coalesced_under_test', path)
        cls.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.m)

    def setup_source(self, blob):
        requests = []
        def fetch(client, count):
            start = client.position
            requests.append((start, count))
            result = blob[start:start+count]
            client.position += len(result)
            client.requested += count
            client.transferred += len(result)
            return result
        rf = self.m.RangeFile(None, 'https://zenodo.org/synthetic-fixture', len(blob))
        return rf, requests, patch.object(self.m._base.RangeFile, '_read_chunk', fetch)

    def test_tiny_footer_and_directory_use_single_tail_range(self):
        blob = random.Random(1).randbytes(200333)
        rf, requests, patched = self.setup_source(blob)
        with patched:
            rf.seek(-22, 2); self.assertEqual(rf.read(), blob[-22:])
            rf.seek(-42, 2); self.assertEqual(rf.read(20), blob[-42:-22])
            rf.seek(-34400, 2); self.assertEqual(rf.read(34315), blob[-34400:-85])
        self.assertEqual(requests, [(len(blob)-65536, 65536)])

    def test_sequential_and_overlap_reads_match_original_bytes(self):
        blob = random.Random(2).randbytes(400333)
        rf, requests, patched = self.setup_source(blob)
        with patched:
            rf.seek(317); self.assertEqual(rf.read(150000), blob[317:150317])
            count = len(requests)
            rf.seek(17000); self.assertEqual(rf.read(100000), blob[17000:117000])
            self.assertEqual(len(requests), count)
            self.assertLessEqual(len(rf.cache), 8)
            self.assertTrue(all(n <= 65536 for _, n in requests))

    def test_clear_cache_forces_fresh_physical_replay(self):
        blob = random.Random(3).randbytes(200333)
        rf, requests, patched = self.setup_source(blob)
        with patched:
            rf.seek(-22, 2); first = rf.read()
            rf.clear_cache()
            rf.seek(-22, 2); self.assertEqual(rf.read(), first)
        self.assertEqual(len(requests), 2)
        self.assertEqual(rf.requested, 2*65536)

    def test_full_synthetic_zip_member_CRC_and_fresh_directory(self):
        member = random.Random(4).randbytes(200333)
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr('data/selected.bin', member)
        blob = buffer.getvalue()
        rf, requests, patched = self.setup_source(blob)
        with patched:
            with zipfile.ZipFile(rf) as archive:
                self.assertEqual(archive.read('data/selected.bin'), member)
                crc = archive.getinfo('data/selected.bin').CRC
            rf.clear_cache()
            with zipfile.ZipFile(rf) as archive:
                self.assertEqual(archive.getinfo('data/selected.bin').CRC, crc)
        self.assertTrue(all(n <= 65536 for _, n in requests))

    def test_no_full_archive_or_pointer_change_on_failure(self):
        blob = b'x'*100000
        rf, _, _ = self.setup_source(blob)
        rf.size = 32*1024**2
        with self.assertRaises(AssertionError):rf.read()
        def fail(client, count):
            client.position = 999
            raise TimeoutError()
        rf.size = len(blob); rf.seek(123)
        with patch.object(self.m._base.RangeFile, '_read_chunk', fail):
            with self.assertRaises(TimeoutError):rf.read(20)
        self.assertEqual(rf.tell(), 123)


if __name__ == '__main__':
    unittest.main()
