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
        path = Path(__file__).resolve().parents[1]/'scripts'/'melanoma_IMC_coalesced_ranges_v4.py'
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
        blob = random.Random(1).randbytes(800333)
        rf, requests, patched = self.setup_source(blob)
        with patched:
            rf.seek(-22, 2); self.assertEqual(rf.read(), blob[-22:])
            rf.seek(-42, 2); self.assertEqual(rf.read(20), blob[-42:-22])
            rf.seek(-34400, 2); self.assertEqual(rf.read(34315), blob[-34400:-85])
        self.assertEqual(requests, [(len(blob)-262144, 262144)])

    def test_sequential_and_overlap_reads_match_original_bytes(self):
        blob = random.Random(2).randbytes(1600333)
        rf, requests, patched = self.setup_source(blob)
        with patched:
            rf.seek(317); self.assertEqual(rf.read(150000), blob[317:150317])
            count = len(requests)
            rf.seek(17000); self.assertEqual(rf.read(100000), blob[17000:117000])
            self.assertEqual(len(requests), count)
            self.assertLessEqual(len(rf.cache), 8)
            self.assertTrue(all(n <= 262144 for _, n in requests))

    def test_clear_cache_forces_fresh_physical_replay(self):
        blob = random.Random(3).randbytes(800333)
        rf, requests, patched = self.setup_source(blob)
        with patched:
            rf.seek(-22, 2); first = rf.read()
            rf.clear_cache()
            rf.seek(-22, 2); self.assertEqual(rf.read(), first)
        self.assertEqual(len(requests), 2)
        self.assertEqual(rf.requested, 2*262144)

    def test_full_synthetic_zip_member_CRC_and_fresh_directory(self):
        member = random.Random(4).randbytes(800333)
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
        self.assertTrue(all(n <= 262144 for _, n in requests))

    def test_no_full_archive_or_pointer_change_on_failure(self):
        blob = b'x'*100000
        rf, _, _ = self.setup_source(blob)
        rf.size = 32*1024**2
        with self.assertRaises(AssertionError):rf.read()
        def fail(client, count):
            client.position = 999
            raise TimeoutError()
        rf.size = len(blob); rf.seek(123)
        with patch.object(self.m._base.RangeFile, '_read_chunk', fail), patch.object(self.m.time, 'sleep'):
            with self.assertRaises(TimeoutError):rf.read(20)
        self.assertEqual(rf.tell(), 123)

    def test_transient_group_failure_keeps_logical_progress_and_recovers(self):
        blob=b'x'*800333
        rf,requests,patched=self.setup_source(blob)
        original_fetch=patched.new
        calls=[]
        def transient(client,count):
            calls.append(count)
            if len(calls)==1:raise TimeoutError()
            return original_fetch(client,count)
        rf.seek(17)
        with patch.object(self.m._base.RangeFile,'_read_chunk',transient),patch.object(self.m.time,'sleep') as sleep:
            self.assertEqual(rf.read(2),b'xx')
        self.assertEqual(rf.tell(),19);sleep.assert_called_once_with(30)

    def test_exhausted_groups_restore_entire_logical_read_position(self):
        blob=b'x'*800333
        rf,requests,patched=self.setup_source(blob)
        calls=[]
        def fail(client,count):calls.append(count);raise TimeoutError()
        rf.seek(17)
        with patch.object(self.m._base.RangeFile,'_read_chunk',fail),patch.object(self.m.time,'sleep') as sleep:
            with self.assertRaises(TimeoutError):rf.read(2)
        self.assertEqual(rf.tell(),17);self.assertEqual(len(calls),3)
        self.assertEqual([call.args[0] for call in sleep.call_args_list],[30,60])

    def test_TLS_and_integrity_errors_are_not_outer_retried(self):
        import ssl
        for error in [ssl.SSLCertVerificationError(1,'synthetic'),AssertionError('identity_change')]:
            rf,_,_=self.setup_source(b'x'*800333)
            with patch.object(self.m._base.RangeFile,'_read_chunk',side_effect=error) as fetch,patch.object(self.m.time,'sleep') as sleep:
                with self.assertRaises(type(error)):rf.read(2)
            self.assertEqual(fetch.call_count,1);sleep.assert_not_called()

    def test_real_inner_attempt_limit_composes_to_nine_and_permissions_stop(self):
        path=Path(__file__).with_name('test_melanoma_IMC_HTTP_status_recovery.py')
        spec=importlib.util.spec_from_file_location('status_fixture',path)
        fixtures=importlib.util.module_from_spec(spec);spec.loader.exec_module(fixtures)
        fixture=fixtures.HTTPStatusRecovery();fixture.m=self.m._base._base
        for code,retry_after,expected_attempts in [(503,None,9),(403,None,1),(429,'31',1)]:
            _,state,factory=fixture.run_fixture([code],retry_after=retry_after)
            rf=self.m.RangeFile(None,'https://zenodo.org/synthetic',800333)
            with patch.object(fixture.m.http.client,'HTTPSConnection',factory),patch.object(self.m.time,'sleep') as sleep:
                with self.assertRaises(fixture.m.HTTPRangeStatusError):rf.read(1)
            self.assertEqual(state.requests,expected_attempts)
            self.assertEqual(rf.tell(),0)
            self.assertEqual(rf.requested,expected_attempts*262144)
            if code!=503:sleep.assert_not_called()


if __name__ == '__main__':
    unittest.main()
