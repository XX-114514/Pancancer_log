"""Single frozen own step closes RAM on success and on a rejected code hash."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]


class SinglePipeline(unittest.TestCase):
    def run_fixture(self,bad_hash=False):
        path=ROOT/'scripts/qualify_independent_melanoma_processed_tar_v6.py'
        spec=importlib.util.spec_from_file_location('single_pipeline_fixture',path)
        m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
        closed=[]
        class Archive:
            def close(self):closed.append(True)
        class Range(io.BytesIO):
            def __init__(self,*args,**kwargs):super().__init__(b'synthetic_input')
        original_load=m.load
        def loader(path,expected,name):
            if path.name=='melanoma_IMC_verified_dns_transport_v2.py':return SimpleNamespace(RangeFile=Range)
            return original_load(path,expected,name)
        ordinary_sha=m.sha
        def digest(raw):
            if raw==b'synthetic_input':return 'd11cd7a5a1672ff002e40658f8d22bea0898733386d24197d817de9adca56112'
            return ordinary_sha(raw)
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);(root/'scripts').mkdir();out=root/'private_output'
            code=b'def analyze(bundle, parameters):\n    return {"status":"SYNTHETIC_COMPLETE"}\n'
            (root/'scripts/own.py').write_bytes(code)
            request=root/'request.json';request.write_text(json.dumps(dict(step='synthetic',code_file='own.py',
                code_SHA256='0'*64 if bad_hash else hashlib.sha256(code).hexdigest(),parameters={})))
            args=['qualify','--output',str(out),'--scripts-root',str(ROOT/'scripts'),'--pipeline-request',str(request)]
            with patch.object(m,'load',side_effect=loader),patch.object(m,'sha',side_effect=digest),patch.object(m,'SIZE',len(b'synthetic_input')),patch.object(m,'qualify',return_value=({'status':'SYNTHETIC_INPUT_ACCEPTED'},Archive())),patch('sys.argv',args):
                if bad_hash:
                    with self.assertRaises(SystemExit):m.main()
                else:m.main()
            status=json.loads((out/'private_RAM_session_status.json').read_text())
            if bad_hash:
                failure=json.loads((out/'private_failed_terminal.json').read_text())
                self.assertEqual(failure['status'],'OWN_PIPELINE_NOT_ACCEPTED_ARCHIVE_INPUT_REMAINS_SEPARATE')
                self.assertEqual(status['status'],'RAM_SESSION_CLOSED_AFTER_FAILURE')
                self.assertFalse((out/'private_single_pipeline_source_summary.json').exists())
            else:
                self.assertEqual(status['status'],'RAM_SESSION_CLOSED_AFTER_SINGLE_PIPELINE')
                self.assertEqual(json.loads((out/'private_single_pipeline_source_summary.json').read_text())['status'],'SYNTHETIC_COMPLETE')
        self.assertTrue(closed)
    def test_success_closes_archive_and_writes_only_summary(self):self.run_fixture()
    def test_bad_code_digest_preserves_input_acceptance_and_closes_archive(self):self.run_fixture(True)


if __name__=='__main__':unittest.main()
