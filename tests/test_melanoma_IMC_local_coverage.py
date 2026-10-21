"""Apply existing identity/ROI contracts to the local exact-neighbor reader."""
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


fixtures = load(ROOT/'tests'/'test_melanoma_IMC_patient_coverage.py', 'existing_coverage_fixtures')


class SamePatientSemantics(fixtures.PatientCoverageSemantics):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.module = load(ROOT/'scripts'/'summarize_melanoma_IMC_patient_coverage_local_v1.py', 'local_semantics')


class ExactLocalNativeGeometry(fixtures.AliasAndNativeGeometry):
    __unittest_skip__ = False

    @classmethod
    def setUpClass(cls):
        cls.module = load(ROOT/'scripts'/'summarize_melanoma_IMC_patient_coverage_local_v1.py', 'local_geometry')


class LocalCompleteInputGate(unittest.TestCase):
    def test_ten_scientific_functions_unchanged(self):
        original = ROOT/'scripts'/'summarize_melanoma_IMC_patient_coverage_v2.py'
        local = ROOT/'scripts'/'summarize_melanoma_IMC_patient_coverage_local_v1.py'
        names = ['explicit_margin_map', 'summarize', 'scalar_replay', 'author_ROI_infix',
            'canonical_ROIs', 'coordinate_array', 'coordinate_summary', 'native_knn_edges',
            'graph_coverage', 'source_summary_only']
        def definitions(path):
            return {x.name: ast.dump(x, include_attributes=False) for x in ast.parse(path.read_text()).body
                if isinstance(x, ast.FunctionDef)}
        first, second = definitions(original), definitions(local)
        for name in names:
            self.assertEqual(first[name], second[name])

    def test_partial_or_wrong_hash_stops_before_any_new_output_or_network(self):
        module = load(ROOT/'scripts'/'summarize_melanoma_IMC_patient_coverage_local_v1.py', 'local_input_gate')
        complete = 'COMPLETE_SELECTED_PROCESSED_CELL_MEMBER_MEASUREMENT_QUALIFICATION'
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            accepted = root/'accepted.json'
            terminal = root/'private_terminal.json'
            accepted.write_text(json.dumps(dict(status='PARTIAL')))
            terminal.write_text(json.dumps(dict(status='PARTIAL')))
            with patch.dict(os.environ, {'PANCANCER_IMC_ACCEPTED_MEASUREMENT': str(accepted)}, clear=True):
                with self.assertRaises(AssertionError):module.main()
            accepted.write_text(json.dumps(dict(status=complete)))
            terminal.write_text(json.dumps(dict(status=complete,
                producer_SHA256='eb142bcdeb172f068f9d182a6ba89866d3d07e852c75f0b976ec4c5b43cd436d',
                files={accepted.name: dict(SHA256='incorrect')})))
            with patch.dict(os.environ, {'PANCANCER_IMC_ACCEPTED_MEASUREMENT': str(accepted)}, clear=True):
                with self.assertRaises(AssertionError):module.main()
            self.assertEqual(len(list(root.iterdir())), 2)


if __name__ == '__main__':
    unittest.main()
