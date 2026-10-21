"""Counter-gate change preserves original mathematics and full source wiring."""
import ast
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def load(filename, name):
    spec = importlib.util.spec_from_file_location(name, ROOT/'scripts'/filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class CounterNativeStates(unittest.TestCase):
    def test_only_analyze_gate_and_documentation_changed(self):
        old = ast.parse((ROOT/'scripts/analyze_melanoma_IMC_continuous_state_neighbors_v2.py').read_text())
        new = ast.parse((ROOT/'scripts/analyze_melanoma_IMC_counter_native_states_v3.py').read_text())
        def functions(tree):
            return {x.name: ast.dump(x) for x in tree.body if isinstance(x, ast.FunctionDef) and x.name != 'analyze'}
        self.assertEqual(functions(old), functions(new))
        a = next(x for x in old.body if isinstance(x, ast.FunctionDef) and x.name == 'analyze')
        b = next(x for x in new.body if isinstance(x, ast.FunctionDef) and x.name == 'analyze')
        for node in ast.walk(a):
            if isinstance(node, ast.Constant) and node.value == 'private_frozen_step_author_node_native_geometry_v1_source_summary.json':
                node.value = 'private_author_node_counter_qualification.json'
            if isinstance(node, ast.Constant) and node.value == 'COMPLETE_AUTHOR_NODE_BIJECTION_AND_EDGE_PREFIX_NATIVE_COORDINATE_REVIEW':
                node.value = 'COMPLETE_EXPLICIT_COUNTER_BASE_AND_LABEL_BIJECTION'
        self.assertEqual(ast.dump(a), ast.dump(b))

    def test_full_2000_draws_18_slots_with_independent_numeric_replay(self):
        core = load('analyze_melanoma_IMC_counter_native_states_v3.py', 'native_full_fixture')
        counts = load('melanoma_IMC_numeric_replay_audit_v1.py', 'native_full_audit').attach(core)
        rows, mapping = [], {}
        for case in range(5):
            roi = 'synthetic_roi_'+str(case)
            mapping[roi] = 'synthetic_case_'+str(case)
            for typ in ['R', 'C', 'E']:
                for cell in range(20):
                    value = 19-cell if typ == 'C' else cell
                    rows.append(dict(sample_id=roi, celltype20=typ,
                        **{marker: float(value) for marker in set(core.R_MARKERS+core.E_MARKERS)}))
        frame = pd.DataFrame(rows)
        coords = np.column_stack([np.tile(np.arange(60), 5), np.zeros(300)])
        def edges(coordinates, k):
            self.assertEqual((len(coordinates), k), (60, 10))
            return np.array([(i, 40+i) for i in range(20)]+[(20+i, 40+i) for i in range(20)])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'scripts').mkdir(); (root/'private_output').mkdir()
            module = root/'scripts'/'fixture.py'; module.write_text('synthetic')
            (root/'private_output'/'private_author_node_counter_qualification.json').write_text(json.dumps(dict(
                status='COMPLETE_EXPLICIT_COUNTER_BASE_AND_LABEL_BIJECTION',
                independent_CSV_Pandas_node_replay=True, complete_node_member_CRC_checked=True)))
            core.__file__ = str(module)
            result = core.analyze(dict(frame=frame, native_coordinates=coords,
                ROI_to_source_patient=mapping, native_knn_edges=edges), dict(
                regulatory_author_label='R', conventional_CD4_author_label='C', effector_author_label='E'))
        self.assertEqual(result['family_size'], 18)
        main = result['family']['ROI_parent_shuffle']
        self.assertTrue(all(row['effect'] > 0 for row in main[:9]))
        self.assertTrue(all(row['effect'] < 0 for row in main[9:]))
        self.assertTrue(all(row['source_case_units'] == 5 and row['case_bootstrap_95_CI'] is None for row in main))
        self.assertEqual(counts['family_positions'], 36)
        self.assertGreater(counts['covariance_calls'], 40000)
        encoded = json.dumps(result, allow_nan=False)
        self.assertNotIn('synthetic_roi_', encoded)
        self.assertNotIn('synthetic_case_', encoded)


if __name__ == '__main__':
    unittest.main()
