"""Apply complete original node/label/ambiguity fixtures to local transport."""
import ast
import importlib.util
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('node_fixtures',
    ROOT/'tests'/'test_melanoma_IMC_native_graph_key_bridge.py')
fixtures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fixtures)


class LocalNodeBridge(fixtures.AuthorGraphIdentity):
    @classmethod
    def setUpClass(cls):
        path = ROOT/'scripts'/'qualify_melanoma_IMC_native_graph_key_bridge_local_v2.py'
        spec = importlib.util.spec_from_file_location('local_node_bridge_under_test', path)
        cls.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.m)


class UnchangedNodeSemantics(unittest.TestCase):
    def test_complete_identity_and_geometry_contracts_unchanged(self):
        def definitions(path):
            return {x.name:ast.dump(x,include_attributes=False) for x in ast.parse(path.read_text()).body
                if isinstance(x,ast.FunctionDef)}
        first = definitions(ROOT/'scripts'/'qualify_melanoma_IMC_native_graph_key_bridge_v1.py')
        second = definitions(ROOT/'scripts'/'qualify_melanoma_IMC_native_graph_key_bridge_local_v2.py')
        for name in ['canonical_images', 'bind_nodes', 'edge_prefix_summary']:
            self.assertEqual(first[name], second[name])


if __name__ == '__main__':
    unittest.main()
