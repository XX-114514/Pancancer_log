"""Run original complete synthetic contracts against the current local backend.

Only test loaders are rebound. Production functions are never mocked. AST
parity additionally checks all scientific reductions against the server code.
"""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def definitions(path):
    return {x.name: ast.dump(x, include_attributes=False) for x in ast.parse(path.read_text()).body
        if isinstance(x, ast.FunctionDef)}


def main():
    reader = ROOT/'scripts'/'summarize_melanoma_IMC_patient_coverage_local_v5.py'
    node = ROOT/'scripts'/'qualify_melanoma_IMC_native_graph_key_bridge_local_v4.py'
    original = definitions(ROOT/'scripts'/'summarize_melanoma_IMC_patient_coverage_v2.py')
    current = definitions(reader)
    names = ['explicit_margin_map', 'summarize', 'scalar_replay', 'author_ROI_infix',
        'canonical_ROIs', 'coordinate_array', 'coordinate_summary', 'native_knn_edges',
        'graph_coverage', 'source_summary_only']
    assert all(original[name] == current[name] for name in names)
    original_node = definitions(ROOT/'scripts'/'qualify_melanoma_IMC_native_graph_key_bridge_v1.py')
    current_node = definitions(node)
    assert all(original_node[name] == current_node[name] for name in ['canonical_images', 'bind_nodes', 'edge_prefix_summary'])
    coverage_tests = load(ROOT/'tests'/'test_melanoma_IMC_local_coalesced_coverage.py', 'current_reader_fixtures')
    ordinary_test_loader = coverage_tests.load
    def current_test_loader(path, name):
        if path.name == 'summarize_melanoma_IMC_patient_coverage_local_v2.py':
            path = reader
        return ordinary_test_loader(path, name)
    coverage_tests.load = current_test_loader
    node_tests = load(ROOT/'tests'/'test_melanoma_IMC_local_coalesced_node_bridge.py', 'current_node_fixtures')
    def node_setup(cls):
        cls.m = load(node, 'current_actual_node_under_test')
    node_tests.LocalNodeBridge.setUpClass = classmethod(node_setup)
    loader = unittest.TestLoader()
    suite = unittest.TestSuite([loader.loadTestsFromModule(coverage_tests), loader.loadTestsFromModule(node_tests)])
    # Reuse complete transport fixtures against the actual dependency chain.
    cache_tests = load(ROOT/'tests'/'test_melanoma_IMC_coalesced_ranges.py', 'current_cache_fixtures')
    def cache_setup(cls):
        cls.m = load(ROOT/'scripts'/'melanoma_IMC_coalesced_ranges_v3.py', 'current_actual_cache')
    cache_tests.CoalescedRanges.setUpClass = classmethod(cache_setup)
    dns_tests = load(ROOT/'tests'/'test_melanoma_IMC_verified_dns_transport.py', 'current_dns_fixtures')
    def dns_setup(cls):
        cls.m = load(ROOT/'scripts'/'melanoma_IMC_verified_dns_transport_v2.py', 'current_actual_dns')
    dns_tests.VerifiedAddressSelection.setUpClass = classmethod(dns_setup)
    http_tests = load(ROOT/'tests'/'test_melanoma_IMC_HTTP_status_recovery.py', 'current_HTTP_fixtures')
    for fixtures in [cache_tests,dns_tests,http_tests]:
        suite.addTests(loader.loadTestsFromModule(fixtures))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    assert result.wasSuccessful() and result.testsRun == 36 and not result.skipped
    print(json.dumps(dict(status='ACTUAL_LOCAL_BACKEND_FULL_SYNTHETIC_CONTRACTS_PASSED',
        tests=result.testsRun, skipped=0, scientific_AST_parity_functions=13,
        source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [reader, node]})))


if __name__ == '__main__':
    main()
