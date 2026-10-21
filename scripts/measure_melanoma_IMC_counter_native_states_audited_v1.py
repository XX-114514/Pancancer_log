"""Separate fixed native k10 exploration in retained, qualified RAM.

No network calls, full-author-graph equivalence or clinical claims. The complete
counter bridge and a partial geometry check precede all state effect reading.
"""
from pathlib import Path
import hashlib
import importlib.util
import json


def load(path, expected, name):
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def analyze(bundle, parameters):
    folder = Path(__file__).parent
    root = folder.parent/'private_output'
    geometry_path = root/'private_frozen_step_counter_geometry_prefix_v1_source_summary.json'
    geometry = json.loads(geometry_path.read_text())
    assert geometry['status'] == 'COMPLETE_COUNTER_GEOMETRY_PREFIX_QC_NOT_FULL_GEOMETRIC_VALIDATION'
    assert geometry['native20_compatibility_in_prefix']
    assert geometry['ROI_extent_exceeds_test_radius']
    assert geometry['complete_counter_proof']['all_regions_and_types_exact']
    core = load(folder/'analyze_melanoma_IMC_counter_native_states_v3.py',
        'acf9136480a1258d474595adb286123a9c2f4ca6e7046373ea41539eb678608c', 'counter_native_core')
    auditor = load(folder/'melanoma_IMC_numeric_replay_audit_v1.py',
        '0b8d9d82e5b353a16e368b2b49775c4b675afcbd3804e85926a800425c823edb',
        'counter_native_numeric_audit')
    counts = auditor.attach(core)
    result = core.analyze(bundle, parameters)
    result['independent_numeric_replay'] = dict(counts, all_checks_passed=True,
        no_statistical_values_or_random_stream_changed=True)
    result['graph_provenance'] = dict(
        definition='Native coordinate within ROI k10 undirected union including boundary ties',
        distinct_from_unavailable_full_author_graph=True,
        author_graph_transport_failed_before_state_effects=True,
        source_geometry_prefix_SHA256=hashlib.sha256(geometry_path.read_bytes()).hexdigest(),
        physical_micron_calibration_not_established=True,
        no_contact_distance_or_exact_source_software_version_claim=True,
        thresholds_family_and_random_seed_unchanged=True)
    result['wrapper_SHA256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return result
