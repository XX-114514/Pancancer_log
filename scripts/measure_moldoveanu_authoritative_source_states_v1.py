"""One frozen RAM pipeline: exact source qualification then final-label audit.

Publishes source aggregates only. Old raw label disagreement never silently
overrides the published final table; each declared probability recipe is tested
before any state effect is computed. Clinical patient inference remains absent.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import importlib.util


def load(path,expected,name):
    assert hashlib.sha256(path.read_bytes()).hexdigest()==expected
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def analyze(bundle,parameters):
    assert parameters==dict(scope='fixed_author_materials_no_patient_inference',permutations=2000,seed=2026102102)
    folder=Path(__file__).parent
    qualifier=load(folder/'qualify_moldoveanu_ICI_cells_in_RAM_v1.py',
        'a498faeb15053036e72003585ba82d5c9b44b2f8f9d3c16d392c3fcf8e693125','exact_author_cells')
    source=qualifier.analyze(bundle,dict(coordinate_tolerance=1e-6,minimum_parent_cells=20))
    core=load(folder/'analyze_moldoveanu_fixed_Ki67_source_organization_v2.py',
        '7bd6561955a40115a014d729bc491d7f863cda9b562f8230a8c3d253b67c871f','final_table_state_core')
    result=core.analyze(bundle,parameters)
    return dict(status=result['status'],completed_UTC=datetime.now(timezone.utc).isoformat(),
        source_qualification=source,analysis=result,
        wrapper_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        no_patient_CI_clinical_model_or_direct_communication_claim=True,
        arrays_only_in_RAM=True)
