#!/usr/bin/env python3
from __future__ import annotations
import importlib.util, os, sys
from pathlib import Path
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[variable] = "1"
os.environ["LOKY_MAX_CPU_COUNT"] = "8"
import common
cfg = common.load_config()
base = Path(cfg["legacy_engine_run"]) / "scripts" / "02_run_copykat_python_per_sample.py"
spec = importlib.util.spec_from_file_location("validated_relaxed_copykat", base)
engine = importlib.util.module_from_spec(spec); sys.modules[spec.name] = engine; spec.loader.exec_module(engine)
engine.MANIFEST_ROOT = common.RUN_ROOT / "attempts" / "manifest_v7_relaxed_candidate" / "manifests"
original_write = common.write_json_atomic
def write_versioned(path, payload):
    payload = dict(payload)
    if "/status/copykat/" in str(Path(path).resolve()):
        payload.update({"annotation_parent": "uncertainty_resolution_v2", "role_definition": "relaxed_candidate_strict_reference", "malignancy_evidence_version": "v7_cnv_rerun", "engine_identity_note": "Pinned Copykat_python; not R CopyKAT"})
    original_write(path, payload)
engine.write_json_atomic = write_versioned
engine.main()

