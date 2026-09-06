#!/usr/bin/env python3
"""Run the pinned local Copykat_python engine on V7 manifests."""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

for variable in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMBA_NUM_THREADS",
):
    os.environ[variable] = "1"
os.environ["LOKY_MAX_CPU_COUNT"] = "8"

import common

cfg = common.load_config()
base_script = Path(cfg["legacy_engine_run"]) / "scripts" / "02_run_copykat_python_per_sample.py"
spec = importlib.util.spec_from_file_location("validated_v7_copykat_engine", base_script)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Unable to load {base_script}")
engine = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = engine
spec.loader.exec_module(engine)
engine.MANIFEST_ROOT = common.RUN_ROOT / "attempts" / "manifest_v7" / "manifests"

write_original = common.write_json_atomic


def write_versioned(path: Path, payload: dict) -> None:
    if "/status/copykat/" in str(Path(path).resolve()):
        payload = dict(payload)
        payload["annotation_parent"] = "uncertainty_resolution_v2"
        payload["malignancy_evidence_version"] = "v7_cnv_rerun"
        payload["engine_identity_note"] = "Pinned Copykat_python; not R CopyKAT"
    write_original(path, payload)


engine.write_json_atomic = write_versioned
engine.main()

