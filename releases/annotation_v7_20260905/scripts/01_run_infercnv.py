#!/usr/bin/env python3
"""Run the validated per-sample inferCNVpy engine on V7 manifests."""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

for variable in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMBA_NUM_THREADS", "LOKY_MAX_CPU_COUNT",
):
    os.environ[variable] = "1"

import common

cfg = common.load_config()
base_script = Path(cfg["legacy_engine_run"]) / "scripts" / "01_run_infercnvpy_per_sample.py"
pipeline_dir = common.PROJECT_ROOT.parent / "Project_v1" / "Data" / "scripts" / "pan_cancer_pipeline"
if str(pipeline_dir) not in sys.path:
    sys.path.append(str(pipeline_dir))

spec = importlib.util.spec_from_file_location("validated_v7_infercnv_engine", base_script)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Unable to load {base_script}")
engine = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = engine
spec.loader.exec_module(engine)
engine.MANIFEST_ROOT = common.RUN_ROOT / "attempts" / "manifest_v7" / "manifests"

write_original = common.write_json_atomic


def write_versioned(path: Path, payload: dict) -> None:
    if "/status/infercnv/" in str(Path(path).resolve()):
        payload = dict(payload)
        payload["annotation_parent"] = "uncertainty_resolution_v2"
        payload["malignancy_evidence_version"] = "v7_cnv_rerun"
        payload["gene_annotation_version"] = "gencode_v27_HGNC18790_to_NSG1_v1"
        payload["gene_annotation_exact_mapping"] = "10432/10432"
    write_original(path, payload)


engine.write_json_atomic = write_versioned
import infercnvpy as cnv

infer_original = cnv.tl.infercnv


def infer_single_job(*args, **kwargs):
    kwargs.setdefault("n_jobs", 1)
    return infer_original(*args, **kwargs)


cnv.tl.infercnv = infer_single_job
engine.main()

