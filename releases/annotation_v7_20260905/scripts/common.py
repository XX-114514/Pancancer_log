from __future__ import annotations

import hashlib
import json
import os
import re
import tempfile
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import zarr


RUN_ROOT = Path(__file__).resolve().parent.parent
PROJECT_ROOT = RUN_ROOT.parents[3]
CONFIG_PATH = RUN_ROOT / "config" / "run_config.json"


def load_config() -> dict:
    cfg = json.loads(CONFIG_PATH.read_text())
    for key in (
        "counts_zarr", "v2_annotation_pickle", "v5_final_checkpoint",
        "fixed_gene_annotation", "legacy_engine_run",
    ):
        cfg[key] = str((PROJECT_ROOT / cfg[key]).resolve())
    for key in ("repo_dir", "python"):
        cfg["copykat"][key] = str((PROJECT_ROOT / cfg["copykat"][key]).resolve())
    cfg["gene_annotation"] = cfg["fixed_gene_annotation"]
    cfg["refined_annotation_pickle"] = cfg["v2_annotation_pickle"]
    return cfg


def safe_token(value: str) -> str:
    digest = hashlib.sha1(value.encode("utf-8")).hexdigest()[:12]
    prefix = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)[:64]
    return f"{prefix}_{digest}"


def index_hash(index) -> str:
    import pandas as pd
    values = pd.util.hash_pandas_object(index, index=False).to_numpy(dtype="uint64")
    return hashlib.sha256(values.tobytes()).hexdigest()


def npz_unicode(values) -> np.ndarray:
    import pandas as pd
    text = pd.Series(values, copy=False).astype("string").fillna("")
    return np.asarray(text.tolist(), dtype=np.str_)


def write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        os.replace(temp_name, path)
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)


def read_csr_rows(group: zarr.Group, start: int, end: int, n_vars: int) -> sp.csr_matrix:
    indptr = np.asarray(group["indptr"][start:end + 1], dtype=np.int64)
    data_start, data_end = int(indptr[0]), int(indptr[-1])
    return sp.csr_matrix(
        (
            np.asarray(group["data"][data_start:data_end], dtype=np.float32),
            np.asarray(group["indices"][data_start:data_end], dtype=np.int32),
            indptr - data_start,
        ),
        shape=(end - start, n_vars), dtype=np.float32,
    )


def read_csr_positions(group: zarr.Group, positions: np.ndarray, n_vars: int) -> sp.csr_matrix:
    positions = np.asarray(positions, dtype=np.int64)
    if positions.size == 0:
        return sp.csr_matrix((0, n_vars), dtype=np.float32)
    if np.any(np.diff(positions) <= 0):
        raise ValueError("Sample global positions must be strictly increasing")
    splits = np.flatnonzero(np.diff(positions) > 1) + 1
    runs = np.split(positions, splits)
    matrices = [read_csr_rows(group, int(run[0]), int(run[-1]) + 1, n_vars) for run in runs]
    return matrices[0] if len(matrices) == 1 else sp.vstack(matrices, format="csr")


def open_counts(cfg: dict):
    root = zarr.open_group(cfg["counts_zarr"], mode="r")
    group = root["layers"]["counts"] if "counts" in root["layers"] else root["X"]
    return root, group, tuple(group.attrs["shape"])


def terminal_status(status: str) -> bool:
    return status in {
        "success", "skipped_insufficient_cells", "skipped_insufficient_quality_cells",
        "skipped_algorithm_data_limitation", "skipped_no_candidates",
    }

