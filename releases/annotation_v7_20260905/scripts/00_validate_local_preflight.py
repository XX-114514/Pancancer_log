#!/usr/bin/env python3
"""Fail closed unless all V7 inputs, engines, mappings and scripts are local."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import zarr

from common import RUN_ROOT, load_config


def main() -> None:
    cfg = load_config()
    required = [
        cfg["counts_zarr"], cfg["v2_annotation_pickle"], cfg["v5_final_checkpoint"],
        cfg["fixed_gene_annotation"], cfg["legacy_engine_run"], cfg["copykat"]["repo_dir"],
        cfg["copykat"]["python"],
    ]
    missing = [path for path in required if not Path(path).exists()]
    if missing:
        raise FileNotFoundError(f"Missing local dependencies: {missing}")

    forbidden = ("requests.", "http://", "https://", "wget ", "curl ")
    remote_hits = []
    for path in sorted((RUN_ROOT / "scripts").glob("*")):
        if path.suffix not in {".py", ".sh"}:
            continue
        text = path.read_text()
        for term in forbidden:
            if term in text:
                remote_hits.append({"file": str(path), "term": term})
    if remote_hits:
        raise RuntimeError(f"Remote dependency tokens detected: {remote_hits}")

    root = zarr.open_group(cfg["counts_zarr"], mode="r")
    group = root["layers"]["counts"] if "counts" in root["layers"] else root["X"]
    shape = tuple(group.attrs["shape"])
    if shape != (cfg["expected_cells"], cfg["expected_genes"]):
        raise ValueError(f"Unexpected counts shape: {shape}")
    data = np.asarray(group["data"][: min(100000, group["data"].shape[0])])
    if data.size and (data.min() < 0 or not np.allclose(data, np.round(data), atol=1e-6)):
        raise ValueError("Counts sample is not non-negative integer-valued")
    genes = ad.io.read_elem(root["var"]).index.astype(str)
    annotation = pd.read_csv(
        cfg["fixed_gene_annotation"], sep=r"\s+", engine="python", header=None,
        names=["gene_symbol", "chromosome", "start", "end"],
    )
    mapped = genes.isin(set(annotation["gene_symbol"].astype(str)))
    if int(mapped.sum()) != len(genes):
        raise ValueError(f"Gene annotation mapping is not exact: {int(mapped.sum())}/{len(genes)}")

    import infercnvpy  # noqa: F401
    import scanpy  # noqa: F401
    copy_probe = subprocess.run(
        [cfg["copykat"]["python"], "-c", "import anndata,pandas,zarr,copykat_py.copykat; print('ok')"],
        capture_output=True, text=True, check=False,
    )
    if copy_probe.returncode != 0 or "ok" not in copy_probe.stdout:
        raise RuntimeError(f"Copykat_python environment probe failed: {copy_probe.stderr}")
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=cfg["copykat"]["repo_dir"], text=True).strip()
    if commit != cfg["copykat"]["expected_commit"]:
        raise RuntimeError(f"Copykat_python commit mismatch: {commit}")

    report = {
        "status": "PASS", "all_dependencies_local": True, "remote_dependency_hits": [],
        "counts_shape": list(shape), "counts_sample_integer_nonnegative": True,
        "gene_annotation_mapping": f"{int(mapped.sum())}/{len(genes)}",
        "infercnvpy_import": "success", "copykat_python_import": "success",
        "copykat_engine": cfg["copykat"]["engine"], "copykat_commit": commit,
    }
    out = RUN_ROOT / "reports" / "00_local_preflight.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

