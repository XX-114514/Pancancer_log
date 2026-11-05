#!/usr/bin/env python3
"""Rebuild kNN and Leiden independently inside each evidence-supported lineage.

Harmony is used provisionally to satisfy the requested sequential audit. The
selected integration basis is decided later by the stratified scIB benchmark;
if another basis wins, this stage is rerun to a new method-specific output.
"""

from __future__ import annotations

import logging
import os
import socket
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import zarr

from common import ROOT, append_status, ensure_layout, load_config, now_iso, refuse_existing, resolve_project_path, write_json_new


BASE_RESOLUTION = {
    "Tcell": 1.0, "NK": 0.8, "ILC": 0.6, "Bcell": 0.8, "Plasma": 0.6,
    "Monocyte": 0.8, "Macrophage": 0.9, "Dendritic": 0.8, "Granulocyte": 0.6,
    "Mast": 0.6, "Endothelial": 0.8, "Fibroblast": 0.9, "Erythroid": 0.5,
    "Platelet": 0.4, "HSPC": 0.7, "Neural": 0.6,
}


def resolution_for(n_cells: int, base: float) -> float:
    # Rare lineages need coarser partitions; lowering thresholds would create
    # attractive but irreproducible single-digit subtypes.
    if n_cells < 5_000:
        return 0.0
    if n_cells < 20_000:
        return min(base, 0.30)
    if n_cells < 100_000:
        return min(base, 0.50)
    return base


def setup_logger(path: Path) -> logging.Logger:
    logger = logging.getLogger("lineage_neighbors_leiden")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    for handler in (logging.StreamHandler(), logging.FileHandler(path)):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def main() -> None:
    started = time.time()
    started_iso = now_iso()
    ensure_layout()
    config = load_config()
    input_path = ROOT / "checkpoints" / "major_evidence_assignments.pkl.gz"
    part_dir = ROOT / "checkpoints" / "lineage_parts_harmony_provisional"
    part_dir.mkdir(parents=True, exist_ok=True)
    output = ROOT / "checkpoints" / "lineage_assignments_harmony_provisional.pkl.gz"
    table_path = ROOT / "tables" / "lineage_graph_summary_harmony_provisional.tsv"
    summary_path = ROOT / "reports" / "03_lineage_neighbors_harmony_provisional_summary.json"
    refuse_existing(output, table_path, summary_path)
    if not input_path.exists():
        raise FileNotFoundError(f"Stage 02 output missing: {input_path}")
    logger = setup_logger(ROOT / "logs" / "03_lineage_neighbors_harmony_provisional.log")
    append_status("03", "lineage_neighbors_leiden_harmony_provisional", "running", started_iso,
                  input_ref=str(input_path), job_id=str(os.getpid()))

    mapping = pd.read_pickle(input_path)
    harmony_root = zarr.open_group(str(resolve_project_path(config["harmony_zarr"])), mode="r")
    embedding = np.asarray(ad.io.read_elem(harmony_root["obsm"]["X_pca_harmony"]), dtype=np.float32)
    if embedding.shape != (len(mapping), 50):
        raise ValueError(f"Unexpected Harmony shape {embedding.shape}; expected {(len(mapping), 50)}")
    major = mapping["major_celltype_evidence"].astype(str).to_numpy()
    rows: list[dict[str, object]] = []
    parts: list[pd.DataFrame] = []

    ordered = sorted(np.unique(major), key=lambda name: int((major == name).sum()))
    for name in ordered:
        pos = np.flatnonzero(major == name)
        n_cells = len(pos)
        part_path = part_dir / f"{name}_assignments.pkl.gz"
        if part_path.exists():
            part = pd.read_pickle(part_path)
            if len(part) != n_cells or not part.index.equals(mapping.index[pos]):
                raise ValueError(f"Existing resumable part failed validation: {part_path}")
            parts.append(part)
            rows.append({"major_celltype": name, "n_cells": n_cells, "status": "resumed_existing_part",
                         "resolution": np.nan, "n_clusters": int(part["cluster_uid"].nunique()),
                         "n_neighbors": int(config["lineage_neighbors"]), "duration_sec": 0.0})
            logger.info("resumed %s cells=%s", name, f"{n_cells:,}")
            continue

        lineage_started = time.time()
        if name in {"Unknown", "Ambiguous"}:
            status = "not_clustered_evidence_unresolved"
            labels = np.repeat("unresolved", n_cells)
            resolution = 0.0
        elif name == "Epithelial":
            # Pan-cancer epithelial clusters are dominated by organ/cancer and
            # malignancy. Final subdivision is deferred to cancer-stratified
            # CNV-aware annotation instead of forcing a global epithelial atlas.
            status = "deferred_cancer_and_malignancy_stratification"
            labels = np.repeat("coarse", n_cells)
            resolution = 0.0
        else:
            base = BASE_RESOLUTION.get(name, 0.5)
            resolution = resolution_for(n_cells, base)
            if resolution == 0.0:
                status = "coarse_rare_lineage"
                labels = np.repeat("coarse", n_cells)
            else:
                logger.info("neighbors start major=%s cells=%s resolution=%.2f", name, f"{n_cells:,}", resolution)
                tmp = ad.AnnData(X=np.zeros((n_cells, 1), dtype=np.float32),
                                 obs=pd.DataFrame(index=mapping.index[pos].copy()))
                tmp.obsm["X_lineage_harmony"] = np.asarray(embedding[pos], dtype=np.float32, order="C")
                sc.pp.neighbors(tmp, n_neighbors=int(config["lineage_neighbors"]), use_rep="X_lineage_harmony",
                                transformer="pynndescent", random_state=int(config["seed"]), key_added="lineage")
                sc.tl.leiden(tmp, adjacency=tmp.obsp["lineage_connectivities"], resolution=resolution,
                             random_state=int(config["seed"]), key_added="lineage_leiden",
                             directed=False, use_weights=True, flavor="igraph", n_iterations=2)
                labels = tmp.obs["lineage_leiden"].astype(str).to_numpy()
                status = "success_recomputed_within_lineage_graph"
                del tmp
        cluster_uid = np.asarray([f"{name}_{label}" for label in labels], dtype=object)
        part = pd.DataFrame({
            "global_position": pos.astype(np.int64), "major_celltype_evidence": name,
            "lineage_cluster_raw": labels, "cluster_uid": cluster_uid,
            "integration_basis": "harmony_provisional",
        }, index=mapping.index[pos].copy())
        part.index.name = "cell_id"
        part.to_pickle(part_path, compression="gzip")
        parts.append(part)
        duration = time.time() - lineage_started
        rows.append({"major_celltype": name, "n_cells": n_cells, "status": status,
                     "resolution": resolution, "n_clusters": int(pd.Series(cluster_uid).nunique()),
                     "n_neighbors": int(config["lineage_neighbors"]), "duration_sec": duration})
        logger.info("lineage complete major=%s clusters=%s elapsed_min=%.1f", name,
                    int(pd.Series(cluster_uid).nunique()), duration / 60)

    combined = pd.concat(parts, axis=0).sort_values("global_position")
    if not combined.index.equals(mapping.index):
        raise ValueError("Combined lineage assignments do not restore global cell order")
    combined.to_pickle(output, compression="gzip")
    summary_table = pd.DataFrame(rows)
    summary_table.to_csv(table_path, sep="\t", index=False)
    ended_iso = now_iso()
    duration = time.time() - started
    payload = {
        "run_id": config["run_id"], "stage": "03_lineage_neighbors_leiden_harmony_provisional", "status": "success",
        "start_time": started_iso, "end_time": ended_iso, "duration_sec": duration,
        "params": {"basis": "X_pca_harmony", "n_neighbors": int(config["lineage_neighbors"]),
                   "rare_lineage_policy": "coarse_below_5000", "epithelial_policy": "defer_cancer_malignancy_stratification"},
        "inputs": [str(input_path), str(resolve_project_path(config["harmony_zarr"]))],
        "outputs": [str(output), str(table_path), str(part_dir)],
        "warnings": ["Harmony is provisional until Stage 05 scIB comparison", "graphs are not persisted to limit disk use"],
        "error": None, "source_scripts": [str(Path(__file__).resolve())], "logic_version": config["logic_version"],
        "hostname": socket.gethostname(), "pid": os.getpid(), "n_cells": int(len(mapping)),
        "n_clusters": int(combined["cluster_uid"].nunique()),
    }
    write_json_new(summary_path, payload)
    append_status("03", "lineage_neighbors_leiden_harmony_provisional", "success", started_iso, ended_iso,
                  round(duration, 3), str(input_path), str(output), str(os.getpid()), "")
    logger.info("stage complete cells=%s clusters=%s elapsed_hr=%.2f", f"{len(mapping):,}",
                payload["n_clusters"], duration / 3600)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logging.getLogger("lineage_neighbors_leiden").exception("stage failed")
        append_status("03", "lineage_neighbors_leiden_harmony_provisional", "failed",
                      error_log=repr(exc), job_id=str(os.getpid()))
        raise
