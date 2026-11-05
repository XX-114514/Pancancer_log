#!/usr/bin/env python3
"""v3: same as 05c_rebuild_selected_lineages_v2.py, plus a fix to
`resolution_for()`'s cell-count-based resolution cap, based on a post-hoc
audit (MARKER_COVERAGE_AUDIT_20260825.md) that found the true bottleneck
behind Endothelial/Plasma/Erythroid/ILC/HSPC being 100% "unresolved" at
Stage06: NOT marker-panel sparsity (their SPECS panels are rich -- 3-8
subtypes each -- and gene coverage across the 54 GSEs is excellent, 89-100%
for every canonical marker checked) but Leiden clustering being forced too
coarse relative to how many literature-defined subtypes each lineage needs
to separate:

    lineage      n_cells  old resolution  n_clusters  SPECS non-contam subtypes
    Erythroid      6,192  0.30 (capped)   1           3  (Early/Late/Cycling)
    HSPC           9,509  0.30 (capped)   1           5  (HSC_MPP/GMP/MEP/Lymphoid_progenitor/Cycling)
    ILC           16,259  0.30 (capped)   1           5  (ILC1/ILC2/ILC3/NK_like/Cycling)
    Endothelial   99,779  0.50 (capped)   3           8  (Blood/Lymphatic/Tip/Venous/Arterial/Capillary/IFN/Cycling)

With only 1 Leiden cluster, a lineage's Stage06 subtype call is mechanically
impossible to resolve regardless of marker quality: there is no "cluster vs.
rest of lineage" contrast to score, because every cell in the lineage IS the
one cluster. This is a pure clustering-granularity bug, not a marker-panel or
gene-coverage problem (that hypothesis was checked and rejected).

Fix: raise the resolution floor for the two low-cell-count brackets (the
<5,000 "too few to subcluster at all" bracket is left untouched -- it only
affects Platelet, 1,130 cells, correctly still coarse). This is deliberately
a moderate bump, not jumping straight to each lineage's SPECS base
resolution, because Stage06's existing `min_cluster_cells` floor per tier
will discard any spurious tiny fragment clusters this produces -- there is
already a downstream safety net, so this errs toward more granularity rather
than less:

    bracket          old cap   new cap
    n < 5,000         0.0       0.0    (unchanged)
    5,000-20,000      0.30      0.55
    20,000-100,000    0.50      0.65
    >= 100,000        base      base   (unchanged; Neural at 107,667 cells is
                                        already uncapped at 0.6 and still only
                                        produced 3 clusters -- that lineage's
                                        under-resolution should be re-examined
                                        AFTER the 02_major_evidence_v3.py
                                        Melanocyte-exclusion fix removes the
                                        ~24,620 contaminating melanoma cells,
                                        not by blindly raising its base
                                        resolution here)

Reads major_evidence_assignments_v3.pkl.gz (02_major_evidence_v3.py's Neural/
HSPC exclusion-marker fix) instead of v2's. Otherwise structurally identical
to 05c_rebuild_selected_lineages_v2.py; all outputs go to *_v3 paths so v1/v2
are untouched (append-only contract)."""

from __future__ import annotations

import logging
import os
import time
from pathlib import Path

import anndata as ad
import bbknn
import h5py
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
    if n_cells < 5_000:
        return 0.0
    if n_cells < 20_000:
        return min(base, 0.55)  # v3: was 0.30
    if n_cells < 100_000:
        return min(base, 0.65)  # v3: was 0.50
    return base


def main() -> None:
    started = time.time()
    started_iso = now_iso()
    ensure_layout()
    config = load_config()
    decision = pd.read_csv(ROOT / "tables" / "integration_basis_decision.tsv", sep="\t").iloc[0]
    method = str(decision["selected_method"])
    selection_path = ROOT / "tables" / "selected_lineage_mapping_v3.tsv"
    summary_path = ROOT / "reports" / "05c_selected_lineage_graph_summary_v3.json"
    refuse_existing(selection_path, summary_path)
    provisional = ROOT / "checkpoints" / "lineage_assignments_harmony_provisional.pkl.gz"  # unused on this run (bbknn selected)
    if method == "harmony":
        if not provisional.exists():
            raise FileNotFoundError(provisional)
        pd.DataFrame([{"selected_method": method, "mapping_path": str(provisional),
                       "rerun_performed": False, "reason": "Harmony won formal scib-metrics benchmark"}]).to_csv(
            selection_path, sep="\t", index=False)
        ended_iso = now_iso()
        write_json_new(summary_path, {
            "run_id": config["run_id"], "stage": "05c_selected_lineage_graph", "status": "success",
            "start_time": started_iso, "end_time": ended_iso, "duration_sec": time.time() - started,
            "params": {"selected_method": method}, "inputs": [str(provisional)],
            "outputs": [str(selection_path)], "warnings": [], "error": None,
            "source_scripts": [str(Path(__file__).resolve())], "logic_version": config["logic_version"],
        })
        append_status("05c", "select_lineage_graph", "success", started_iso, ended_iso,
                      round(time.time() - started, 3), str(provisional), str(selection_path), str(os.getpid()), "")
        return
    if method not in {"scvi", "bbknn"}:
        raise ValueError(f"Unsupported selected method: {method}")

    input_path = ROOT / "checkpoints" / "major_evidence_assignments_v3.pkl.gz"
    output = ROOT / "checkpoints" / f"lineage_assignments_selected_{method}_v3.pkl.gz"
    part_dir = ROOT / "checkpoints" / f"lineage_parts_selected_{method}_v3"
    part_dir.mkdir(parents=True, exist_ok=True)
    refuse_existing(output)
    append_status("05cv3", f"rebuild_lineages_selected_{method}_v3", "running", started_iso,
                  input_ref=str(input_path), job_id=str(os.getpid()))
    mapping = pd.read_pickle(input_path)
    if method == "scvi":
        with h5py.File(resolve_project_path(config["scvi_h5ad"]), "r") as handle:
            embedding = np.asarray(handle["obsm"]["X_scVI"][:], dtype=np.float32)
    else:
        root = zarr.open_group(str(resolve_project_path(config["harmony_zarr"])), mode="r")
        embedding = np.asarray(root["obsm"]["X_pca"][:], dtype=np.float32)
    major = mapping["major_celltype_evidence"].astype(str).to_numpy()
    batches = mapping["gse_id"].astype(str).to_numpy()
    parts: list[pd.DataFrame] = []
    rows: list[dict[str, object]] = []
    for name in sorted(np.unique(major), key=lambda x: int((major == x).sum())):
        pos = np.flatnonzero(major == name)
        part_path = part_dir / f"{name}_assignments.pkl.gz"
        if part_path.exists():
            part = pd.read_pickle(part_path)
            if not part.index.equals(mapping.index[pos]):
                raise ValueError(f"Invalid resumable part {part_path}")
            parts.append(part)
            rows.append({"major_celltype": name, "n_cells": len(pos),
                         "resolution": resolution_for(len(pos), BASE_RESOLUTION.get(name, 0.5)),
                         "n_clusters": int(part["cluster_uid"].nunique()),
                         "status": "resumed_existing_checkpoint", "knn_backend": "not_applicable",
                         "min_batch_cells": int(pd.Series(batches[pos]).value_counts().min())})
            continue
        n_cells = len(pos)
        resolution = resolution_for(n_cells, BASE_RESOLUTION.get(name, 0.5))
        if name in {"Unknown", "Ambiguous"}:
            labels, status = np.repeat("unresolved", n_cells), "not_clustered_evidence_unresolved"
        elif name == "Epithelial":
            labels, status = np.repeat("coarse", n_cells), "deferred_cancer_and_malignancy_stratification"
        elif resolution == 0:
            labels, status = np.repeat("coarse", n_cells), "coarse_rare_lineage"
        else:
            tmp = ad.AnnData(X=np.zeros((n_cells, 1), dtype=np.float32),
                             obs=pd.DataFrame(index=mapping.index[pos].copy()))
            if method == "scvi":
                tmp.obsm["X_selected"] = np.asarray(embedding[pos], dtype=np.float32, order="C")
                sc.pp.neighbors(tmp, n_neighbors=int(config["lineage_neighbors"]), use_rep="X_selected",
                                transformer="pynndescent", random_state=int(config["seed"]), key_added="selected")
                adjacency = tmp.obsp["selected_connectivities"]
            else:
                local_batches = batches[pos]
                min_batch = int(pd.Series(local_batches).value_counts().min())
                n_within = 2 if min_batch >= 2 else 1
                # PyNNDescent in bbknn hard-fails when any batch has <11 cells.
                # Annoy has no such minimum and lets us preserve the real GSE labels
                # instead of pooling rare datasets into an artificial batch.
                knn_backend = "pynndescent" if min_batch >= 11 else "annoy"
                _, adjacency, _ = bbknn.matrix.bbknn(
                    np.asarray(embedding[pos], dtype=np.float32), local_batches,
                    neighbors_within_batch=n_within, n_pcs=embedding.shape[1], trim=200,
                    computation=knn_backend, annoy_n_trees=10, pynndescent_n_neighbors=30,
                    pynndescent_random_state=int(config["seed"]), metric="euclidean")
            sc.tl.leiden(tmp, adjacency=adjacency, resolution=resolution, random_state=int(config["seed"]),
                         key_added="selected_leiden", directed=False, use_weights=True,
                         flavor="igraph", n_iterations=2)
            labels = tmp.obs["selected_leiden"].astype(str).to_numpy()
            status = f"success_recomputed_within_lineage_{method}"
            del tmp, adjacency
        part = pd.DataFrame({
            "global_position": pos.astype(np.int64), "major_celltype_evidence": name,
            "lineage_cluster_raw": labels,
            "cluster_uid": [f"{name}_{label}" for label in labels], "integration_basis": method,
        }, index=mapping.index[pos].copy())
        part.index.name = "cell_id"
        part.to_pickle(part_path, compression="gzip")
        parts.append(part)
        rows.append({"major_celltype": name, "n_cells": n_cells, "resolution": resolution,
                     "n_clusters": int(part["cluster_uid"].nunique()), "status": status,
                     "knn_backend": knn_backend if method == "bbknn" and resolution > 0
                     and name not in {"Unknown", "Ambiguous", "Epithelial"} else "not_applicable",
                     "min_batch_cells": int(pd.Series(batches[pos]).value_counts().min())})
        print(rows[-1], flush=True)
    combined = pd.concat(parts).sort_values("global_position")
    if not combined.index.equals(mapping.index):
        raise ValueError("Selected-basis lineage mapping order mismatch")
    combined.to_pickle(output, compression="gzip")
    pd.DataFrame([{"selected_method": method, "mapping_path": str(output),
                   "rerun_performed": True, "reason": "non-Harmony method won formal scib-metrics benchmark"}]).to_csv(
        selection_path, sep="\t", index=False)
    pd.DataFrame(rows).to_csv(ROOT / "tables" / f"selected_{method}_lineage_graph_summary_v3.tsv", sep="\t", index=False)
    ended_iso = now_iso()
    duration = time.time() - started
    write_json_new(summary_path, {
        "run_id": config["run_id"], "stage": "05c_selected_lineage_graph", "status": "success",
        "start_time": started_iso, "end_time": ended_iso, "duration_sec": duration,
        "params": {"selected_method": method}, "inputs": [str(input_path)],
        "outputs": [str(output), str(selection_path)],
        "warnings": ["Epithelial remains deferred to cancer/CNV-aware final annotation",
                     "v3 raises the resolution cap for the 5,000-100,000 cell brackets "
                     "(see script docstring) to fix Erythroid/HSPC/ILC/Endothelial "
                     "under-clustering found by MARKER_COVERAGE_AUDIT_20260825.md"],
        "error": None,
        "source_scripts": [str(Path(__file__).resolve())], "logic_version": config["logic_version"],
    })
    append_status("05cv3", f"rebuild_lineages_selected_{method}_v3", "success", started_iso, ended_iso,
                  round(duration, 3), str(input_path), str(output), str(os.getpid()), "")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logging.exception("selected lineage graph stage failed")
        append_status("05cv3", "selected_lineage_graph_v3", "failed", error_log=repr(exc), job_id=str(os.getpid()))
        raise
