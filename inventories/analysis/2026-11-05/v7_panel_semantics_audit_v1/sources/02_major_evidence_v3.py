#!/usr/bin/env python3
"""v3: same as 02_major_evidence_v2.py, plus two developmental-lineage-collision
exclusion fixes found by a post-hoc marker/cluster-coverage audit on 2026-08-25:

1. Neural vs Melanocyte/Melanoma (confirmed at 97.9% concentration):
   `Neural_Melanocytic_PMEL` (a subtype *inside* Neural's own SPECS panel,
   pancancer_evidence_marker_specs.py) turned out to be almost entirely Uveal
   melanoma tumor cells (24,102/24,620 cells), not neuroglial tissue -- SOX10
   and S100B are shared between genuine Schwann/neural-crest-derived glia and
   melanocytes (both are neural-crest lineages), so Neural's MAJOR_MARKERS
   panel (S100B, SOX10, PLP1, MPZ, SLC1A3) was never guarded against strong
   melanocyte-specific signal. Adding PMEL/MLANA/TYR/DCT (all >=94% GSE
   coverage per wide_marker_gse_gene_availability.tsv, so this is not a
   coverage-limited fix) as Neural's exclusion markers stops melanocytic
   cells from being confidently called "Neural" -- they fall back to
   Ambiguous/Unknown/original-label evidence like any other unresolved case,
   which is a more honest state than a wrong confident label. This does NOT
   create a new major class or touch malignancy determination (explicitly
   out of scope for this pass); it only stops a specific mislabeling.

2. HSPC vs CD34+ Endothelial/Fibroblast (lower-confidence, precautionary):
   CD34 is shared between hematopoietic stem/progenitor cells and vascular
   endothelium (esp. tip/proliferating cells) and some fibroblast/stromal
   progenitor populations -- a well-documented marker collision. HSPC's
   MAJOR_EXCLUDE_MARKERS previously only excluded CD3D/CD3E/MS4A1 (lymphocyte
   contamination), with no stromal/vascular guard at all. ~39% of the 9,509
   HSPC-labeled cells come from solid-tumor cancer_types with no expected
   substantial bone-marrow infiltration (vs. the clearly-genuine 60.6% from
   Acute myeloid leukemia), consistent with -- but not proof of -- some
   contamination from CD34+ non-hematopoietic cells. Adding PECAM1/VWF
   (Endothelial) and COL1A1 (Fibroblast) as HSPC exclusion markers is a
   low-risk precaution in the same spirit as Tcell's existing PECAM1/COL1A1
   exclusion.

Both additions only ever *remove* a false-positive major-class call (routing
those cells to Ambiguous/Unknown/other-evidence-source instead); they cannot
manufacture a new positive call. Writes to *_v3 output paths; v1/v2 outputs
are untouched (append-only contract). See MARKER_COVERAGE_AUDIT_20260825.md
for the full audit this fix is based on (marker-panel sparsity was RULED OUT
as the cause of the Endothelial/Plasma/Erythroid/ILC/HSPC unresolved problem;
Leiden clustering resolution being capped too low for those lineages'
cell counts is the actual cause -- see 05c_rebuild_selected_lineages_v3.py,
a separate fix, not this one).
"""

from __future__ import annotations

import os
import socket
import time
import logging
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr

from common import ROOT, append_status, ensure_layout, load_config, now_iso, refuse_existing, resolve_project_path, write_json_new


# Conservative broad-lineage panels. Highly state-dependent genes are avoided
# where possible; a single marker is never sufficient for a confident call.
MAJOR_MARKERS: dict[str, tuple[str, ...]] = {
    "Tcell": ("CD3D", "CD3E", "TRAC", "LCK", "CD247"),
    "NK": ("NKG7", "GNLY", "KLRD1", "KLRF1", "FCGR3A"),
    "Bcell": ("CD79A", "MS4A1", "CD37", "CD74", "CD22"),
    "Plasma": ("MZB1", "JCHAIN", "SDC1", "TNFRSF17", "DERL3"),
    "Monocyte": ("S100A8", "S100A9", "FCN1", "LILRB1", "CTSD"),
    "Macrophage": ("C1QA", "C1QB", "C1QC", "APOC1", "MARCO"),
    "Dendritic": ("FCER1A", "CD1C", "CLEC10A", "CLEC9A", "GZMB"),
    "Mast": ("KIT", "TPSAB1", "TPSB2", "CPA3", "MS4A2"),
    "Granulocyte": ("FCGR3B", "CSF3R", "S100A8", "S100A9", "MNDA"),
    "Epithelial": ("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"),
    "Endothelial": ("PECAM1", "VWF", "KDR", "EMCN", "RAMP2"),
    "Fibroblast": ("DCN", "LUM", "COL1A1", "COL1A2", "COL3A1"),
    "Erythroid": ("HBB", "HBA1", "HBA2", "AHSP", "GYPA"),
    "Platelet": ("PPBP", "PF4", "GP9", "NRGN", "ITGA2B"),
    "HSPC": ("CD34", "SPINK2", "GATA2", "HLF", "MEIS1"),
    "Neural": ("S100B", "SOX10", "PLP1", "MPZ", "SLC1A3"),
    "ILC": ("IL7R", "KIT", "KLRB1", "RORC", "IL1R1"),
}

# v3 adds two entries (Neural, HSPC) vs v2; everything else identical.
MAJOR_EXCLUDE_MARKERS: dict[str, tuple[str, ...]] = {
    "Tcell": ("LST1", "TYROBP", "MS4A1", "COL1A1", "PECAM1"),
    "NK": ("CD3D", "CD3E", "MS4A1"),
    "Bcell": ("CD3D", "CD3E", "LST1", "COL1A1"),
    "Plasma": ("CD3D", "CD3E", "MS4A1"),
    "Monocyte": ("CD3D", "CD3E", "MS4A1"),
    "Macrophage": ("CD3D", "CD3E", "MS4A1", "EPCAM"),
    "Dendritic": ("CD3D", "CD3E", "MS4A1"),
    "Mast": ("CD3D", "CD3E", "MS4A1"),
    "Granulocyte": ("CD3D", "CD3E", "MS4A1"),
    "Epithelial": ("PTPRC", "COL1A1", "PECAM1"),
    "Endothelial": ("PTPRC", "COL1A1"),
    "Fibroblast": ("PTPRC", "EPCAM", "PECAM1"),
    "Erythroid": ("PTPRC",),
    "Platelet": ("PTPRC",),
    "HSPC": ("CD3D", "CD3E", "MS4A1", "PECAM1", "VWF", "COL1A1"),  # v3: + vascular/stromal guard
    "Neural": ("PTPRC", "EPCAM", "PMEL", "MLANA", "TYR", "DCT"),  # v3: + melanocyte guard
    "ILC": ("CD3D", "CD3E", "MS4A1"),
}


def setup_logger(path: Path) -> logging.Logger:
    logger = logging.getLogger("major_evidence_v3")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    for handler in (logging.StreamHandler(), logging.FileHandler(path)):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def source_gene_symbols(path: Path) -> set[str]:
    with h5py.File(path, "r") as handle:
        var = ad.io.read_elem(handle["var"])
    names = pd.Index(var.index.astype(str))
    looks_ensembl = names.str.match(r"^ENSG\d+(\.\d+)?$").mean() > 0.5
    symbol_col = next((c for c in ("gene_symbol", "gene_symbols", "gene_name", "symbol", "SYMBOL") if c in var.columns), None)
    if looks_ensembl and symbol_col is not None:
        names = pd.Index(var[symbol_col].astype(str))
    return set(names.str.strip().str.upper())


def normalize_marker_counts(counts: sp.csr_matrix, full_totals: np.ndarray, target_sum: float = 1e4) -> sp.csr_matrix:
    matrix = counts.astype(np.float32, copy=True).tocsr()
    scales = np.divide(target_sum, full_totals, out=np.zeros_like(full_totals, dtype=np.float32), where=full_totals > 0)
    matrix.data *= np.repeat(scales, np.diff(matrix.indptr))
    np.log1p(matrix.data, out=matrix.data)
    matrix.eliminate_zeros()
    return matrix


def main() -> None:
    started = time.time()
    started_iso = now_iso()
    ensure_layout()
    config = load_config()
    logger = setup_logger(ROOT / "logs" / "02_major_evidence_v3.log")
    originals_path = ROOT / "checkpoints" / "per_cell_original_annotations.pkl.gz"
    output = ROOT / "checkpoints" / "major_evidence_assignments_v3.pkl.gz"
    counts_path = ROOT / "tables" / "major_evidence_counts_v3.tsv"
    availability_path = ROOT / "tables" / "major_marker_availability_by_gse_v3.tsv"
    l1_path = ROOT / "tables" / "l1_major_marker_effect_sizes_v3.tsv"
    gse_path = ROOT / "tables" / "gse_major_marker_support_v3.tsv"
    exclude_path = ROOT / "tables" / "exclusion_marker_flips_v3.tsv"
    summary_path = ROOT / "reports" / "02_major_evidence_v3_summary.json"
    refuse_existing(output, counts_path, availability_path, l1_path, gse_path, exclude_path, summary_path)
    if not originals_path.exists():
        raise FileNotFoundError(f"Stage 01 output missing: {originals_path}")
    append_status("02v3", "major_evidence_v3", "running", started_iso,
                  input_ref=str(originals_path), job_id=str(os.getpid()))

    logger.info("reading recovered originals")
    originals = pd.read_pickle(originals_path)
    logger.info("reading legacy mapping for global L1 and CellTypist support only")
    legacy = pd.read_pickle(resolve_project_path(config["legacy_lineage_mapping"]))
    if not originals.index.equals(legacy.index):
        raise ValueError("Recovered originals and legacy mapping indices differ")
    evidence = originals.copy()
    evidence["l1_cluster"] = legacy["l1_cluster"].astype(str).to_numpy()
    evidence["celltypist_major"] = legacy["celltypist_major"].astype(str).to_numpy()
    evidence["celltypist_confidence"] = pd.to_numeric(legacy["celltypist_confidence"], errors="coerce").fillna(0).astype(np.float32).to_numpy()

    logger.info("reading merge obs and full-library totals")
    merge_root = zarr.open_group(str(resolve_project_path(config["merge_zarr"])), mode="r")
    merge_obs = ad.io.read_elem(merge_root["obs"])
    if not evidence.index.equals(merge_obs.index):
        raise ValueError("Evidence and merge obs indices differ")
    totals = pd.to_numeric(merge_obs["total_counts"], errors="coerce").fillna(0).to_numpy(dtype=np.float32)

    logger.info("reading full retry marker counts")
    wide = ad.read_h5ad(resolve_project_path(config["legacy_wide_panel_full_retry"]))
    wide_pos = wide.obs_names.get_indexer(evidence.index)
    if (wide_pos < 0).any() or len(np.unique(wide_pos)) != len(wide_pos):
        raise ValueError("Full retry marker panel is not a one-to-one cell-ID match")
    if not evidence.index.equals(wide.obs_names):
        logger.warning("retry marker rows differ from merge order; applying explicit cell-ID permutation")
    if "counts" not in wide.layers:
        raise KeyError("Full retry marker panel lacks layers['counts']")
    counts = wide.layers["counts"][wide_pos, :]
    counts = counts.tocsr() if sp.issparse(counts) else sp.csr_matrix(counts)
    logx = normalize_marker_counts(counts, totals)
    logger.info("marker counts normalized with full-library totals; nnz=%s", f"{logx.nnz:,}")
    genes = pd.Index(wide.var_names.astype(str).str.upper())
    gene_to_col = {gene: idx for idx, gene in enumerate(genes)}
    classes = list(MAJOR_MARKERS)
    scores = np.full((len(evidence), len(classes)), np.nan, dtype=np.float32)
    exclude_scores = np.full((len(evidence), len(classes)), np.nan, dtype=np.float32)

    availability_rows: list[dict[str, object]] = []
    gse_rows: list[dict[str, object]] = []
    source_by_gse = (merge_obs.assign(_gse=merge_obs["gse_id"].astype(str))
                     .groupby("_gse", observed=True)["phase07_source_h5ad"].first())
    gse_values = evidence["gse_id"].astype(str).to_numpy()
    for gse_id in sorted(np.unique(gse_values)):
        pos = np.flatnonzero(gse_values == gse_id)
        available = source_gene_symbols(Path(str(source_by_gse.loc[gse_id])))
        local = logx[pos, :]
        for class_index, class_name in enumerate(classes):
            requested = MAJOR_MARKERS[class_name]
            present = [gene for gene in requested if gene in available and gene in gene_to_col]
            availability_rows.append({
                "gse_id": gse_id, "major_celltype": class_name,
                "n_requested": len(requested), "n_available": len(present),
                "available_genes": ",".join(present),
                "missing_genes": ",".join(gene for gene in requested if gene not in present),
            })
            if len(present) < 2:
                continue
            cols = [gene_to_col[gene] for gene in present]
            local_score = np.asarray(local[:, cols].mean(axis=1)).ravel().astype(np.float32)
            scores[pos, class_index] = local_score

            excl_requested = MAJOR_EXCLUDE_MARKERS.get(class_name, ())
            excl_present = [gene for gene in excl_requested if gene in available and gene in gene_to_col]
            if excl_present:
                excl_cols = [gene_to_col[gene] for gene in excl_present]
                exclude_scores[pos, class_index] = np.asarray(local[:, excl_cols].mean(axis=1)).ravel().astype(np.float32)
            gse_rows.append({
                "gse_id": gse_id, "major_celltype": class_name, "n_cells": len(pos),
                "n_markers_available": len(present), "mean_marker_score": float(local_score.mean()),
                "median_marker_score": float(np.median(local_score)),
                "fraction_score_ge_0_25": float((local_score >= 0.25).mean()),
            })
        logger.info("marker modules complete gse=%s cells=%s", gse_id, f"{len(pos):,}")

    ranking = np.where(np.isfinite(scores), scores, -np.inf)
    top_idx = np.argmax(ranking, axis=1)
    top_score = ranking[np.arange(len(ranking)), top_idx]
    second_matrix = ranking.copy()
    second_matrix[np.arange(len(ranking)), top_idx] = -np.inf
    second_idx = np.argmax(second_matrix, axis=1)
    second_score = second_matrix[np.arange(len(ranking)), second_idx]
    top_name = np.asarray(classes, dtype=object)[top_idx]
    second_name = np.asarray(classes, dtype=object)[second_idx]
    margin = top_score - second_score
    top_name[~np.isfinite(top_score)] = "Unknown"
    second_name[~np.isfinite(second_score)] = "Unknown"
    top_score[~np.isfinite(top_score)] = np.nan
    second_score[~np.isfinite(second_score)] = np.nan
    margin[~np.isfinite(margin)] = np.nan

    excl_score_top = exclude_scores[np.arange(len(ranking)), top_idx]
    exclusion_ok = ~np.isfinite(excl_score_top) | (excl_score_top < float(config["exclusion_score_max"]))

    marker_strong_no_excl = np.isfinite(top_score) & (top_score >= float(config["marker_score_min"])) & (margin >= float(config["marker_margin_min"]))
    marker_strong = marker_strong_no_excl & exclusion_ok
    excl_flipped = marker_strong_no_excl & ~exclusion_ok
    logger.info("exclusion-marker gate: %s cells that passed v2 positive-marker criteria were disqualified for high exclusion-marker expression", f"{int(excl_flipped.sum()):,}")

    original = evidence["original_major"].astype("string").fillna("").to_numpy(dtype=object)
    ct = evidence["celltypist_major"].astype("string").fillna("").to_numpy(dtype=object)
    conf = evidence["celltypist_confidence"].to_numpy(dtype=np.float32)
    ct_confident = conf >= float(config["celltypist_confidence_min"])
    final = np.full(len(evidence), "Unknown", dtype=object)
    source = np.full(len(evidence), "insufficient_evidence", dtype=object)

    for i in range(len(evidence)):
        orig = str(original[i])
        ctm = str(ct[i])
        marker = str(top_name[i])
        has_orig = bool(orig)
        marker_ok = bool(marker_strong[i])
        ct_ok = bool(ct_confident[i]) and bool(ctm)
        if has_orig:
            if marker_ok and marker != orig:
                final[i], source[i] = "Ambiguous", "original_marker_conflict"
            elif ct_ok and ctm != orig and (not marker_ok or marker != orig):
                final[i], source[i] = "Ambiguous", "original_celltypist_conflict"
            else:
                final[i] = orig
                source[i] = "original_plus_marker" if marker_ok else ("original_plus_celltypist" if ct_ok and ctm == orig else "original_per_cell")
        elif marker_ok:
            if ct_ok and ctm != marker:
                final[i], source[i] = "Ambiguous", "marker_celltypist_conflict"
            else:
                final[i], source[i] = marker, "marker_plus_celltypist" if ct_ok else "marker_panel"
        elif ct_ok and conf[i] >= 0.80:
            final[i], source[i] = ctm, "celltypist_high_confidence_only"

    evidence["marker_major_top1"] = top_name
    evidence["marker_major_top2"] = second_name
    evidence["marker_score_top1"] = top_score.astype(np.float32)
    evidence["marker_score_top2"] = second_score.astype(np.float32)
    evidence["marker_score_margin"] = margin.astype(np.float32)
    evidence["marker_exclusion_score_top1"] = excl_score_top.astype(np.float32)
    evidence["marker_strong"] = marker_strong
    evidence["excluded_by_negative_marker"] = excl_flipped
    evidence["major_celltype_evidence"] = final
    evidence["major_evidence_source"] = source

    l1_rows: list[dict[str, object]] = []
    l1_values = evidence["l1_cluster"].astype(str).to_numpy()
    for l1 in sorted(np.unique(l1_values)):
        inside = l1_values == l1
        outside = ~inside
        for j, class_name in enumerate(classes):
            valid_in = inside & np.isfinite(scores[:, j])
            valid_out = outside & np.isfinite(scores[:, j])
            mean_in = float(np.mean(scores[valid_in, j])) if valid_in.any() else np.nan
            mean_out = float(np.mean(scores[valid_out, j])) if valid_out.any() else np.nan
            l1_rows.append({
                "l1_cluster": l1, "major_celltype": class_name, "n_cells": int(inside.sum()),
                "n_valid_in": int(valid_in.sum()), "mean_score_in": mean_in, "mean_score_out": mean_out,
                "delta_mean_score": mean_in - mean_out if np.isfinite(mean_in) and np.isfinite(mean_out) else np.nan,
                "fraction_marker_top1_in": float((top_name[inside] == class_name).mean()),
                "fraction_final_label_in": float((final[inside] == class_name).mean()),
            })

    evidence.to_pickle(output, compression="gzip")
    logger.info("wrote major evidence checkpoint %s", output)
    (evidence.groupby(["major_celltype_evidence", "major_evidence_source"], observed=True)
     .size().rename("n_cells").reset_index().to_csv(counts_path, sep="\t", index=False))
    pd.DataFrame(availability_rows).to_csv(availability_path, sep="\t", index=False)
    pd.DataFrame(l1_rows).to_csv(l1_path, sep="\t", index=False)
    pd.DataFrame(gse_rows).to_csv(gse_path, sep="\t", index=False)
    flip_summary = (pd.DataFrame({"top_candidate": top_name[excl_flipped], "gse_id": gse_values[excl_flipped]})
                     .value_counts().rename("n_cells").reset_index().sort_values("n_cells", ascending=False))
    flip_summary.to_csv(exclude_path, sep="\t", index=False)

    ended_iso = now_iso()
    duration = time.time() - started
    value_counts = pd.Series(final).value_counts().to_dict()
    payload = {
        "run_id": config["run_id"], "stage": "02_major_evidence_v3", "status": "success",
        "start_time": started_iso, "end_time": ended_iso, "duration_sec": duration,
        "params": {k: config[k] for k in ("celltypist_confidence_min", "marker_score_min", "marker_margin_min", "exclusion_score_max")},
        "inputs": [str(originals_path), str(resolve_project_path(config["legacy_wide_panel_full_retry"])), str(resolve_project_path(config["celltypist_predictions"]))],
        "outputs": [str(output), str(counts_path), str(availability_path), str(l1_path), str(gse_path), str(exclude_path)],
        "warnings": [
            "CellTypist Immune_All_Low is supporting evidence, not ground truth",
            "L1 tables report effect sizes and replication support; pooled cell-level p-values are intentionally not used",
            "Unknown and Ambiguous cells are retained for later review",
            "v3 adds two exclusion-marker fixes vs v2: Neural gets PMEL/MLANA/TYR/DCT "
            "(melanocyte/melanoma guard, confirmed cause of Neural_Melanocytic_PMEL "
            "being 97.9% Uveal melanoma), HSPC gets PECAM1/VWF/COL1A1 "
            "(CD34+ endothelial/fibroblast guard, precautionary). See "
            "MARKER_COVERAGE_AUDIT_20260825.md for the full audit.",
        ],
        "error": None, "source_scripts": [str(Path(__file__).resolve())], "logic_version": config["logic_version"],
        "hostname": socket.gethostname(), "pid": os.getpid(), "n_cells": int(len(evidence)),
        "major_counts": {str(k): int(v) for k, v in value_counts.items()},
        "n_flipped_by_exclusion_marker": int(excl_flipped.sum()),
    }
    write_json_new(summary_path, payload)
    append_status("02v3", "major_evidence_v3", "success", started_iso, ended_iso, round(duration, 3),
                  str(originals_path), str(output), str(os.getpid()), "")
    print(payload)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logging.getLogger("major_evidence_v3").exception("stage failed")
        append_status("02v3", "major_evidence_v3", "failed", error_log=repr(exc), job_id=str(os.getpid()))
        raise
