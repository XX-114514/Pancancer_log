#!/usr/bin/env python3
"""v5: same inputs as v4 (major_evidence_assignments_v3.pkl.gz + the BBKNN
lineage clusters selected_lineage_mapping_v3.tsv points to -- 05c is NOT
rerun, its clusters are reused as-is), fixing one bug found by manually
auditing v4's output (MARKER_COVERAGE_AUDIT_20260825.md section 4b):

    margin = top[1] - second[1] if np.isfinite(top[1]) and np.isfinite(second[1]) else np.nan

When a cluster's runner-up subtype candidate has ZERO marker genes clearing
the pct_in/log2FC/delta_pct gates, its score is -inf (by construction, a few
lines above: `float(np.mean(values)) if values else -np.inf`). The old code
required BOTH top and second scores to be finite before computing a margin,
so a cluster where exactly one subtype has any passing evidence and every
other candidate has none -- arguably the single LEAST ambiguous case
possible -- got margin=NaN, which `passes = ... and np.isfinite(margin) and
margin >= tier["min_margin"] ...` then treats as failure alongside "no
information available".

Concretely this hid a real signal: v4's Neural_1 cluster (62,127 cells, the
real Glioblastoma NPC-like population, Neftel et al. 2019) scored top_score=
0.31 with n_pass_genes=3 (tier requires >=2) and 42/49 samples supporting
(tier requires >=2) -- comfortably clearing every other gate -- yet was
marked unresolved purely because no other Neural subtype had any passing
gene in that cluster, so second was -inf and margin was NaN. (In v3 this was
masked: leftover melanocytic cells sharing that cluster gave the Melanocytic
candidate a nonzero score, which is what let margin compute at all -- fixing
the melanocyte contamination is what exposed this latent bug, not what
caused it.)

Fix: treat "no gene passed for the runner-up" as a runner-up score of 0.0
(no signal), not -inf, for the margin computation only (ranking/sorting
logic is untouched -- a candidate with zero passing genes still ranks below
any candidate with a positive score). This makes margin = top_score itself
when there is no viable competitor, which is the correct behavior: an
unopposed candidate is a strong signal, not an unscoreable one."""

from __future__ import annotations

import logging
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import zarr

from common import PROJECT, ROOT, append_status, ensure_layout, load_config, now_iso, refuse_existing, resolve_project_path, write_json_new


LEGACY_SCRIPTS = PROJECT / "runs" / "20260725_120105_unified_pancancer_covarnet" / "scripts"
sys.path.insert(0, str(LEGACY_SCRIPTS))
from pancancer_evidence_marker_specs import SPECS, is_blacklisted_marker  # noqa: E402

ALIASES = {"KIM1": "HAVCR1"}


def tier_for(n_lineage: int) -> dict[str, object]:
    if n_lineage < 5_000:
        return {"tier": "rare_coarse", "subtype_enabled": False, "min_cluster_cells": 0,
                "min_support_sample": 0, "min_pass_genes": 0, "min_margin": 0.0}
    if n_lineage < 20_000:
        return {"tier": "low_abundance", "subtype_enabled": True, "min_cluster_cells": 80,
                "min_support_sample": 2, "min_pass_genes": 2, "min_margin": 0.010}
    if n_lineage < 100_000:
        return {"tier": "medium_abundance", "subtype_enabled": True, "min_cluster_cells": 60,
                "min_support_sample": 2, "min_pass_genes": 2, "min_margin": 0.010}
    return {"tier": "high_abundance", "subtype_enabled": True, "min_cluster_cells": 100,
            "min_support_sample": 3, "min_pass_genes": 2, "min_margin": 0.015}


def boolish(series: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    text = series.astype("string").str.lower()
    valid = text.isin(["true", "false"]).to_numpy()
    value = text.eq("true").fillna(False).to_numpy()
    return value, valid


def main() -> None:
    started = time.time()
    started_iso = now_iso()
    ensure_layout()
    config = load_config()
    selection = pd.read_csv(ROOT / "tables" / "selected_lineage_mapping_v3.tsv", sep="\t").iloc[0]
    mapping_path = Path(str(selection["mapping_path"]))
    if not mapping_path.is_absolute():
        mapping_path = resolve_project_path(mapping_path)
    wide_path = ROOT / "objects" / "wide_marker_panel_full_library_masked_v2.h5ad"
    major_path = ROOT / "checkpoints" / "major_evidence_assignments_v3.pkl.gz"
    final_path = ROOT / "checkpoints" / "final_evidence_annotations_v5.pkl.gz"
    obs_h5ad = ROOT / "objects" / "final_annotation_obs_only_v5.h5ad"
    cluster_path = ROOT / "tables" / "cluster_annotation_evidence_v5.tsv"
    gene_path = ROOT / "tables" / "cluster_marker_effect_sizes_v5.tsv.gz"
    threshold_path = ROOT / "tables" / "lineage_specific_thresholds_v5.tsv"
    counts_path = ROOT / "tables" / "final_annotation_counts_v5.tsv"
    summary_path = ROOT / "reports" / "06_final_evidence_annotation_v5_summary.json"
    refuse_existing(final_path, obs_h5ad, cluster_path, gene_path, threshold_path, counts_path, summary_path)
    append_status("06v5", "final_evidence_annotation_v5", "running", started_iso,
                  input_ref=str(mapping_path), job_id=str(os.getpid()))

    major = pd.read_pickle(major_path)
    lineage = pd.read_pickle(mapping_path)
    if not major.index.equals(lineage.index):
        raise ValueError("Major and selected lineage mappings differ")
    wide = ad.read_h5ad(wide_path)
    if not major.index.equals(wide.obs_names):
        raise ValueError("Corrected wide marker object is not in merge order")
    X = wide.X.tocsr() if sp.issparse(wide.X) else sp.csr_matrix(wide.X)
    gene_names = pd.Index(wide.var_names.astype(str))
    gene_to_col = {gene: i for i, gene in enumerate(gene_names)}
    availability_df = wide.varm["availability_by_gse"]
    gse_order = list(availability_df.columns.astype(str))
    availability = availability_df.to_numpy(dtype=bool).T

    cluster_cat = pd.Categorical(lineage["cluster_uid"].astype(str))
    cluster_names = np.asarray(cluster_cat.categories.astype(str))
    cluster_codes = cluster_cat.codes.astype(np.int32)
    n_clusters = len(cluster_names)
    gse_cat = pd.Categorical(major["gse_id"].astype(str), categories=gse_order)
    if (gse_cat.codes < 0).any():
        raise ValueError("Some cells have a GSE absent from availability mask")

    sample_cat = pd.Categorical(major["sampleID"].astype(str))
    sample_names = np.asarray(sample_cat.categories.astype(str))
    sample_codes = sample_cat.codes.astype(np.int32)
    n_sample = len(sample_names)
    sample_to_gse = (major.assign(_s=major["sampleID"].astype(str), _g=major["gse_id"].astype(str))
                     .drop_duplicates("_s").set_index("_s")["_g"].reindex(sample_names))
    if sample_to_gse.isna().any():
        raise ValueError("Some samples have no resolvable GSE for availability lookup")
    gse_index_of_sample = pd.Index(gse_order).get_indexer(sample_to_gse.to_numpy())
    if (gse_index_of_sample < 0).any():
        raise ValueError("Some samples map to a GSE absent from the availability mask")
    availability_sample = availability[gse_index_of_sample]  # (n_sample, n_genes)

    rows = np.arange(len(major), dtype=np.int64)
    selector = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cluster_codes)),
                             shape=(len(rows), n_clusters))
    sums = np.asarray((selector.T @ X).toarray(), dtype=np.float64)
    binary = X.copy()
    binary.data.fill(1.0)
    nnz = np.asarray((selector.T @ binary).toarray(), dtype=np.float64)
    cluster_counts = np.bincount(cluster_codes, minlength=n_clusters).astype(np.int64)
    cg_codes = cluster_codes.astype(np.int64) * n_sample + sample_codes
    cg_selector = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cg_codes)),
                                shape=(len(rows), n_clusters * n_sample))
    cg_sums = np.asarray((cg_selector.T @ X).toarray(), dtype=np.float32).reshape(n_clusters, n_sample, X.shape[1])
    cg_counts = np.bincount(cg_codes, minlength=n_clusters * n_sample).reshape(n_clusters, n_sample).astype(np.int64)
    eligible = np.asarray(cg_counts @ availability_sample, dtype=np.float64)
    means = np.divide(sums, eligible, out=np.zeros_like(sums), where=eligible > 0)
    fractions = np.divide(nnz, eligible, out=np.zeros_like(nnz), where=eligible > 0)

    cluster_major = (lineage.groupby("cluster_uid", observed=True)["major_celltype_evidence"].first()
                     .reindex(cluster_names).astype(str).to_numpy())
    threshold_rows: list[dict[str, object]] = []
    cluster_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    cluster_final: dict[str, str] = {}
    cluster_conf: dict[str, str] = {}

    for major_name in sorted(np.unique(cluster_major)):
        cluster_idx = np.flatnonzero(cluster_major == major_name)
        n_lineage = int(cluster_counts[cluster_idx].sum())
        tier = tier_for(n_lineage)
        threshold_rows.append({"major_celltype": major_name, "n_cells": n_lineage, **tier})
        if major_name in {"Unknown", "Ambiguous", "Epithelial"} or major_name not in SPECS or not tier["subtype_enabled"]:
            for ci in cluster_idx:
                if major_name in {"Unknown", "Ambiguous"}:
                    label, confidence, reason = major_name, "unresolved", "broad_evidence_unresolved"
                elif major_name == "Epithelial":
                    label, confidence, reason = "Epithelial_CNV_pending", "provisional", "resolved_per_cell_below"
                else:
                    label, confidence, reason = f"{major_name}_coarse", "coarse", "rare_lineage_no_forced_subtype"
                cluster_final[cluster_names[ci]] = label
                cluster_conf[cluster_names[ci]] = confidence
                cluster_rows.append({"cluster_uid": cluster_names[ci], "major_celltype": major_name,
                                     "n_cells": int(cluster_counts[ci]), "tier": tier["tier"],
                                     "top_subtype": "", "top_score": np.nan, "second_score": np.nan,
                                     "margin": np.nan, "n_pass_genes": 0, "n_support_sample": 0,
                                     "final_cluster_label": label, "confidence": confidence, "reason": reason})
            continue

        lineage_sums = sums[cluster_idx].sum(axis=0)
        lineage_nnz = nnz[cluster_idx].sum(axis=0)
        lineage_eligible = eligible[cluster_idx].sum(axis=0)
        lineage_cg_sums = cg_sums[cluster_idx].sum(axis=0)
        lineage_cg_counts = cg_counts[cluster_idx].sum(axis=0)
        spec = SPECS[major_name]
        for ci in cluster_idx:
            out_eligible = lineage_eligible - eligible[ci]
            mean_out = np.divide(lineage_sums - sums[ci], out_eligible,
                                 out=np.zeros_like(lineage_sums), where=out_eligible > 0)
            pct_out = np.divide(lineage_nnz - nnz[ci], out_eligible,
                                out=np.zeros_like(lineage_nnz), where=out_eligible > 0)
            log2fc = (means[ci] - mean_out) / np.log(2)
            delta = fractions[ci] - pct_out
            subtype_results: list[tuple[str, float, int, list[int], list[int]]] = []
            for subtype, panel in spec.panels.items():
                requested = spec.canonical.get(subtype, panel)
                resolved = [ALIASES.get(gene, gene) for gene in requested]
                cols = [gene_to_col[gene] for gene in resolved if gene in gene_to_col and not is_blacklisted_marker(gene)]
                passed: list[int] = []
                values: list[float] = []
                for col in cols:
                    is_pass = (fractions[ci, col] >= spec.min_pct_in and log2fc[col] >= spec.min_log2fc
                               and delta[col] >= spec.min_delta_pct and eligible[ci, col] > 0)
                    score = fractions[ci, col] * max(log2fc[col], 0) * max(delta[col], 0)
                    if is_pass:
                        passed.append(col)
                        values.append(float(score))
                    gene_rows.append({
                        "cluster_uid": cluster_names[ci], "major_celltype": major_name, "subtype": subtype,
                        "gene": gene_names[col], "n_eligible_cells": int(eligible[ci, col]),
                        "pct_in": float(fractions[ci, col]), "pct_out": float(pct_out[col]),
                        "mean_in": float(means[ci, col]), "mean_out": float(mean_out[col]),
                        "log2fc": float(log2fc[col]), "delta_pct": float(delta[col]),
                        "passes_threshold": bool(is_pass), "gene_evidence_score": float(score),
                    })
                subtype_results.append((subtype, float(np.mean(values)) if values else -np.inf,
                                        len(passed), cols, passed))
            subtype_results.sort(key=lambda item: item[1], reverse=True)
            top = subtype_results[0]
            second = subtype_results[1] if len(subtype_results) > 1 else ("", -np.inf, 0, [], [])
            # v5 fix: a runner-up with zero passing genes (-inf) means "no
            # competing signal detected", not "margin unscoreable" -- treat
            # it as a baseline score of 0.0 for the margin gap only. Ranking
            # above is untouched (a -inf-scored candidate still never wins
            # `top`); this only changes what "distance from the runner-up"
            # means when there IS no runner-up.
            if np.isfinite(top[1]):
                second_for_margin = second[1] if np.isfinite(second[1]) else 0.0
                margin = top[1] - second_for_margin
            else:
                margin = np.nan

            support = 0
            evaluable_sample = 0
            for sj in range(n_sample):
                in_n = int(cg_counts[ci, sj])
                out_n = int(lineage_cg_counts[sj] - in_n)
                cols = [col for col in top[4] if availability_sample[sj, col]]
                if in_n < 10 or out_n < 10 or len(cols) < int(tier["min_pass_genes"]):
                    continue
                in_mean = cg_sums[ci, sj, cols] / in_n
                out_mean = (lineage_cg_sums[sj, cols] - cg_sums[ci, sj, cols]) / out_n
                evaluable_sample += 1
                support += int(float(np.mean(in_mean - out_mean)) > 0.10)

            passes = (cluster_counts[ci] >= int(tier["min_cluster_cells"])
                      and top[2] >= int(tier["min_pass_genes"])
                      and np.isfinite(top[1]) and np.isfinite(margin)
                      and margin >= float(tier["min_margin"])
                      and support >= int(tier["min_support_sample"]))
            if passes:
                marker_gene = gene_names[top[4][0]] if len(top[4]) > 0 else "NA"
                label = f"{major_name}_{top[0]}_{marker_gene}"
                confidence = "high" if support >= max(3, int(tier["min_support_sample"])) and margin >= 0.03 else "medium"
                reason = "coverage_aware_markers_plus_sample_replication"
            else:
                label, confidence, reason = f"{major_name}_unresolved", "unresolved", "threshold_or_replication_not_met"
            cluster_final[cluster_names[ci]] = label
            cluster_conf[cluster_names[ci]] = confidence
            cluster_rows.append({
                "cluster_uid": cluster_names[ci], "major_celltype": major_name,
                "n_cells": int(cluster_counts[ci]), "tier": tier["tier"],
                "top_subtype": top[0], "top_score": top[1] if np.isfinite(top[1]) else np.nan,
                "second_score": second[1] if np.isfinite(second[1]) else np.nan, "margin": margin,
                "n_pass_genes": top[2], "n_evaluable_sample": evaluable_sample, "n_support_sample": support,
                "final_cluster_label": label, "confidence": confidence, "reason": reason,
            })

    final = major.copy()
    final["integration_basis_selected"] = str(selection["selected_method"])
    final["cluster_uid"] = lineage["cluster_uid"].astype(str).to_numpy()
    final["final_annotation"] = final["cluster_uid"].map(cluster_final).astype(str)
    final["annotation_confidence"] = final["cluster_uid"].map(cluster_conf).astype(str)

    malignancy_path = resolve_project_path(
        "runs/20260818_pancancer_5M_merge_prep/06_malignancy_new_cohort/manifests/combined_malignancy_evidence.tsv")
    malignancy = pd.read_csv(malignancy_path, sep="\t", index_col=0, low_memory=False).reindex(final.index)
    intersection, intersection_valid = boolish(malignancy["malignant_intersection"])
    union, union_valid = boolish(malignancy["malignant_union"])
    evaluable, evaluable_valid = boolish(malignancy["copykat_evaluable"])
    final["malignant_intersection"] = pd.array(np.where(intersection_valid, intersection, pd.NA), dtype="boolean")
    final["malignant_union"] = pd.array(np.where(union_valid, union, pd.NA), dtype="boolean")
    final["copykat_evaluable"] = pd.array(np.where(evaluable_valid, evaluable, pd.NA), dtype="boolean")
    epi = final["major_celltype_evidence"].astype(str).eq("Epithelial").to_numpy()
    labels = final["final_annotation"].to_numpy(dtype=object)
    conf = final["annotation_confidence"].to_numpy(dtype=object)
    high_specificity = epi & intersection_valid & intersection
    discordant = epi & union_valid & union & ~high_specificity
    labels[high_specificity] = "Malignant_Epithelial_high_specificity"
    conf[high_specificity] = "high"
    labels[discordant] = "Epithelial_malignancy_discordant"
    conf[discordant] = "provisional"
    negative = epi & evaluable_valid & evaluable & union_valid & ~union
    labels[negative] = "Nonmalignant_Epithelial_CNV_negative"
    conf[negative] = "medium"
    pending = epi & ~(high_specificity | discordant | negative)
    labels[pending] = "Epithelial_malignancy_not_evaluable"
    conf[pending] = "unresolved"
    final["final_annotation"] = labels
    final["annotation_confidence"] = conf

    merge_root = zarr.open_group(str(resolve_project_path(config["merge_zarr"])), mode="r")
    merge_obs = ad.io.read_elem(merge_root["obs"])
    final["cancer_type"] = merge_obs["cancer_type"].astype(str).to_numpy()
    final.to_pickle(final_path, compression="gzip")
    pd.DataFrame(cluster_rows).to_csv(cluster_path, sep="\t", index=False)
    pd.DataFrame(gene_rows).to_csv(gene_path, sep="\t", index=False, compression="gzip")
    pd.DataFrame(threshold_rows).to_csv(threshold_path, sep="\t", index=False)
    final["final_annotation"].value_counts().rename_axis("final_annotation").rename("n_cells").to_csv(counts_path, sep="\t")

    obs_cols = ["gse_id", "sampleID", "cohort", "cancer_type", "major_celltype_evidence",
                "major_evidence_source", "cluster_uid", "integration_basis_selected",
                "final_annotation", "annotation_confidence", "malignant_intersection",
                "malignant_union", "copykat_evaluable"]
    obs_only = ad.AnnData(X=sp.csr_matrix((len(final), 0), dtype=np.float32), obs=final[obs_cols].copy())
    ad.settings.allow_write_nullable_strings = True
    obs_only.write_h5ad(obs_h5ad, compression="gzip", compression_opts=4)
    ended_iso = now_iso()
    duration = time.time() - started
    write_json_new(summary_path, {
        "run_id": config["run_id"], "stage": "06_final_evidence_annotation_v5", "status": "success",
        "start_time": started_iso, "end_time": ended_iso, "duration_sec": duration,
        "params": {"coverage_masked": True, "replication_unit": "sample_patient_biopsy", "malignancy_rule": "inferCNV_intersection_CopyKAT_high_specificity"},
        "inputs": [str(mapping_path), str(wide_path), str(malignancy_path)],
        "outputs": [str(final_path), str(obs_h5ad), str(cluster_path), str(gene_path), str(threshold_path), str(counts_path)],
        "warnings": ["malignancy union is discordant/provisional", "rare lineages remain coarse", "no pooled million-cell p-values used",
                     "v5 = v4 inputs (major_evidence_assignments_v3.pkl.gz + the v3 BBKNN lineage clusters, NOT rerun) "
                     "+ margin-computation fix: a runner-up subtype with zero passing genes now contributes a "
                     "baseline score of 0.0 (no signal) instead of -inf to the margin gap, instead of forcing "
                     "margin=NaN/failure. See MARKER_COVERAGE_AUDIT_20260825.md section 4b."],
        "error": None, "source_scripts": [str(Path(__file__).resolve())], "logic_version": config["logic_version"],
        "n_cells": int(len(final)), "n_final_labels": int(final["final_annotation"].nunique()),
        "malignant_epithelial_high_specificity": int((labels == "Malignant_Epithelial_high_specificity").sum()),
    })
    append_status("06v5", "final_evidence_annotation_v5", "success", started_iso, ended_iso,
                  round(duration, 3), str(mapping_path), str(final_path), str(os.getpid()), "")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logging.exception("final evidence annotation failed")
        append_status("06v5", "final_evidence_annotation_v5", "failed", error_log=repr(exc), job_id=str(os.getpid()))
        raise
