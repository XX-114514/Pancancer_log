#!/usr/bin/env python3
"""Rescore known-major V5 unresolved clusters with the complete subtype panels.

This is append-only and cluster-level.  The complete panel nominates a subtype,
but a rescue additionally requires a canonical anchor and sample replication.
Contamination panels are reported as conflicts and never become mature labels.
"""

from __future__ import annotations

import json
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from common_v2 import ROOT, load_config, now_iso, project_path, refuse_existing, write_json


LEGACY_SCRIPTS = project_path("runs/20260725_120105_unified_pancancer_covarnet/scripts")
sys.path.insert(0, str(LEGACY_SCRIPTS))
from pancancer_evidence_marker_specs import SPECS, is_blacklisted_marker  # noqa: E402


ALIASES = {"KIM1": "HAVCR1"}


def tier_for(n_lineage: int) -> dict[str, float | int | str]:
    if n_lineage < 5_000:
        return {"tier": "rare_coarse", "enabled": 0, "min_cells": 0,
                "min_samples": 0, "min_full": 0, "min_margin": 0.0}
    if n_lineage < 20_000:
        return {"tier": "low_abundance", "enabled": 1, "min_cells": 80,
                "min_samples": 2, "min_full": 2, "min_margin": 0.010}
    if n_lineage < 100_000:
        return {"tier": "medium_abundance", "enabled": 1, "min_cells": 60,
                "min_samples": 2, "min_full": 2, "min_margin": 0.010}
    return {"tier": "high_abundance", "enabled": 1, "min_cells": 100,
            "min_samples": 3, "min_full": 2, "min_margin": 0.015}


def main() -> None:
    started = time.time()
    cfg = load_config()
    cluster_out = ROOT / "tables" / "expanded_panel_cluster_rescue_v1.tsv"
    gene_out = ROOT / "tables" / "expanded_panel_gene_evidence_v1.tsv.gz"
    counts_out = ROOT / "tables" / "expanded_panel_rescue_counts_v1.tsv"
    summary_out = ROOT / "reports" / "02_expanded_panel_rescore_summary_v1.json"
    refuse_existing(cluster_out, gene_out, counts_out, summary_out)

    final = pd.read_pickle(project_path(cfg["final_v5"]))
    wide = ad.read_h5ad(project_path(cfg["wide_marker_panel"]))
    if not final.index.equals(wide.obs_names):
        raise ValueError("Final V5 and wide marker object differ in cell order")
    required = {"cluster_uid", "major_celltype_final", "annotation_confidence", "gse_id", "sampleID"}
    missing = required.difference(final.columns)
    if missing:
        raise KeyError(f"Final V5 lacks columns: {sorted(missing)}")

    X = wide.X.tocsr() if sp.issparse(wide.X) else sp.csr_matrix(wide.X)
    genes = pd.Index(wide.var_names.astype(str))
    gene_to_col = {gene: i for i, gene in enumerate(genes)}
    availability_df = wide.varm["availability_by_gse"]
    gse_order = list(availability_df.columns.astype(str))
    availability_gse = availability_df.to_numpy(dtype=bool).T

    cluster_cat = pd.Categorical(final["cluster_uid"].astype(str))
    cluster_names = np.asarray(cluster_cat.categories.astype(str))
    cluster_codes = cluster_cat.codes.astype(np.int32)
    n_clusters = len(cluster_names)
    cluster_counts = np.bincount(cluster_codes, minlength=n_clusters).astype(np.int64)

    # sampleID is not assumed globally unique; GSE is part of the replication unit.
    sample_uid = final["gse_id"].astype(str) + "::" + final["sampleID"].astype(str)
    sample_cat = pd.Categorical(sample_uid)
    sample_names = np.asarray(sample_cat.categories.astype(str))
    sample_codes = sample_cat.codes.astype(np.int32)
    n_samples = len(sample_names)
    sample_gse = (pd.DataFrame({"sample_uid": sample_uid, "gse": final["gse_id"].astype(str)})
                  .drop_duplicates("sample_uid").set_index("sample_uid")["gse"].reindex(sample_names))
    if sample_gse.isna().any():
        raise ValueError("Some sample units have no GSE")
    sample_gse_codes = pd.Index(gse_order).get_indexer(sample_gse.to_numpy())
    if (sample_gse_codes < 0).any():
        raise ValueError("Some samples refer to an unavailable GSE")
    availability_sample = availability_gse[sample_gse_codes]

    rows = np.arange(len(final), dtype=np.int64)
    selector = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cluster_codes)),
                             shape=(len(rows), n_clusters))
    sums = np.asarray((selector.T @ X).toarray(), dtype=np.float64)
    binary = X.copy()
    binary.data.fill(1.0)
    nnz = np.asarray((selector.T @ binary).toarray(), dtype=np.float64)

    cs_codes = cluster_codes.astype(np.int64) * n_samples + sample_codes
    cs_selector = sp.csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cs_codes)),
                                shape=(len(rows), n_clusters * n_samples))
    cs_sums = np.asarray((cs_selector.T @ X).toarray(), dtype=np.float32).reshape(
        n_clusters, n_samples, X.shape[1]
    )
    cs_counts = np.bincount(cs_codes, minlength=n_clusters * n_samples).reshape(
        n_clusters, n_samples
    ).astype(np.int64)
    eligible = np.asarray(cs_counts @ availability_sample, dtype=np.float64)
    means = np.divide(sums, eligible, out=np.zeros_like(sums), where=eligible > 0)
    fractions = np.divide(nnz, eligible, out=np.zeros_like(nnz), where=eligible > 0)

    cluster_frame = (final.groupby("cluster_uid", observed=True)
                     .agg(major=("major_celltype_final", "first"),
                          confidence=("annotation_confidence", "first"),
                          old_label=("final_annotation", "first"))
                     .reindex(cluster_names))
    if cluster_frame.isna().any().any():
        raise ValueError("Could not map all clusters to V5 labels")
    cluster_major = cluster_frame["major"].astype(str).to_numpy()
    unresolved_cluster = cluster_frame["confidence"].astype(str).eq("unresolved").to_numpy()

    cluster_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    for major in sorted(set(cluster_major).intersection(SPECS)):
        lineage_idx = np.flatnonzero(cluster_major == major)
        target_idx = lineage_idx[unresolved_cluster[lineage_idx]]
        if len(target_idx) == 0:
            continue
        n_lineage = int(cluster_counts[lineage_idx].sum())
        tier = tier_for(n_lineage)
        lineage_sums = sums[lineage_idx].sum(axis=0)
        lineage_nnz = nnz[lineage_idx].sum(axis=0)
        lineage_eligible = eligible[lineage_idx].sum(axis=0)
        lineage_cs_sums = cs_sums[lineage_idx].sum(axis=0)
        lineage_cs_counts = cs_counts[lineage_idx].sum(axis=0)
        spec = SPECS[major]

        for ci in target_idx:
            out_eligible = lineage_eligible - eligible[ci]
            mean_out = np.divide(lineage_sums - sums[ci], out_eligible,
                                 out=np.zeros_like(lineage_sums), where=out_eligible > 0)
            pct_out = np.divide(lineage_nnz - nnz[ci], out_eligible,
                                out=np.zeros_like(lineage_nnz), where=out_eligible > 0)
            log2fc = (means[ci] - mean_out) / np.log(2)
            delta = fractions[ci] - pct_out
            results: list[dict[str, object]] = []

            for subtype, requested_panel in spec.panels.items():
                full_requested = [ALIASES.get(g, g) for g in requested_panel]
                canonical_requested = [ALIASES.get(g, g) for g in spec.canonical.get(subtype, [])]
                full_cols = [gene_to_col[g] for g in full_requested
                             if g in gene_to_col and not is_blacklisted_marker(g)]
                canonical_set = {gene_to_col[g] for g in canonical_requested
                                 if g in gene_to_col and not is_blacklisted_marker(g)}
                passed: list[int] = []
                canonical_passed: list[int] = []
                values: list[float] = []
                for col in full_cols:
                    is_pass = bool(fractions[ci, col] >= spec.min_pct_in
                                   and log2fc[col] >= spec.min_log2fc
                                   and delta[col] >= spec.min_delta_pct
                                   and eligible[ci, col] > 0)
                    score = float(fractions[ci, col] * max(log2fc[col], 0.0) * max(delta[col], 0.0))
                    if is_pass:
                        passed.append(col)
                        values.append(score)
                        if col in canonical_set:
                            canonical_passed.append(col)
                    gene_rows.append({
                        "cluster_uid": cluster_names[ci], "major_celltype": major,
                        "subtype": subtype, "gene": genes[col],
                        "is_canonical": col in canonical_set,
                        "n_eligible_cells": int(eligible[ci, col]),
                        "pct_in": float(fractions[ci, col]), "pct_out": float(pct_out[col]),
                        "mean_in": float(means[ci, col]), "mean_out": float(mean_out[col]),
                        "log2fc": float(log2fc[col]), "delta_pct": float(delta[col]),
                        "passes_threshold": is_pass, "gene_evidence_score": score,
                    })
                results.append({
                    "subtype": subtype,
                    "score": float(np.mean(values)) if values else -np.inf,
                    "passed": passed,
                    "canonical_passed": canonical_passed,
                    "full_available": full_cols,
                    "canonical_available": list(canonical_set),
                })

            results.sort(key=lambda item: float(item["score"]), reverse=True)
            top = results[0]
            second = results[1] if len(results) > 1 else {"subtype": "", "score": -np.inf}
            top_score = float(top["score"])
            second_score = float(second["score"])
            margin = (top_score - (second_score if np.isfinite(second_score) else 0.0)
                      if np.isfinite(top_score) else np.nan)
            top_passed = list(top["passed"])
            top_canonical = list(top["canonical_passed"])

            support = 0
            evaluable_sample = 0
            for sj in range(n_samples):
                in_n = int(cs_counts[ci, sj])
                out_n = int(lineage_cs_counts[sj] - in_n)
                available_passed = [col for col in top_passed if availability_sample[sj, col]]
                available_canonical = [col for col in top_canonical if availability_sample[sj, col]]
                if in_n < 10 or out_n < 10 or len(available_passed) < int(tier["min_full"]) \
                        or len(available_canonical) < 1:
                    continue
                in_mean = cs_sums[ci, sj, available_passed] / in_n
                out_mean = (lineage_cs_sums[sj, available_passed]
                            - cs_sums[ci, sj, available_passed]) / out_n
                evaluable_sample += 1
                support += int(float(np.mean(in_mean - out_mean)) > 0.10)

            contamination = "contamination" in str(top["subtype"]).lower()
            passes = bool(int(tier["enabled"])
                          and cluster_counts[ci] >= int(tier["min_cells"])
                          and len(top_passed) >= int(tier["min_full"])
                          and len(top_canonical) >= 1
                          and np.isfinite(margin) and margin >= float(tier["min_margin"])
                          and support >= int(tier["min_samples"])
                          and not contamination)
            if passes:
                best_canonical = max(top_canonical, key=lambda c: float(
                    fractions[ci, c] * max(log2fc[c], 0.0) * max(delta[c], 0.0)
                ))
                proposed = f"{major}_{top['subtype']}_{genes[best_canonical]}"
                reason = "expanded_panel_plus_canonical_anchor_plus_sample_replication"
            elif contamination:
                proposed = f"{major}_unresolved"
                reason = "top_expanded_panel_is_contamination_conflict"
            elif len(top_canonical) == 0:
                proposed = f"{major}_unresolved"
                reason = "expanded_signal_without_canonical_anchor"
            else:
                proposed = f"{major}_unresolved"
                reason = "expanded_panel_threshold_or_replication_not_met"

            cluster_rows.append({
                "cluster_uid": cluster_names[ci], "major_celltype": major,
                "n_cells": int(cluster_counts[ci]), "n_lineage_cells": n_lineage,
                "tier": tier["tier"], "old_label": cluster_frame.iloc[ci]["old_label"],
                "top_subtype": top["subtype"],
                "top_score": top_score if np.isfinite(top_score) else np.nan,
                "second_subtype": second["subtype"],
                "second_score": second_score if np.isfinite(second_score) else np.nan,
                "margin": margin, "n_full_panel_available": len(top["full_available"]),
                "n_full_panel_pass": len(top_passed),
                "n_canonical_available": len(top["canonical_available"]),
                "n_canonical_pass": len(top_canonical),
                "full_panel_pass_genes": ",".join(genes[top_passed]),
                "canonical_pass_genes": ",".join(genes[top_canonical]),
                "n_evaluable_sample": evaluable_sample, "n_support_sample": support,
                "rescued": passes, "proposed_label": proposed, "reason": reason,
            })

    clusters = pd.DataFrame(cluster_rows)
    gene_evidence = pd.DataFrame(gene_rows)
    clusters.to_csv(cluster_out, sep="\t", index=False)
    gene_evidence.to_csv(gene_out, sep="\t", index=False, compression="gzip")
    count_table = (clusters.groupby(["major_celltype", "rescued", "reason"], observed=True)
                   .agg(n_clusters=("cluster_uid", "size"), n_cells=("n_cells", "sum"))
                   .reset_index().sort_values("n_cells", ascending=False))
    count_table.to_csv(counts_out, sep="\t", index=False)
    payload = {
        "run_id": cfg["run_id"], "stage": "02_rescore_unresolved_with_expanded_panels",
        "status": "success", "start_time": now_iso(),
        "duration_sec": round(time.time() - started, 3),
        "n_unresolved_clusters_evaluated": int(len(clusters)),
        "n_unresolved_cells_evaluated": int(clusters["n_cells"].sum()),
        "n_clusters_rescued": int(clusters["rescued"].sum()),
        "n_cells_rescued": int(clusters.loc[clusters["rescued"], "n_cells"].sum()),
        "scoring_change": "full subtype panel nominates; >=1 canonical anchor and sample replication required",
        "warnings": [
            "Candidate rescue table only; final V5 is unchanged.",
            "Contamination panels are conflict evidence and cannot name a mature subtype.",
            "The score uses the existing 359-gene coverage-masked object; genome-wide de novo markers are a separate stage.",
        ],
        "outputs": [str(cluster_out), str(gene_out), str(counts_out)],
    }
    write_json(summary_out, payload)
    print(json.dumps(payload, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
