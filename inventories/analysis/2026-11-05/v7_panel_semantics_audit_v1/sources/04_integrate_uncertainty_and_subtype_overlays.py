#!/usr/bin/env python3
"""Integrate major-refinement and subtype-rescue overlays into a compact 5M table."""

from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

from common_v2 import ROOT, load_config, now_iso, project_path, refuse_existing, write_json


def main() -> None:
    started = time.time()
    cfg = load_config()
    output = ROOT / "checkpoints" / "annotation_uncertainty_resolution_v2.pkl.gz"
    counts = ROOT / "tables" / "integrated_annotation_counts_v2.tsv"
    transitions = ROOT / "tables" / "final_v5_to_uncertainty_v2_transitions.tsv"
    summary = ROOT / "reports" / "04_integrated_uncertainty_annotation_summary_v2.json"
    refuse_existing(output, counts, transitions, summary)

    final = pd.read_pickle(project_path(cfg["final_v5"]))
    major_overlay = pd.read_pickle(ROOT / "checkpoints" / "major_uncertainty_overlay_v2.pkl.gz")
    expanded = pd.read_csv(ROOT / "tables" / "expanded_panel_cluster_rescue_v1.tsv", sep="\t")
    denovo = pd.read_csv(ROOT / "tables" / "de_novo_unresolved_labels_v1.tsv", sep="\t")
    if not major_overlay.index.is_unique or not major_overlay.index.isin(final.index).all():
        raise ValueError("Major overlay index is not a unique V5 subset")

    keep = ["gse_id", "sampleID", "cancer_type", "cohort", "cluster_uid",
            "major_celltype_final", "final_annotation", "annotation_confidence"]
    out = final[keep].copy()
    out = out.rename(columns={
        "major_celltype_final": "major_celltype_v5",
        "final_annotation": "final_annotation_v5",
        "annotation_confidence": "annotation_confidence_v5",
    })
    out["major_celltype_uncertainty_v2"] = out["major_celltype_v5"].astype(str)
    out["major_confidence_uncertainty_v2"] = "unchanged_v5"
    out["major_resolution_reason_v2"] = "unchanged_v5"
    out["major_boundary_parent_v2"] = ""
    out["operational_major_candidate_v2"] = out["major_celltype_v5"].astype(str)
    out["strict_major_eligible_v2"] = ~out["major_celltype_v5"].isin(["Unknown", "Ambiguous"])
    out["artifact_or_state_review_v2"] = "not_reviewed"
    out["exclude_from_primary_major_v2"] = False
    out["original_provenance_class_v2"] = "not_applicable"
    out["conflict_left_v2"] = ""
    out["conflict_right_v2"] = ""
    out["conflict_relation_v2"] = ""
    out["negative_genes_ge_threshold_v2"] = ""

    idx = major_overlay.index
    mapping = {
        "major_celltype_refined_v2": "major_celltype_uncertainty_v2",
        "major_confidence_refined_v2": "major_confidence_uncertainty_v2",
        "major_resolution_reason_v2": "major_resolution_reason_v2",
        "major_boundary_parent_v2": "major_boundary_parent_v2",
        "operational_major_candidate_v2": "operational_major_candidate_v2",
        "strict_major_eligible_v2": "strict_major_eligible_v2",
        "artifact_or_state_review_v2": "artifact_or_state_review_v2",
        "exclude_from_primary_major_v2": "exclude_from_primary_major_v2",
        "original_provenance_class_v2": "original_provenance_class_v2",
        "conflict_left_v2": "conflict_left_v2",
        "conflict_right_v2": "conflict_right_v2",
        "conflict_relation_v2": "conflict_relation_v2",
        "negative_genes_ge_threshold_v2": "negative_genes_ge_threshold_v2",
    }
    for source, target in mapping.items():
        out.loc[idx, target] = major_overlay[source].to_numpy()

    out["final_annotation_uncertainty_v2"] = out["final_annotation_v5"].astype(str)
    out["annotation_confidence_uncertainty_v2"] = out["annotation_confidence_v5"].astype(str)
    out["subtype_resolution_method_v2"] = "unchanged_v5"
    out["subtype_descriptive_marker_v2"] = ""
    out["subtype_top_robust_markers_v2"] = ""
    out["subtype_known_panel_hint_v2"] = ""
    out["strict_subtype_eligible_v2"] = ~out["annotation_confidence_v5"].astype(str).isin(
        ["unresolved", "provisional"]
    )
    out["requires_recluster_for_subtype_v2"] = False

    # Cells whose old broad class was uncertain receive a major-level result;
    # their old Unknown/Ambiguous clusters cannot be reused for exact subtype.
    new_major = major_overlay["major_celltype_refined_v2"].astype(str)
    major_conf = major_overlay["major_confidence_refined_v2"].astype(str)
    boundary = major_conf.eq("hierarchical_boundary")
    exact_unresolved = new_major.isin(["Unknown", "Ambiguous"])
    strict_major = major_overlay["strict_major_eligible_v2"].astype(bool)
    provisional = ~(boundary | exact_unresolved | strict_major)

    out.loc[idx[exact_unresolved], "final_annotation_uncertainty_v2"] = new_major.loc[exact_unresolved].to_numpy()
    out.loc[idx[boundary], "final_annotation_uncertainty_v2"] = new_major.loc[boundary].to_numpy()
    strict_idx = idx[strict_major]
    strict_labels = new_major.loc[strict_major].astype(str) + "_coarse_subtype_pending"
    strict_labels.loc[new_major.loc[strict_major].eq("Epithelial")] = \
        "Epithelial_coarse_CNV_not_recomputed"
    out.loc[strict_idx, "final_annotation_uncertainty_v2"] = strict_labels.to_numpy()
    provisional_idx = idx[provisional]
    out.loc[provisional_idx, "final_annotation_uncertainty_v2"] = (
        new_major.loc[provisional].astype(str) + "_coarse_provisional"
    ).to_numpy()
    out.loc[idx, "annotation_confidence_uncertainty_v2"] = major_conf.to_numpy()
    out.loc[idx, "subtype_resolution_method_v2"] = "major_refined_subtype_requires_recluster"
    out.loc[idx, "strict_subtype_eligible_v2"] = False
    out.loc[idx[~exact_unresolved & ~boundary], "requires_recluster_for_subtype_v2"] = True

    # Expanded-panel rescues require a full panel, >=1 canonical anchor and
    # replication across samples.  These are mature candidate subtype calls.
    rescued = expanded.loc[expanded["rescued"].astype(bool)].copy()
    rescue_label = rescued.set_index("cluster_uid")["proposed_label"].to_dict()
    rescue_marker = rescued.set_index("cluster_uid")["canonical_pass_genes"].to_dict()
    cell_rescue = out["cluster_uid"].astype(str).isin(rescue_label)
    out.loc[cell_rescue, "final_annotation_uncertainty_v2"] = out.loc[
        cell_rescue, "cluster_uid"].astype(str).map(rescue_label).to_numpy()
    out.loc[cell_rescue, "annotation_confidence_uncertainty_v2"] = "medium_expanded_panel"
    out.loc[cell_rescue, "subtype_resolution_method_v2"] = "expanded_panel_canonical_sample_replicated"
    out.loc[cell_rescue, "subtype_descriptive_marker_v2"] = out.loc[
        cell_rescue, "cluster_uid"].astype(str).map(rescue_marker).fillna("").to_numpy()
    out.loc[cell_rescue, "strict_subtype_eligible_v2"] = True

    # Genome-wide FindMarkers suffixes remain explicitly unresolved at subtype
    # level; they make clusters identifiable without overstating biology.
    denovo_map = denovo.set_index("cluster_uid")["proposed_label"].to_dict()
    denovo_marker = denovo.set_index("cluster_uid")["descriptive_marker"].fillna("").to_dict()
    denovo_robust = denovo.set_index("cluster_uid")["top_robust_markers"].fillna("").to_dict()
    denovo_hint = denovo.set_index("cluster_uid")["known_panel_hint"].fillna("").to_dict()
    cell_denovo = out["cluster_uid"].astype(str).isin(denovo_map) & ~cell_rescue
    out.loc[cell_denovo, "final_annotation_uncertainty_v2"] = out.loc[
        cell_denovo, "cluster_uid"].astype(str).map(denovo_map).to_numpy()
    has_marker = out.loc[cell_denovo, "cluster_uid"].astype(str).map(denovo_marker).fillna("").ne("")
    denovo_idx = out.index[cell_denovo]
    out.loc[denovo_idx, "annotation_confidence_uncertainty_v2"] = np.where(
        has_marker.to_numpy(), "descriptive_marker_only", "unresolved"
    )
    out.loc[denovo_idx, "subtype_resolution_method_v2"] = np.where(
        has_marker.to_numpy(), "genomewide_findmarkers_sample_replicated", "no_robust_de_novo_marker"
    )
    out.loc[denovo_idx, "subtype_descriptive_marker_v2"] = out.loc[
        denovo_idx, "cluster_uid"].astype(str).map(denovo_marker).fillna("").to_numpy()
    out.loc[denovo_idx, "subtype_top_robust_markers_v2"] = out.loc[
        denovo_idx, "cluster_uid"].astype(str).map(denovo_robust).fillna("").to_numpy()
    out.loc[denovo_idx, "subtype_known_panel_hint_v2"] = out.loc[
        denovo_idx, "cluster_uid"].astype(str).map(denovo_hint).fillna("").to_numpy()
    out.loc[denovo_idx, "strict_subtype_eligible_v2"] = False

    out["strict_downstream_subtype_eligible_v2"] = (
        out["strict_major_eligible_v2"].astype(bool)
        & out["strict_subtype_eligible_v2"].astype(bool)
        & ~out["exclude_from_primary_major_v2"].astype(bool)
    )
    out.to_pickle(output, compression="gzip")
    (out.groupby(["major_celltype_uncertainty_v2", "annotation_confidence_uncertainty_v2",
                  "subtype_resolution_method_v2"], observed=True)
     .size().rename("n_cells").reset_index().sort_values("n_cells", ascending=False)
     .to_csv(counts, sep="\t", index=False))
    (out.groupby(["final_annotation_v5", "final_annotation_uncertainty_v2"], observed=True)
     .size().rename("n_cells").reset_index().sort_values("n_cells", ascending=False)
     .to_csv(transitions, sep="\t", index=False))

    payload = {
        "run_id": cfg["run_id"], "stage": "04_integrate_uncertainty_and_subtype_overlays",
        "status": "success", "start_time": now_iso(),
        "duration_sec": round(time.time() - started, 3), "n_cells": int(len(out)),
        "major_old_uncertain": int(out["major_celltype_v5"].isin(["Unknown", "Ambiguous"]).sum()),
        "major_new_exact_uncertain": int(out["major_celltype_uncertainty_v2"].isin(["Unknown", "Ambiguous"]).sum()),
        "major_hierarchical_boundary": int(out["major_confidence_uncertainty_v2"].eq("hierarchical_boundary").sum()),
        "major_strictly_resolved_from_old_uncertain": int((out["major_celltype_v5"].isin(["Unknown", "Ambiguous"])
                                                           & out["strict_major_eligible_v2"]).sum()),
        "subtype_cells_rescued_expanded_panel": int(cell_rescue.sum()),
        "subtype_cells_with_de_novo_descriptive_suffix": int((cell_denovo & out[
            "subtype_descriptive_marker_v2"].ne("")).sum()),
        "subtype_cells_still_without_robust_marker": int((cell_denovo & out[
            "subtype_descriptive_marker_v2"].eq("")).sum()),
        "strict_downstream_subtype_eligible": int(out["strict_downstream_subtype_eligible_v2"].sum()),
        "excluded_high_doublet_cross_compartment": int(out["exclude_from_primary_major_v2"].sum()),
        "warnings": [
            "Compact candidate annotation table; authoritative V5 is preserved.",
            "Newly reassigned major classes require lineage-specific reclustering before exact subtype calls.",
            "Newly reassigned epithelial cells have not been rerun through CNV inference.",
            "Labels containing _unresolved_GENE are descriptive cluster identifiers and are excluded from strict subtype analysis.",
        ],
        "outputs": [str(output), str(counts), str(transitions)],
    }
    write_json(summary, payload)
    print(json.dumps(payload, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
