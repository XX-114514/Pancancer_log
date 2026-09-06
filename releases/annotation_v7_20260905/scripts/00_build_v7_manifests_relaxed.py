#!/usr/bin/env python3
"""Build final V7 roles: hard-exclude doublets from candidates, stricter references."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from common import RUN_ROOT, index_hash, load_config, npz_unicode, open_counts, safe_token


MANIFEST_ROOT = RUN_ROOT / "attempts" / "manifest_v7_relaxed_candidate" / "manifests"


def choose_policy(cancer: str, policies: dict) -> str:
    text = str(cancer).lower()
    for name, policy in policies.items():
        if name != "epithelial_solid_default" and any(pattern.lower() in text for pattern in policy.get("patterns", [])):
            return name
    return "epithelial_solid_default"


def main() -> None:
    cfg = load_config()
    per_sample = MANIFEST_ROOT / "per_sample"
    audit_path = MANIFEST_ROOT / "manifest_audit.json"
    sample_path = MANIFEST_ROOT / "samples.tsv"
    if audit_path.exists():
        audit = json.loads(audit_path.read_text())
        if audit.get("status") == "PASS" and sample_path.exists() and len(list(per_sample.glob("*.npz"))) == cfg["expected_samples"]:
            print(json.dumps({"event": "manifest_reuse", **audit}, ensure_ascii=False), flush=True)
            return
        raise FileExistsError("Existing relaxed V7 manifest is incomplete")
    per_sample.mkdir(parents=True, exist_ok=True)

    ann = pd.read_pickle(cfg["v2_annotation_pickle"])
    v5 = pd.read_pickle(cfg["v5_final_checkpoint"])[
        ["global_position", "sampleID", "gse_id", "cohort", "cancer_type", "qc_annotation_final"]
    ]
    if not ann.index.equals(v5.index) or index_hash(ann.index) != cfg["index_hash_sha256"]:
        raise ValueError("V2/V5 index contract failed")
    _, _, counts_shape = open_counts(cfg)
    if counts_shape != (cfg["expected_cells"], cfg["expected_genes"]):
        raise ValueError(f"Counts shape mismatch: {counts_shape}")

    cells = v5.copy()
    cells.rename(columns={"sampleID": "sample_uid"}, inplace=True)
    cells["major_v2"] = ann["major_celltype_uncertainty_v2"].astype("string").fillna("")
    cells["operational_v2"] = ann["operational_major_candidate_v2"].astype("string").fillna("")
    cells["strict_major_v2"] = ann["strict_major_eligible_v2"].fillna(False).astype(bool)
    cells["artifact_review_v2"] = ann["artifact_or_state_review_v2"].astype("string").fillna("")
    cells.index.name = "cell_id"
    policies = cfg["policies"]
    cancer_to_policy = {cancer: choose_policy(cancer, policies) for cancer in cells["cancer_type"].astype(str).unique()}
    cells["candidate_policy"] = cells["cancer_type"].astype(str).map(cancer_to_policy).astype("category")

    major = cells["major_v2"].astype(str)
    operational = cells["operational_v2"].astype(str)
    uncertain_major = major.isin(set(cfg["uncertain_major_values"])).to_numpy()
    hard_candidate_exclude = cells["artifact_review_v2"].eq("doublet_likely").to_numpy()
    reference_exclude = (
        cells["qc_annotation_final"].astype("string").fillna("").ne("").to_numpy()
        | cells["artifact_review_v2"].isin({"doublet_likely", "ambient_or_doublet_review"}).to_numpy()
    )
    expected = np.zeros(len(cells), dtype=bool)
    context = np.zeros(len(cells), dtype=bool)
    operational_candidate = np.zeros(len(cells), dtype=bool)
    primary = np.zeros(len(cells), dtype=bool)
    fallback = np.zeros(len(cells), dtype=bool)
    for policy_name, positions in cells.groupby("candidate_policy", observed=True).indices.items():
        pos = np.asarray(positions, dtype=np.int64)
        policy = policies[str(policy_name)]
        expected_values = set(policy["candidate_expected"])
        context_values = set(policy["candidate_context"])
        allowed = expected_values | context_values
        expected[pos] = major.iloc[pos].isin(expected_values).to_numpy()
        context[pos] = major.iloc[pos].isin(context_values).to_numpy()
        operational_candidate[pos] = uncertain_major[pos] & operational.iloc[pos].isin(allowed).to_numpy()
        primary[pos] = major.iloc[pos].isin(policy["reference_primary"]).to_numpy()

    candidate = (expected | context | operational_candidate) & ~hard_candidate_exclude
    reference_clean = cells["strict_major_v2"].to_numpy() & ~reference_exclude
    primary &= ~candidate & reference_clean
    groups = cells.groupby("sample_uid", sort=True, observed=True).indices
    for positions in groups.values():
        pos = np.asarray(positions, dtype=np.int64)
        if int(primary[pos].sum()) < int(cfg["min_reference_cells"]):
            policy = policies[str(cells.iloc[pos[0]]["candidate_policy"])]
            fallback[pos] = major.iloc[pos].isin(policy["reference_fallback"]).to_numpy() & ~candidate[pos] & reference_clean[pos]
    reference = primary | fallback

    source = np.full(len(cells), "not_candidate", dtype="U34")
    source[expected] = "expected_major_v2"
    source[context] = "context_major_v2"
    source[operational_candidate] = "uncertain_operational_v2"
    source[hard_candidate_exclude & (expected | context | operational_candidate)] = "excluded_doublet_likely"
    reference_source = np.full(len(cells), "not_reference", dtype="U24")
    reference_source[primary] = "strict_primary_v2"
    reference_source[fallback] = "strict_fallback_v2"

    rows = []
    for ordinal, (sample_uid, positions) in enumerate(groups.items(), start=1):
        pos = np.asarray(positions, dtype=np.int64)
        token = safe_token(str(sample_uid))
        npz_path = per_sample / f"{token}.npz"
        np.savez_compressed(
            npz_path, global_position=pos, cell_id=npz_unicode(cells.index[pos]),
            tumor_lineage_candidate_v2=candidate[pos], tumor_reference_cell_v2=reference[pos],
            tumor_lineage_candidate_v7=candidate[pos], tumor_reference_cell_v7=reference[pos],
            candidate_source_v7=npz_unicode(source[pos]), reference_source_v7=npz_unicode(reference_source[pos]),
        )
        first = int(pos[0])
        n_candidate = int(candidate[pos].sum())
        n_reference = int(reference[pos].sum())
        rows.append({
            "sample_uid": str(sample_uid), "sample_token": token,
            "gse_id": str(cells.iloc[first]["gse_id"]), "cohort": str(cells.iloc[first]["cohort"]),
            "cancer_type": str(cells.iloc[first]["cancer_type"]), "candidate_policy": str(cells.iloc[first]["candidate_policy"]),
            "n_cells": int(len(pos)), "n_candidate": n_candidate,
            "n_candidate_expected": int((expected[pos] & ~hard_candidate_exclude[pos]).sum()),
            "n_candidate_context": int((context[pos] & ~hard_candidate_exclude[pos]).sum()),
            "n_candidate_operational": int((operational_candidate[pos] & ~hard_candidate_exclude[pos]).sum()),
            "n_reference": n_reference, "n_reference_primary": int(primary[pos].sum()),
            "n_reference_fallback": int(fallback[pos].sum()),
            "n_candidate_excluded_doublet": int((hard_candidate_exclude[pos] & (expected[pos] | context[pos] | operational_candidate[pos])).sum()),
            "n_reference_excluded_review": int(reference_exclude[pos].sum()),
            "infercnv_eligible": n_candidate >= int(cfg["min_candidate_cells"]) and n_reference >= int(cfg["min_reference_cells"]),
            "copykat_eligible": n_candidate >= int(cfg["min_copykat_cells"]),
            "positions_contiguous": bool(np.all(np.diff(pos) == 1)) if len(pos) > 1 else True,
            "manifest_npz": str(npz_path.resolve()),
        })
        if ordinal % 100 == 0:
            print(f"relaxed manifest heartbeat {ordinal}/{len(groups)}", flush=True)

    samples = pd.DataFrame(rows).sort_values(["candidate_policy", "cancer_type", "gse_id", "sample_uid"])
    if len(samples) != cfg["expected_samples"] or int(samples["n_cells"].sum()) != cfg["expected_cells"]:
        raise ValueError("Sample universe contract failed")
    samples.to_csv(sample_path, sep="\t", index=False)
    policy_summary = samples.groupby("candidate_policy", observed=True).agg(
        n_samples=("sample_uid", "size"), n_cells=("n_cells", "sum"), n_candidate=("n_candidate", "sum"),
        n_candidate_operational=("n_candidate_operational", "sum"), n_reference=("n_reference", "sum"),
        infercnv_eligible_samples=("infercnv_eligible", "sum"), copykat_eligible_samples=("copykat_eligible", "sum"),
        fallback_reference_samples=("n_reference_fallback", lambda values: int((values > 0).sum())),
    ).reset_index()
    policy_summary.to_csv(MANIFEST_ROOT / "policy_summary.tsv", sep="\t", index=False)
    low_candidate = samples["n_candidate"] < int(cfg["min_candidate_cells"])
    low_reference = samples["n_reference"] < int(cfg["min_reference_cells"])
    preflight_skipped = samples.loc[low_candidate | low_reference].copy()
    preflight_skipped["skip_reason"] = np.select(
        [(low_candidate & low_reference)[low_candidate | low_reference], low_candidate[low_candidate | low_reference]],
        ["insufficient_candidate_and_reference_cells", "insufficient_candidate_cells"],
        default="insufficient_reference_cells",
    )
    preflight_skipped.to_csv(MANIFEST_ROOT / "preflight_skipped_samples.tsv", sep="\t", index=False)
    audit = {
        "status": "PASS", "run_id": cfg["run_id"], "role_definition": "relaxed_candidate_strict_reference",
        "annotation_parent": "uncertainty_resolution_v2", "n_cells": int(len(cells)),
        "n_samples": int(len(samples)), "n_gse": int(samples["gse_id"].nunique()),
        "counts_shape": list(counts_shape), "index_hash_sha256": index_hash(cells.index),
        "candidate_cells": int(candidate.sum()), "reference_cells": int(reference.sum()),
        "candidate_operational_cells": int((operational_candidate & ~hard_candidate_exclude).sum()),
        "candidate_doublets_excluded": int((hard_candidate_exclude & (expected | context | operational_candidate)).sum()),
        "infercnv_eligible_samples": int(samples["infercnv_eligible"].sum()),
        "copykat_eligible_samples": int(samples["copykat_eligible"].sum()),
        "noncontiguous_samples": int((~samples["positions_contiguous"]).sum()),
        "cancer_to_policy": cancer_to_policy,
    }
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(audit, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

