#!/usr/bin/env python3
"""Map V7 per-sample CNV results and add exact role/status provenance."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

import common

MANIFEST_ROOT = common.RUN_ROOT / "attempts" / "manifest_v7" / "manifests"
INTERMEDIATE = common.RUN_ROOT / "checkpoints" / "cancer_aware_malignancy_evidence.pkl.gz"
OUTPUT = common.RUN_ROOT / "checkpoints" / "cancer_aware_malignancy_evidence_v7.pkl.gz"

RENAME = {
    "candidate_policy_v2": "candidate_policy_v7", "cnv_role_v2": "cnv_role_v7",
    "tumor_lineage_candidate_v2": "tumor_lineage_candidate_v7",
    "tumor_reference_cell_v2": "tumor_reference_cell_v7",
    "sample_infercnv_evaluable_v2": "sample_infercnv_evaluable_v7",
    "sample_infercnv_positive_v2": "sample_infercnv_positive_v7",
    "sample_cnv_burden_v2": "sample_cnv_burden_v7",
    "sample_copykat_evaluable_v2": "sample_copykat_evaluable_v7",
    "sample_copykat_aneuploid_v2": "sample_copykat_aneuploid_v7",
    "malignant_intersection_v2": "malignant_intersection_v7",
    "malignant_union_v2": "malignant_union_v7",
    "sample_cnv_cluster_v2": "sample_cnv_cluster_v7",
    "sample_cluster_label_v2": "sample_cluster_label_v7",
    "malignant_consensus_v2": "malignant_consensus_v7",
    "malignancy_evidence_version_v2": "malignancy_evidence_version_v7",
}


def read_status(method: str, token: str) -> dict:
    path = common.RUN_ROOT / "status" / method / f"{token}.json"
    if not path.exists():
        return {"status": "missing", "error": "status_file_missing"}
    return json.loads(path.read_text())


def explain(method: str, status: dict) -> str:
    value = str(status.get("status", "missing"))
    if value == "success":
        return "success"
    if value == "skipped_insufficient_cells":
        return f"{method}: candidate/reference below method threshold"
    if value == "skipped_insufficient_quality_cells":
        return f"{method}: fewer than 200 candidates with >200 detected genes"
    if value == "skipped_algorithm_data_limitation":
        limitation = status.get("limitation") or {}
        return f"{method}: algorithm data limitation: {limitation.get('code', 'unspecified')}"
    return f"{method}: {value}: {status.get('error', '')}".rstrip(": ")


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite versioned output: {OUTPUT}")
    cfg = common.load_config()
    base = Path(cfg["legacy_engine_run"]) / "scripts" / "03_integrate_cancer_aware_evidence.py"
    spec = importlib.util.spec_from_file_location("validated_v7_evidence_integrator", base)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {base}")
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    engine.MANIFEST_ROOT = MANIFEST_ROOT
    engine.main()

    evidence = pd.read_pickle(INTERMEDIATE).rename(columns=RENAME)
    samples = pd.read_csv(MANIFEST_ROOT / "samples.tsv", sep="\t")
    index_values = evidence.index.to_numpy(dtype=str)
    n_cells = len(evidence)
    candidate_source = np.full(n_cells, "", dtype=object)
    reference_source = np.full(n_cells, "", dtype=object)
    infer_status = np.full(n_cells, "", dtype=object)
    copy_status = np.full(n_cells, "", dtype=object)
    skip_reason = np.full(n_cells, "", dtype=object)
    skipped_rows = []
    for row in samples.itertuples(index=False):
        manifest = np.load(row.manifest_npz, allow_pickle=False)
        positions = manifest["global_position"].astype(np.int64)
        ids = manifest["cell_id"].astype(str)
        if not np.array_equal(index_values[positions], ids):
            raise ValueError(f"V7 mapping mismatch for {row.sample_uid}")
        candidate_source[positions] = manifest["candidate_source_v7"].astype(str)
        reference_source[positions] = manifest["reference_source_v7"].astype(str)
        inf = read_status("infercnv", row.sample_token)
        copy = read_status("copykat", row.sample_token)
        inf_value = str(inf.get("status", "missing"))
        copy_value = str(copy.get("status", "missing"))
        infer_status[positions] = inf_value
        copy_status[positions] = copy_value
        non_success = [item for item in (explain("inferCNVpy", inf), explain("Copykat_python", copy)) if item != "success"]
        reason = "; ".join(non_success)
        skip_reason[positions] = reason
        if non_success:
            skipped_rows.append({
                "sample_uid": row.sample_uid, "sample_token": row.sample_token,
                "gse_id": row.gse_id, "cancer_type": row.cancer_type,
                "candidate_policy": row.candidate_policy, "n_cells": int(row.n_cells),
                "n_candidate": int(row.n_candidate), "n_reference": int(row.n_reference),
                "infercnv_status": inf_value, "copykat_status": copy_value, "explanation": reason,
            })

    evidence["candidate_source_v7"] = pd.Categorical(candidate_source)
    evidence["reference_source_v7"] = pd.Categorical(reference_source)
    evidence["infercnv_run_status_v7"] = pd.Categorical(infer_status)
    evidence["copykat_run_status_v7"] = pd.Categorical(copy_status)
    evidence["cnv_skip_reason_v7"] = pd.Categorical(skip_reason)
    evidence["copykat_engine_v7"] = pd.Categorical.from_codes(
        np.zeros(n_cells, dtype=np.int8), ["Copykat_python_1.0.0_c2390e8a46e2"]
    )
    evidence["malignancy_evidence_version_v7"] = pd.Categorical.from_codes(
        np.zeros(n_cells, dtype=np.int8), ["v7_per_sample_counts_annotation_uncertainty_v2_20260905"]
    )
    evidence.to_pickle(OUTPUT, compression="gzip")

    table_dir = common.RUN_ROOT / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    samples.to_csv(table_dir / "per_sample_cnv_targets.tsv", sep="\t", index=False)
    pd.DataFrame(skipped_rows).to_csv(table_dir / "skipped_samples_with_explanations.tsv", sep="\t", index=False)
    pd.crosstab(evidence["infercnv_run_status_v7"].astype(str), evidence["copykat_run_status_v7"].astype(str)).to_csv(
        table_dir / "method_status_cell_counts_v7.tsv", sep="\t"
    )
    summary = {
        "status": "PASS", "n_cells": n_cells, "n_samples": int(len(samples)),
        "candidate_cells": int(evidence["tumor_lineage_candidate_v7"].sum()),
        "reference_cells": int(evidence["tumor_reference_cell_v7"].sum()),
        "skipped_or_partial_samples": int(len(skipped_rows)),
        "mapping_conflicts": 0, "unmapped_cells": 0, "output": str(OUTPUT.resolve()),
    }
    report_dir = common.RUN_ROOT / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    (report_dir / "03_v7_integration_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
