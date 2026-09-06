#!/usr/bin/env python3
"""Freeze V2 identity plus rerun CNV evidence into the immutable V7 schema."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from common import RUN_ROOT, index_hash, load_config


EVIDENCE = RUN_ROOT / "checkpoints" / "cancer_aware_malignancy_evidence_v7.pkl.gz"
CHECKPOINT = RUN_ROOT / "checkpoints" / "final_annotation_v7.pkl.gz"
OBS_ONLY = RUN_ROOT / "objects" / "final_annotation_v7_obs_only.h5ad"
TABLES = RUN_ROOT / "tables"
REPORTS = RUN_ROOT / "reports"
FROZEN = RUN_ROOT / "FROZEN_V7.json"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(8 * 1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    cfg = load_config()
    for path in (CHECKPOINT, OBS_ONLY, FROZEN):
        if path.exists():
            raise FileExistsError(f"Refusing to overwrite frozen/versioned output: {path}")
    TABLES.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    OBS_ONLY.parent.mkdir(parents=True, exist_ok=True)

    v5 = pd.read_pickle(cfg["v5_final_checkpoint"])
    v2 = pd.read_pickle(cfg["v2_annotation_pickle"])
    evidence = pd.read_pickle(EVIDENCE)
    if not (v5.index.equals(v2.index) and v5.index.equals(evidence.index) and v5.index.is_unique):
        raise ValueError("V5/V2/V7-evidence cell indices are not exact and unique")
    if len(v5) != cfg["expected_cells"] or index_hash(v5.index) != cfg["index_hash_sha256"]:
        raise ValueError("V7 cell universe/hash contract failed")
    if not v5["sampleID"].astype(str).equals(evidence["sample_uid"].astype(str)):
        raise ValueError("sampleID/sample_uid mismatch")

    final = v5.copy()
    for column in v2.columns:
        if column in {"gse_id", "sampleID", "cancer_type", "cohort", "cluster_uid"}:
            if not v2[column].astype(str).equals(final[column].astype(str)):
                raise ValueError(f"V2 metadata mismatch: {column}")
            continue
        if column not in final.columns:
            final[column] = v2[column]
    for column in evidence.columns:
        if column in {"global_position", "gse_id", "cohort", "cancer_type", "sample_uid"}:
            continue
        if column in final.columns:
            raise ValueError(f"Unexpected V7 evidence collision: {column}")
        final[column] = evidence[column]

    final["major_celltype_v7"] = pd.Categorical(v2["major_celltype_uncertainty_v2"].astype(str))
    final["final_annotation_v7"] = pd.Categorical(v2["final_annotation_uncertainty_v2"].astype(str))
    final["annotation_confidence_v7"] = pd.Categorical(v2["annotation_confidence_uncertainty_v2"].astype(str))
    final["strict_major_eligible_v7"] = v2["strict_major_eligible_v2"].fillna(False).astype(bool)
    final["strict_subtype_eligible_v7"] = v2["strict_subtype_eligible_v2"].fillna(False).astype(bool)
    final["strict_downstream_subtype_eligible_v7"] = v2["strict_downstream_subtype_eligible_v2"].fillna(False).astype(bool)

    raw_label = final["sample_cluster_label_v7"].astype(str)
    label = raw_label.replace({"non_malignant": "cnv_no_malignancy_detected"})
    label_order = ["not_candidate", "not_evaluable", "cnv_no_malignancy_detected", "uncertain", "malignant"]
    final["malignancy_call_v7"] = pd.Categorical(label, categories=label_order)
    candidate = final["tumor_lineage_candidate_v7"].astype(bool).to_numpy()
    intersection = final["malignant_intersection_v7"].fillna(False).astype(bool).to_numpy()
    union = final["malignant_union_v7"].fillna(False).astype(bool).to_numpy()
    malignant = label.eq("malignant").to_numpy()
    uncertain = label.eq("uncertain").to_numpy()
    no_detection = label.eq("cnv_no_malignancy_detected").to_numpy()
    not_evaluable = label.eq("not_evaluable").to_numpy()

    tier = np.full(len(final), "not_candidate", dtype=object)
    tier[candidate] = "not_evaluable"
    tier[no_detection] = "cnv_no_malignancy_detected"
    tier[uncertain] = "uncertain"
    tier[union] = "sensitive_union_only"
    tier[malignant] = "cluster_supported_malignant"
    tier[intersection] = "high_confidence_malignant"
    tier_order = [
        "not_candidate", "not_evaluable", "cnv_no_malignancy_detected", "uncertain",
        "sensitive_union_only", "cluster_supported_malignant", "high_confidence_malignant",
    ]
    final["malignancy_tier_v7"] = pd.Categorical(tier, categories=tier_order, ordered=True)
    final["is_malignant_primary_v7"] = pd.arrays.BooleanArray(
        malignant, ~(malignant | no_detection)
    )
    final["is_malignant_high_confidence_v7"] = final["malignant_intersection_v7"].astype("boolean")
    final["is_malignant_sensitive_v7"] = final["malignant_union_v7"].astype("boolean")
    final["strict_malignant_downstream_eligible_v7"] = intersection

    identity = final["final_annotation_v7"].astype(str).to_numpy(dtype=object)
    display = identity.copy()
    display[malignant] = np.char.add("Malignant::", identity[malignant].astype(str))
    display[uncertain] = np.char.add("Malignancy_uncertain::", identity[uncertain].astype(str))
    display[candidate & not_evaluable] = np.char.add(
        "Malignancy_not_evaluable::", identity[candidate & not_evaluable].astype(str)
    )
    final["final_annotation_with_malignancy_v7"] = pd.Categorical(display)
    final["annotation_schema_version_v7"] = pd.Categorical.from_codes(
        np.zeros(len(final), dtype=np.int8), [cfg["annotation_schema"]]
    )
    final["freeze_version"] = pd.Categorical.from_codes(np.zeros(len(final), dtype=np.int8), ["V7"])
    final["freeze_date"] = pd.Categorical.from_codes(np.zeros(len(final), dtype=np.int8), ["2026-09-05"])

    final["malignancy_tier_v7"].value_counts(dropna=False).rename_axis("malignancy_tier_v7").rename("n_cells").to_csv(
        TABLES / "malignancy_tier_v7_counts.tsv", sep="\t"
    )
    final["major_celltype_v7"].value_counts(dropna=False).rename_axis("major_celltype_v7").rename("n_cells").to_csv(
        TABLES / "major_celltype_v7_counts.tsv", sep="\t"
    )
    final["final_annotation_v7"].value_counts(dropna=False).rename_axis("final_annotation_v7").rename("n_cells").to_csv(
        TABLES / "final_annotation_v7_counts.tsv", sep="\t"
    )
    pd.crosstab(final["cancer_type"].astype(str), final["malignancy_tier_v7"].astype(str)).to_csv(
        TABLES / "malignancy_tier_v7_by_cancer.tsv", sep="\t"
    )
    pd.crosstab(final["major_celltype_final"].astype(str), final["major_celltype_v7"].astype(str)).to_csv(
        TABLES / "major_v5_to_v7_transition.tsv", sep="\t"
    )
    pd.crosstab(final["malignancy_call_v3"].astype(str), final["malignancy_call_v7"].astype(str)).to_csv(
        TABLES / "malignancy_v3_to_v7_transition.tsv", sep="\t"
    )
    changed = final["malignancy_call_v3"].astype(str).ne(final["malignancy_call_v7"].astype(str))
    pd.crosstab(final.loc[changed, "cancer_type"].astype(str), final.loc[changed, "malignancy_call_v7"].astype(str)).to_csv(
        TABLES / "changed_malignancy_cells_by_cancer.tsv", sep="\t"
    )

    temp_pickle = CHECKPOINT.with_name(f".{CHECKPOINT.name}.{os.getpid()}.tmp")
    temp_h5ad = OBS_ONLY.with_name(f".{OBS_ONLY.name}.{os.getpid()}.tmp.h5ad")
    final.to_pickle(temp_pickle, compression="gzip")
    obs_columns = [
        "gse_id", "sampleID", "cancer_type", "cohort", "cluster_uid",
        "major_celltype_v7", "final_annotation_v7", "annotation_confidence_v7",
        "strict_major_eligible_v7", "strict_subtype_eligible_v7",
        "strict_downstream_subtype_eligible_v7", "candidate_policy_v7", "cnv_role_v7",
        "candidate_source_v7", "reference_source_v7", "sample_infercnv_evaluable_v7",
        "sample_infercnv_positive_v7", "sample_cnv_burden_v7", "sample_copykat_evaluable_v7",
        "sample_copykat_aneuploid_v7", "malignant_intersection_v7", "malignant_union_v7",
        "sample_cnv_cluster_v7", "sample_cluster_label_v7", "malignant_consensus_v7",
        "malignancy_call_v7", "malignancy_tier_v7", "is_malignant_primary_v7",
        "is_malignant_high_confidence_v7", "is_malignant_sensitive_v7",
        "strict_malignant_downstream_eligible_v7", "final_annotation_with_malignancy_v7",
        "infercnv_run_status_v7", "copykat_run_status_v7", "cnv_skip_reason_v7",
        "copykat_engine_v7", "malignancy_evidence_version_v7", "annotation_schema_version_v7",
        "freeze_version", "freeze_date",
    ]
    ad.settings.allow_write_nullable_strings = True
    obs = ad.AnnData(obs=final[obs_columns].copy())
    obs.strings_to_categoricals()
    obs.uns.update({
        "annotation_schema_version": cfg["annotation_schema"],
        "identity_parent": cfg["v2_annotation_pickle"],
        "v5_parent": cfg["v5_final_checkpoint"],
        "malignancy_evidence_parent": str(EVIDENCE.resolve()),
        "counts_zarr": cfg["counts_zarr"], "index_hash_sha256": cfg["index_hash_sha256"],
        "primary_malignancy_rule": "sample_cluster_label_v7 == malignant",
        "strict_malignancy_rule": "inferCNVpy positive AND Copykat_python aneuploid",
        "sensitive_malignancy_rule": "inferCNVpy positive OR Copykat_python aneuploid",
        "negative_call_note": "cnv_no_malignancy_detected is not proof of benign identity",
    })
    obs.write_h5ad(temp_h5ad, compression="gzip", compression_opts=4)
    del obs
    check = ad.read_h5ad(temp_h5ad, backed="r")
    required = {"major_celltype_v7", "final_annotation_v7", "malignancy_tier_v7", "final_annotation_with_malignancy_v7"}
    if check.shape != (len(final), 0) or not check.obs_names.equals(final.index) or not required.issubset(check.obs.columns):
        raise ValueError("Written V7 obs-only AnnData validation failed")
    check.file.close()
    os.replace(temp_pickle, CHECKPOINT)
    os.replace(temp_h5ad, OBS_ONLY)

    summary = {
        "status": "FROZEN", "freeze_version": "V7", "run_id": cfg["run_id"],
        "annotation_schema": cfg["annotation_schema"], "n_cells": int(len(final)),
        "n_samples": int(final["sampleID"].nunique()), "n_gse": int(final["gse_id"].nunique()),
        "n_cancers": int(final["cancer_type"].nunique()), "index_hash_sha256": index_hash(final.index),
        "candidate_cells": int(candidate.sum()), "reference_cells": int(final["tumor_reference_cell_v7"].sum()),
        "primary_malignant_cells": int(malignant.sum()), "uncertain_cells": int(uncertain.sum()),
        "high_confidence_malignant_cells": int(intersection.sum()), "sensitive_union_cells": int(union.sum()),
        "not_evaluable_candidate_cells": int((candidate & not_evaluable).sum()),
        "strict_major_cells": int(final["strict_major_eligible_v7"].sum()),
        "strict_subtype_cells": int(final["strict_subtype_eligible_v7"].sum()),
        "changed_malignancy_call_vs_v3": int(changed.sum()),
        "checkpoint": str(CHECKPOINT.resolve()), "obs_only_h5ad": str(OBS_ONLY.resolve()),
        "checkpoint_sha256": sha256_file(CHECKPOINT), "obs_only_h5ad_sha256": sha256_file(OBS_ONLY),
    }
    (REPORTS / "04_v7_freeze_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    REPORTS.joinpath("V7_ANNOTATION_AND_CNV_LOGIC_CN.md").write_text(f"""# V7 注释与 CNV 冻结逻辑

- 细胞宇宙：{len(final):,} 个细胞，索引哈希 `{summary['index_hash_sha256']}`。
- 身份主轴：V2 保守大类与亚型；V5 全部字段保留为历史与审计父版本。
- CNV 输入：5M 合并对象中的逐样本原始非负整数 counts；不使用 Harmony/scVI/整合表达矩阵。
- 候选：按癌种指定肿瘤谱系；Unknown/Ambiguous/边界细胞仅在 operational 候选落入该癌种肿瘤谱系时纳入。
- 参考：V2 strict major，且排除候选、QC、doublet_likely 与 ambient_or_doublet_review。
- 严格恶性：inferCNVpy 阳性与固定版本 Copykat_python 非整倍体的交集。
- 主恶性：逐样本 CNV cluster 规则；并集只作为敏感证据层。
- `cnv_no_malignancy_detected` 仅表示未检出 CNV 恶性证据，不等于证明良性。
- CopyKAT 身份：本次使用 `Copykat_python 1.0.0 @ c2390e8a46e2`，不是 R CopyKAT。
- 跳过样本：见 `tables/skipped_samples_with_explanations.tsv`，绝不把跳过强制改为正常。
""")
    FROZEN.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()

