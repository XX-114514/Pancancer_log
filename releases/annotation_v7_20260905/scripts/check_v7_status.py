#!/usr/bin/env python3
from __future__ import annotations
import json
from collections import Counter
from pathlib import Path
import pandas as pd
from common import RUN_ROOT


def status_counts(method: str) -> dict:
    values = []
    for path in (RUN_ROOT / "status" / method).glob("*.json"):
        try:
            values.append(json.loads(path.read_text()).get("status", "unknown"))
        except Exception:
            values.append("unreadable")
    return dict(Counter(values))


def log_tail(name: str, n: int = 4) -> list[str]:
    path = RUN_ROOT / "logs" / name
    if not path.exists():
        return []
    return path.read_text(errors="replace").splitlines()[-n:]


def main() -> None:
    sample_path = RUN_ROOT / "attempts" / "manifest_v7_relaxed_candidate" / "manifests" / "samples.tsv"
    expected = len(pd.read_csv(sample_path, sep="\t")) if sample_path.exists() else 0
    payload = {
        "run_root": str(RUN_ROOT), "expected_samples": expected,
        "infercnv": status_counts("infercnv"), "copykat": status_counts("copykat"),
        "infercnv_log_tail": log_tail("01_infercnv_cohort.log"),
        "copykat_log_tail": log_tail("02_copykat_cohort.log"),
        "integration_complete": (RUN_ROOT / "checkpoints" / "cancer_aware_malignancy_evidence_v7.pkl.gz").exists(),
        "v7_frozen": (RUN_ROOT / "FROZEN_V7.json").exists(),
    }
    print(json.dumps(payload, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

