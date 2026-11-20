#!/usr/bin/env python3
"""Index known PanCancer result registries without traversing native data.

Detailed paths and status receipts are server-only outputs. Directory counts
are workflow records, never study, sample, patient, or independent-finding counts.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

STATUS_FILES = {
    "complete.json", "independent_audit.json", "acceptance.json",
    "receipt.json", "verification.json", "input_summary.json",
    "upstream_audit_linkage_current.json", "independent_readback.json",
    "qualification.json", "independent_replay.json",
}
REGISTRIES = ["runs", "annotation_rebuild", "reports", "manuscripts", "run_catalog",
              "results", "CoVarNet", "Downstream", "Spatial"]
REPORT_DIRS = ["outputs", "paper", "review", "final_delivery_v1", "checkpoints"]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def utc(stamp):
    return datetime.fromtimestamp(stamp, timezone.utc).isoformat()


def classify(name):
    tags = []
    for tag, pattern in [
        ("annotation", r"annot|identity|lineage|celltypist|reference"),
        ("composition_module", r"covarnet|cm_|nmf|membership"),
        ("RNA_state", r"state|ucell|program|dialogue|mofa"),
        ("spatial", r"spatial|visium|xenium|cosmx|imaging|ccrcc|histology"),
        ("clinical", r"clinical|surviv|bulk|response|ici_|icb_|tcga|lihc"),
        ("communication", r"communication|cpdb|liana|cellchat|_lr_|receiver"),
        ("audit_or_delivery", r"audit|accept|deliver|handoff|review|archive"),
    ]:
        if re.search(pattern, name.lower()):
            tags.append(tag)
    return tags or ["other_or_legacy"]


def index(project_root, handoff_root, output):
    if output.exists():
        raise FileExistsError("Use a new versioned output directory")
    records = []
    seen = set()

    def add(path, role, registry):
        key = str(path.absolute())
        if key in seen or not path.is_dir():
            return
        seen.add(key)
        files = [p for p in path.iterdir() if p.is_file() and not p.name.startswith(".")]
        reports = []
        for folder in [path] + [path / s for s in REPORT_DIRS]:
            if not folder.is_dir():
                continue
            candidates = [p for p in folder.iterdir() if p.is_file() and (
                p.name in STATUS_FILES or (p.suffix == ".json" and (
                    p.name.startswith("accepted_") or "summary" in p.name.lower()
                    or "audit" in p.name.lower() or "receipt" in p.name.lower())))]
            for p in sorted(candidates):
                row = {"relative_path": str(p.relative_to(path)), "bytes": p.stat().st_size,
                       "mtime_UTC": utc(p.stat().st_mtime)}
                if p.stat().st_size <= 2097152:
                    row["sha256"] = sha(p)
                else:
                    row["fingerprint_scope"] = "metadata_only_over_2MiB"
                if p.stat().st_size <= 131072:
                    try:
                        item = json.loads(p.read_text())
                        if isinstance(item, dict):
                            row["indicated_status"] = {
                                k: item[k] for k in ["status", "stage", "scope", "schema"]
                                if isinstance(item.get(k), (str, bool, int))
                            }
                    except (UnicodeError, ValueError):
                        row["parse_status"] = "not_readable_JSON"
                reports.append(row)
        records.append({
            "root_role": role, "registry": registry, "path": key, "run_id": path.name,
            "directory_mtime_UTC": utc(path.stat().st_mtime),
            "tags_from_name_only": classify(path.name),
            "top_level_file_count": len(files),
            "top_level_directory_count": sum(p.is_dir() for p in path.iterdir()),
            "bounded_status_receipts": reports,
        })

    for i in range(1, 6):
        project = project_root / f"Project_v{i}"
        if not project.is_dir():
            continue
        for child in sorted(project.iterdir()):
            if child.is_dir() and not child.name.startswith("."):
                add(child, f"Project_v{i}", "top_level")
        for name in REGISTRIES:
            registry = project / name
            if registry.is_dir():
                for child in sorted(registry.iterdir()):
                    if child.is_dir() and not child.name.startswith("."):
                        add(child, f"Project_v{i}", name)
    for child in sorted(handoff_root.iterdir()):
        if child.is_dir() and not child.name.startswith("."):
            add(child, "WEHI_handoff", "top_level")
    validation = handoff_root / "CM_validation_20261118_v1"
    if validation.is_dir():
        for child in sorted(validation.iterdir()):
            if child.is_dir():
                add(child, "R4_validation", "validation_top_level")
    output.mkdir(parents=True)
    index_path = output / "server_result_registry.json"
    index_path.write_text(json.dumps(records, ensure_ascii=False, indent=2) + "\n")
    receipt = {
        "status": "BOUNDED_RESULT_REGISTRY_INDEX_COMPLETE_NOT_ALL_RAW_RESULTS_REAUDITED",
        "UTC": datetime.now(timezone.utc).isoformat(),
        "registry_records": len(records),
        "records_by_root_role": dict(Counter(r["root_role"] for r in records)),
        "status_receipts_indexed": sum(len(r["bounded_status_receipts"]) for r in records),
        "registry_sha256": sha(index_path), "script_sha256": sha(Path(__file__)),
        "scope": "Known Project_v1..5 registries and WEHI/R4 result roots; bounded metadata and status JSON only. No native data traversal, model arrays, per-cell or per-patient rows read.",
        "limitations": "Name tags are navigation, not scientific classification. Filesystem mtime is not analysis completion time. Old pending statuses require later authoritative handoff review. Directory records overlap scientific inputs and do not represent independent experiments.",
    }
    (output / "inventory_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project-root", type=Path, required=True)
    p.add_argument("--handoff-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    index(a.project_root, a.handoff_root, a.output)
