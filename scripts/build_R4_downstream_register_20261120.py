#!/usr/bin/env python3
"""Assemble existing R4/K15 aggregate evidence; never refit or infer new P values.

Run on the research server. Inputs and scientific outputs stay there. The
module rows count events/settings, not independent patients or discoveries.
"""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

MODULES = {f"CM{i:02d}" for i in range(1, 16)}


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def check_hash(path, expected):
    if digest(path) != expected:
        raise ValueError(f"Changed input: {path.name}")


def check_modules(rows):
    ids = [r["cm"] for r in rows]
    if len(ids) != 15 or set(ids) != MODULES:
        raise ValueError("Expected exactly one aggregate row per frozen CM01..15")


def finite(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None


def write_table(path, rows):
    with path.open("x", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def build(root, output):
    provenance = {}

    def read_json(rel):
        p = root / rel
        provenance[rel] = digest(p)
        return json.loads(p.read_text())

    def read_table(rel):
        p = root / rel
        provenance[rel] = digest(p)
        with p.open() as f:
            return list(csv.DictReader(f, delimiter="\t"))

    cp = "communication_joint15_v2/"
    hl = "highlow_joint15_v1/"
    sp = "spatial_melanoma/N_sensitivity_all24_v1/"
    ca = read_json(cp + "independent_audit.json")
    ha = read_json(hl + "independent_audit.json")
    link = read_json(hl + "upstream_audit_linkage_current.json")
    sc = read_json(sp + "complete.json")
    sr = read_json(sp + "independent_readback.json")
    for audit in [ca, ha, sc, sr]:
        if not audit["status"].startswith("PASS"):
            raise ValueError("Required independent audit not passing")
    if link["status"] != "ALL15_UPSTREAM_AUDITS_LINKED" or link["modules_pass"] != 15:
        raise ValueError("Current upstream audit linkage incomplete")
    check_modules(link["records"])
    checked_files = 0
    for row in link["records"]:
        rel = f"highlow/{row['cm']}/independent_audit.json"
        check_hash(root / rel, row["audit_sha256"])
        audit = read_json(rel)
        if not audit["status"].startswith("PASS"):
            raise ValueError("Upstream audit not passing")
        for name, sha in audit["hashes"].items():
            check_hash(root / "highlow" / row["cm"] / name, sha)
            checked_files += 1
    cc = read_json(cp + "complete.json")
    for name, sha in cc["output_sha256"].items():
        check_hash(root / cp / name, sha)
        checked_files += 1
    for name, sha in sr["hashes"].items():
        check_hash(root / sp / name, sha)
        checked_files += 1
    for name, sha in sc["source_hashes"].items():
        check_hash(root / "spatial_melanoma" / name, sha)
        checked_files += 1
    check_hash(root / cp / "module_summary.tsv", ca["summary_sha256"])
    check_hash(root / hl / "comparison_summary.tsv", ha["summary_sha256"])

    communication = read_table(cp + "module_summary.tsv")
    check_modules(communication)
    expression = read_table(hl + "comparison_summary.tsv")
    edges = read_table(sp + "all_edge_metrics.tsv")
    scores = read_table(sp + "all_author_scores.tsv")
    sensitivity = read_table(sp + "all_module_sensitivity_summary.tsv")
    if len(edges) != sc["edge_rows"] or len(scores) != sc["author_rows"]:
        raise ValueError("Spatial row count differs from audited receipt")
    if len(scores) != 8 * 3 * 15 * 2 * 2:
        raise ValueError("Expected all eight sections, three N settings and 15 modules")
    score_keys = [(r["sample"],r["assumed_N"],r["cm"],r["posterior"],r["statistic"]) for r in scores]
    if len(score_keys) != len(set(score_keys)):
        raise ValueError("Duplicate author-score keys")
    if sum(int(r["tests"]) for r in communication) != ca["tests"]:
        raise ValueError("Communication denominator mismatch")
    if sum(int(r["joint15_events"]) for r in communication) != ca["joint_significant_events"]:
        raise ValueError("Communication event total mismatch")
    if sum(int(r["n_tests"]) for r in expression) != ha["n_tests"]:
        raise ValueError("Expression test total mismatch")

    module_rows = []
    for cm in sorted(MODULES):
        c = next(r for r in communication if r["cm"] == cm)
        x = [r for r in expression if r["cm"] == cm and r["family"] == "unit"]
        members = {r["node"] for r in expression if r["cm"] == cm}
        edge_set = {(r["node1"], r["node2"]) for r in edges if r["cm"] == cm}
        module_rows.append({
            "cm_version": "2026-11-15_R4_K15", "cm": cm,
            "original_network_members": len(members),
            "member_labels": ";".join(sorted(members)),
            "network_edges": len(edge_set),
            "CPDB_joint_tests": int(c["tests"]),
            "CPDB_joint_q05_events_not_independent_discoveries": int(c["joint15_events"]),
            "unit_DE_eligible_member_comparisons": len(x),
            "unit_DE_ineligible_member_comparisons": len(members) - len(x),
            "unit_DE_tests": sum(int(r["n_tests"]) for r in x),
            "unit_DE_model_family_q05_gene_identity_events": sum(int(r["q_model05"]) for r in x),
            "audit_state": "JOINT_AND_CURRENT_UPSTREAM_AUDITS_HASH_LINKED",
            "interpretation": "discovery_data;units_not_all_verified_patients;no_module_specific_LR_claim",
        })

    grouped = defaultdict(list)
    for r in edges:
        if r["cm"] not in MODULES or r["assumed_N"] not in {"5", "30", "60"}:
            raise ValueError("Unexpected CM or N")
        grouped[(r["cm"], r["posterior"], r["assumed_N"])].append(r)
    spatial_rows = []
    for (cm, posterior, n), rows in sorted(grouped.items()):
        keys = [(r["sample"],r["node1"],r["node2"]) for r in rows]
        if len(keys) != len(set(keys)):
            raise ValueError("Duplicate spatial edge keys")
        for metric in ["raw_spearman", "partial_rank_total_melanoma_other_nuisance", "parent_fraction_spearman"]:
            values = [finite(r[metric]) for r in rows]
            good = [v for v in values if v is not None]
            spatial_rows.append({
                "cm": cm, "posterior_abundance_summary": posterior, "assumed_N": n,
                "metric": metric, "edge_section_settings": len(rows),
                "finite_settings": len(good), "missing_settings": len(rows)-len(good),
                "positive_settings": sum(v > 0 for v in good),
                "negative_settings": sum(v < 0 for v in good),
                "zero_settings": sum(v == 0 for v in good),
                "interpretation": "descriptive;shared_sections_edges_and_N_not_independent;no_new_P",
            })
    patient = read_json("differential_LR_three_source_joint_v1/complete.json")
    if not patient["status"].startswith("JOINT_BH_COMPLETE"):
        raise ValueError("Patient LR joint correction not complete")
    receipt = {
        "status": "PASS_EXISTING_AGGREGATE_EVIDENCE_REGISTER_WITH_CURRENT_HASH_CHECKS",
        "UTC": datetime.now(timezone.utc).isoformat(),
        "cm_version": "2026-11-15_R4_K15", "modules": 15,
        "current_linked_upstream_audits": 15, "spatial_fits": 24,
        "hash_checked_artifacts": checked_files + 2,
        "CPDB_events": ca["joint_significant_events"],
        "patient_LR_variable_tests": patient["variable_tests"],
        "patient_LR_q05": patient["q05"],
        "input_sha256": provenance,
        "script_sha256": digest(Path(__file__)),
        "scope": "Existing aggregate results only. No native rows exported; no new fitting, BH family, P value, endpoint selection or causal claim.",
        "supersession": "Current hash-linked upstream audits supersede early pending snapshots; all24 comparison supersedes archive-time18/24.",
        "unverified": ["raw numerical replay beyond upstream audited scope", "universal patient independence", "histology/segmentation validation", "external clinical prediction"],
    }
    output.mkdir(parents=True, exist_ok=False)
    write_table(output / "module_evidence.tsv", module_rows)
    write_table(output / "spatial_edge_setting_coverage.tsv", spatial_rows)
    write_table(output / "spatial_N_sensitivity_existing_summary.tsv", sensitivity)
    receipt["output_sha256"] = {p.name: digest(p) for p in output.glob("*.tsv")}
    (output / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({k: receipt[k] for k in ["status", "UTC", "modules", "spatial_fits", "hash_checked_artifacts", "patient_LR_variable_tests", "patient_LR_q05"]}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build(args.validation_root, args.output)
