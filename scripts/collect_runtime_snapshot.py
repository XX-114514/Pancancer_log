#!/usr/bin/env python3
"""Read-only, aggregate-only PanCancer runtime exporter (Python >=3.10).

Run on the data host; output contains no environment, credentials, cell rows,
patient rows, full process commands, or complete runtime logs.
"""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

EXPANSION = "Project_v5/runs/20261010_clinical_scrna_expansion_v1"
ANNOTATION = "Project_v3/annotation_rebuild"
NAMES = ["CAR-T治疗产品", "卵巢癌SPECTRUM", "大B细胞淋巴瘤", "儿童高级别胶质瘤",
         "转移性乳腺癌", "未治疗TNBC", "CTCL图谱"]
GEO = [("GSE286827", "头颈癌免疫治疗"), ("GSE294109", "肾癌治疗"), ("GSE292046", "胶质瘤疫苗")]
SCALARS = {"status", "schema", "run_id", "version", "samples", "nodes", "CM", "K", "pairs", "cells", "GSE",
           "interactions", "primary_prioritized_CM_rows", "summary_rows", "rows", "full_rows", "included_rows", "n_cells", "candidate_cells", "cells_all", "nodes_all", "candidate_nodes",
           "samples_min100", "samples_before_min100", "samples_before_min100_filter", "cells_after_source_hold",
           "cells_before_downstream_holds", "nodes_before_min100_filter", "all_input_LIANA_rows", "unique_within_CM_rows",
           "prioritized_primary_CM_rows", "prioritized_summary_rows", "all_LIANA_rows", "within_CM_rows", "primary_sample_count",
           "global_nodes", "global_edges", "CM_node_rows", "CM_edge_rows", "new_edges_vs_reference", "removed_edges_vs_reference",
           "release_eligible", "ecology_release_eligible", "scientific_release", "identity_promoted", "identity_or_causal_release",
           "biological_identity_changed", "NMF_recomputed", "NMF_W_H_preserved", "rank_survey_rows", "K_optimality_review_pending",
           "independent_patient_inference_claimed", "reviewed_samples", "reviewed_patient_units", "unique_composition_pairs",
           "CM_edge_records", "n_CM", "candidates_inferential", "disputed_descriptive_candidates", "composition_estimated_strata",
           "primary_model_estimated_strata", "primary_model_BH_significant_strata", "canonical_annotation_rows",
           "canonical_write_status", "visual_review_status", "artifact_task_status", "annotation_run_id", "downstream_run_id",
           "execution_status", "scientific_identity_status", "scientific_identity_decision", "exit_code", "completed_utc",
           "scientific_identity_promotion", "causal_release", "sorting_unknown", "patient_identity_is_not_complete",
           "native_gene_policy", "CM_membership_reused", "no_old_CM_membership", "rank_cutoffs_are_not_FDR"}

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def read(path):
    if not path.is_file():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def scalar_summary(obj):
    return {k: v for k, v in obj.items() if k in SCALARS and isinstance(v, (str, int, float, bool, type(None)))}

def logical(path, root):
    try:
        return "${PANCANCER_ROOT}/" + str(Path(path).relative_to(root))
    except ValueError:
        return "external:" + Path(path).name

def evidence(path, root):
    if not path.is_file():
        return {"path": logical(path, root), "exists": False}
    return {"path": logical(path, root), "exists": True, "bytes": path.stat().st_size,
            "sha256": sha(path), "modified_utc": dt.datetime.fromtimestamp(path.stat().st_mtime, dt.timezone.utc).isoformat()}

def tail(path):
    if not path.is_file():
        return ""
    with path.open("rb") as f:
        f.seek(max(0, path.stat().st_size - 24000))
        return f.read().decode("utf-8", "replace")

def progress(path):
    text = re.sub(r"\x1b\[[0-9;]*m", "", tail(path))
    matches = re.findall(r"\[#[a-z0-9]+\s+([^\]\n]+)\]", text)
    if not matches:
        return {"state": "no_progress_observed"}
    m = re.search(r"(\S+/\S+)\((\d+)%\).*?DL:(\S+)", matches[-1])
    if not m:
        return {"state": "progress_format_unrecognized"}
    age = max(0, dt.datetime.now().timestamp() - path.stat().st_mtime)
    return {"transfer": m[1], "percent": int(m[2]), "rate_per_second": m[3],
            "log_age_seconds": round(age), "stale_over_10min": age > 600,
            "not_a_download_completion_proof": True}

def proof_summary(proof):
    keys = {"status", "sha256", "source_sha256", "expected_bytes", "actual_bytes", "n_raw", "n_qc_nonduplicate", "n_duplicate",
            "source_cells", "n_source_cells", "n_genes", "n_original_library_units", "n_libraries", "counts_sum", "counts_sum_before",
            "counts_sum_after", "duplicate_symbols_collapsed", "floor_pass_flagged", "high_mt_flagged", "baseline_exact_matches",
            "within_source_exact_matches", "n_full_source_cells_audited", "biological_filter_applied", "cells_removed",
            "final_accepted_cells", "final_accepted_new_cells", "dedup_status", "qc_policy", "official_source_metadata_sha256"}
    result = {k: v for k, v in proof.items() if k in keys and isinstance(v, (str, int, float, bool, type(None)))}
    if "tier_counts" in proof:
        result["tier_counts"] = {k: v for k, v in proof["tier_counts"].items() if k in {"A", "B", "C"}}
    return result

def processes(root):
    rows = []
    ps = subprocess.run(["ps", "-eo", "pid,uid,comm,etime,stat,pcpu,args"], capture_output=True, text=True, check=True)
    for line in ps.stdout.splitlines()[1:]:
        fields = line.strip().split(None, 6)
        if len(fields) != 7 or not fields[1].isdigit() or int(fields[1]) != os.getuid():
            continue
        pid, uid, command, elapsed, state, cpu, args = fields
        if command not in {"python", "python3", "R", "Rscript", "aria2c"} or int(pid) == os.getpid():
            continue
        try:
            cwd = Path("/proc") / pid / "cwd"
            cwd = cwd.resolve()
        except OSError:
            cwd = Path(".")
        scope = None
        if str(root / "Project_v3") in args or cwd.is_relative_to(root / "Project_v3"):
            scope = "Project_v3"
        elif str(root / EXPANSION) in args or cwd.is_relative_to(root / EXPANSION):
            scope = "clinical_expansion"
        if scope:
            scripts = re.findall(r"(?:[^\s]+/)?([^/\s]+\.(?:py|R))(?=\s|$)", args)
            script = scripts[-1] if scripts else command
            role = "infrastructure" if any(word in script for word in ["proxy", "bridge", "relay"]) else "monitor" if any(word in script for word in ["readable", "refresh_current", "monitor"]) else "pipeline"
            rows.append({"pid": int(pid), "scope": scope, "role": role, "executable": command, "script_basename": script,
                         "elapsed": elapsed, "state": state, "cpu_percent": float(cpu)})
    return rows

def collect(root):
    run = root / EXPANSION
    ar = root / ANNOTATION
    snapshot = {"schema": "pancancer_cloud_snapshot_v1", "collected_utc": now(), "source_mode": "read_only_live_server",
                "repository_visibility": "public", "observed_processes": processes(root), "evidence": []}
    downloads = []; source_total = 0
    queue = read(run / "metadata/download_queue_v1.json")
    for i, d in enumerate(queue):
        uid = "CXG_" + d["dataset_id"].replace("-", "_")
        p = run / "manifests" / ("download_" + d["dataset_version_id"] + ".json")
        proof = read(p); verified = proof.get("status") == "downloaded_header_verified"
        if verified:
            source_total += d["cell_count"]
        record = {"dataset": uid, "name": NAMES[i], "source_cells": d["cell_count"], "download_verified": verified,
                  "collection_url": "https://cellxgene.cziscience.com/collections/" + d["collection_id"],
                  "download_proof": proof_summary(proof)}
        record["QC"] = proof_summary(read(run / "manifests" / ("accepted_" + uid + ".json")))
        record["native_preprocessing"] = proof_summary(read(run / "artifacts/public_native_preprocessing_20261011" / (uid + "_native_preprocessing.json")))
        record["baseline_overlap"] = proof_summary(read(run / "artifacts/native_baseline_overlap_20261011" / (uid + "_baseline_overlap.json")))
        downloads.append(record); snapshot["evidence"].append(evidence(p, root))
    for g, name in GEO:
        p = run / "manifests" / ("download_public_" + g + ".json")
        proof = read(p); inv = read(run / "artifacts/preintegration_20261011" / (g + "_inventory_summary.json"))
        cells = inv.get("n_source_barcodes", 0)
        if proof.get("files"):
            source_total += cells
        downloads.append({"dataset": g, "name": name, "source_cells": cells, "download_verified": bool(proof.get("files")),
                          "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + g,
                          "QC": proof_summary(read(run / "manifests" / ("accepted_" + g + ".json")))})
        snapshot["evidence"].append(evidence(p, root))
    rmc = read(run / "artifacts/source_audit_20261011/GSE256326_native_counts.json")
    source_total += rmc.get("n_source_cells", 0)
    if rmc:
        downloads.append({"dataset": "GSE256326", "name": "肾髓质癌", "source_cells": rmc["n_source_cells"],
                          "download_verified": True, "native_preprocessing": proof_summary(rmc)})
    qc = [proof_summary(read(p)) for p in (run / "manifests").glob("accepted_*.json")]
    snapshot["expansion"] = {"run_id": run.name, "source_cells_verified_download_or_inventory": source_total,
                             "qc_shards": len(qc), "qc_cells_pending_final_source_acceptance": sum(p.get("n_qc_nonduplicate", 0) for p in qc),
                             "final_accepted_new_cells": 0, "final_acceptance_status": "no_final_acceptance_manifest_registered",
                             "target_after_dedup_QC": [4500000, 5500000], "cross_study_integration_authorized": False,
                             "datasets": downloads, "batch2": []}
    for d in read(run / "metadata/public_clinical_batch2_download_queue.json").get("studies", []):
        files = []
        for url in d["selected_files"]:
            name = url.rsplit("/", 1)[-1]
            proof = read(run / "manifests" / ("batch2_download_" + name + ".json"))
            verified = proof.get("status") == "public_file_downloaded_size_sha256_verified_pending_modality_and_raw_counts"
            files.append({"name": name, "download_verified": verified, "proof": proof_summary(proof),
                          "transfer": {} if verified else progress(run / "logs" / ("batch2_" + name + ".download.log"))})
        snapshot["expansion"]["batch2"].append({"dataset": d["gse_id"], "files": files})
    snapshot["expansion"]["zenodo"] = {"record": "19586163", "source_url": "https://zenodo.org/records/19586163",
         "download_proof": proof_summary(read(run / "manifests/download_zenodo_19586163.json")),
         "transfer": progress(run / "logs/zenodo_19586163_download_v2.log")}
    state_path = ar / "CURRENT_STATE.json"; state = read(state_path)
    snapshot["project_v3"] = {"source_state_generated_utc": state.get("generated_utc"),
        "identity_release_readiness": {k: v for k, v in state.get("identity_release_readiness", {}).items() if k in {
            "status", "full_cells", "current_nonempty_fine_cells", "remaining_GSE299623_groups", "remaining_GSE299623_cells", "approved_remaining_group_decisions", "release_eligible"}},
        "accepted_identity_assembly": {k: v for k, v in state.get("accepted_baseline", {}).items() if k in {"rows", "fine_nonempty", "sha256"}},
        "expression_annotation": scalar_summary(state.get("uniform_expression_annotation_and_downstream", {})),
        "run_evidence": []}
    snapshot["evidence"].append(evidence(state_path, root))
    selections = {
        "20261010_uniform_expression_annotation_v1": ["input_contract.json", "run_register.json", "manifest.json", "independent_audit.json", "final_goal_scope_audit.json", "final_goal_handoff_registration.json"],
        "20261010_expression_annotation_downstream_v1": ["run_contract.json", "run_register.json", "input_independent_audit.json",
           "covarnet_audit_v1/manifest.json", "liana/independent_audit.json",
           "cm_liana_summary_v1/independent_audit.json", "covarnet_network_corrected_v1/network_independent_audit.json",
           "downstream_review_v2/table_and_render_audit.json", "downstream_review_v2/agent_visual_review.json"],
        "20261010_CM_source_replication_v1": ["run_register.json", "results_summary.json", "execution_final.json", "audits/independent_audit.json", "audits/meta_readback.json", "audits/visual_review.json"]}
    for rid, files in selections.items():
        for rel in files:
            p = ar / rid / rel
            meta = evidence(p, root); data = read(p)
            snapshot["project_v3"]["run_evidence"].append({"run_id": rid, "artifact": rel, "file": meta, "summary": scalar_summary(data)})
    clinical = {}
    allowed = {"source_cells_metadata", "n_source_cells_metadata", "n_libraries", "n_patients", "n_exact_unique_GSM_matches",
               "source_patients", "exact_patient_matches", "OS_PFS_duration_event_available", "candidate_barcode_prefixes",
               "orig_ident_groups", "n_source_cells", "n_source_accessions", "n_nonprimary_quarantined_pending_review",
               "n_normal_labels_pending_context", "exact_accession_overlap_count", "clinical_tier", "response_column_present"}
    for name in ["GSE236581_clinical_audit.json", "GSE256326_clinical_audit_v2.json", "SPECTRUM_clinical_audit.json", "TNBC_provenance_hold.json", "CTCL_provenance_audit.json"]:
        p = run / "artifacts/source_audit_20261011" / name
        clinical[name] = {k: v for k, v in read(p).items() if k in allowed}
        snapshot["evidence"].append(evidence(p, root))
    snapshot["expansion"]["clinical_source_audits"] = clinical
    return snapshot

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output already exists: preserve snapshots and choose a new dated filename")
    root = args.project_root.resolve()
    if not (root / EXPANSION).is_dir() or not (root / ANNOTATION).is_dir():
        parser.error("Missing expected expansion or Project_v3 annotation directories")
    result = collect(root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output_basename": args.output.name, "collected_utc": result["collected_utc"]}))

if __name__ == "__main__":
    main()
