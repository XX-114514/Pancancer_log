#!/usr/bin/env python3
"""Read the five fixed CosMx inputs; emit only reviewed aggregate metadata.

Requires h5py, numpy and scipy. The inventory stays outside the log repository.
This is a read-only observation, not a cell annotation or spatial inference run.
"""
import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_matrix


def read_array(node):
    if isinstance(node, h5py.Group):
        categories = node["categories"][()].astype(str)
        codes = node["codes"][()]
        result = np.full(len(codes), "<missing>", dtype=object)
        present = codes >= 0
        result[present] = categories[codes[present]]
        return result
    return node[()]


def scalar(node):
    value = node[()]
    if isinstance(value, bytes):
        return value.decode()
    return value.item() if hasattr(value, "item") else value


def quantiles(values):
    return np.quantile(values, [0, 0.5, 1]).astype(float).tolist()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--converter", type=Path, required=True)
    args = parser.parse_args()
    started = datetime.now(timezone.utc).isoformat()
    inventory = args.inventory.read_bytes()
    rows = list(csv.DictReader(inventory.decode().splitlines(), delimiter="\t"))
    selected = [r for r in rows if r["dataset_id"] in {"GSE287472", "GSE299786"}
                and r["path"].endswith(".cosmx-smi.h5ad")]
    assert len(selected) == 5 and len({r["path"] for r in selected}) == 5
    records = []
    for row in selected:
        path = Path(row["path"])
        before = path.stat()
        with h5py.File(path, "r") as handle:
            x = handle["X"]
            assert x.attrs["encoding-type"] == "csr_matrix"
            shape = tuple(int(v) for v in x.attrs["shape"])
            data, indices, indptr = (x[k][()] for k in ["data", "indices", "indptr"])
            matrix = csr_matrix((data, indices, indptr), shape=shape)
            totals = np.asarray(matrix.sum(axis=1)).ravel()
            obs = handle["obs"]
            coords = handle["obsm"]["spatial"][()]
            px = np.column_stack([read_array(obs[k]) for k in
                                  ["CenterX_global_px", "CenterY_global_px"]])
            qc = {key: quantiles(read_array(obs[key])) for key in
                  ["Area", "Area.um2", "nCount_RNA", "nCount_negprobes",
                   "nCount_falsecode", "unassignedTranscripts"] if key in obs}
            mapping_fields = [k for k in obs if any(token in k.lower() for token in
                              ["patient", "donor", "core", "subject"])]
            result = {
                "dataset_id": row["dataset_id"], "public_file_name": path.name,
                "cancer_label": row["cancer"], "cells": shape[0], "genes": shape[1],
                "input_bytes": before.st_size, "input_mtime_ns": before.st_mtime_ns,
                "matrix_dtype": str(data.dtype), "stored_nonzero_entries": len(data),
                "stored_counts_nonnegative_integer": bool(np.all(data >= 0) and
                                                         np.all(data == np.floor(data))),
                "RNA_row_total_min_median_max": quantiles(totals),
                "RNA_row_total_lt30_cells_descriptive_only": int((totals < 30).sum()),
                "zero_RNA_cells": int((totals == 0).sum()),
                "finite_coordinate_cells": int(np.isfinite(coords).all(axis=1).sum()),
                "global_coordinate_duplicate_rows": int(len(coords) - len(np.unique(coords, axis=0))),
                "stored_spatial_units": scalar(handle["uns"]["spatial_units"]),
                "converter_default_px_to_um": 0.12028,
                "max_coordinate_residual_default_scale": float(np.max(np.abs(coords-px*0.12028))),
                "physical_calibration_independently_verified": False,
                "fov_unique_count": int(len(np.unique(read_array(obs["fov"])))),
                "file_sample_key_unique_count": int(len(np.unique(read_array(obs["sample_id"])))),
                "donor_core_candidate_field_names": mapping_fields,
                "patient_mapping_reviewed": False,
                "control_features_dropped_count": int(scalar(handle["uns"]["control_features_dropped"])),
                "negative_probe_obs_available": "nCount_negprobes" in obs,
                "QC_field_min_median_max": qc,
            }
            if "nCount_RNA" in obs:
                result["stored_total_matches_obs_RNA"] = bool(np.array_equal(totals, read_array(obs["nCount_RNA"])))
            records.append(result)
        after = path.stat()
        assert (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns)
    payload = {
        "schema": "spatial_readonly_metadata_checkpoint_v1",
        "event_started_UTC": started,
        "event_ended_UTC": datetime.now(timezone.utc).isoformat(),
        "inventory_sha256": hashlib.sha256(inventory).hexdigest(),
        "converter_sha256": hashlib.sha256(args.converter.read_bytes()).hexdigest(),
        "collector_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": {"files": 5, "cells": sum(r["cells"] for r in records),
                  "fixed_contact_records_pending_expression_preflight": 13},
        "checks": {"fixed_file_count": "PASS", "input_stat_stability": "PASS",
                   "all_stored_counts_nonnegative_integer": all(r["stored_counts_nonnegative_integer"] for r in records),
                   "coordinates_match_converter_default": all(r["max_coordinate_residual_default_scale"] < 1e-6 for r in records)},
        "files": records,
        "scientific_release": False,
        "spatial_contact_validated": False,
        "limitations": ["File/sample keys and FOVs are not patient identifiers.",
                        "Physical calibration has not been independently verified.",
                        "Count>0 has not been checked against negative-probe background.",
                        "No ligand-receptor endpoint identities or proximity tests were run.",
                        "Read-only metadata checks are not an independent scientific audit."]
    }
    print(json.dumps(payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
