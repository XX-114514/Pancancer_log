#!/usr/bin/env python3
"""Convert one Xenium sample (author/10x processed outputs) into a counts AnnData plus a QC summary.

Raw counts are kept unmodified in X (gene-expression features only); control features are summed
into obs. Cells are flagged, never removed.
Usage: python -I scripts/ingest_xenium.py SAMPLES_TSV ROW_INDEX(1-based) OUT_ROOT
"""
import csv
import gzip
import io
import json
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scanpy as sc

MIN_TRANSCRIPTS = 10  # 10x default QC threshold; used as a flag only
CONTROL_TYPES = {
    "Negative Control Probe": "neg_probe_counts",
    "Negative Control Codeword": "neg_codeword_counts",
    "Unassigned Codeword": "unassigned_codeword_counts",
    "Blank Codeword": "unassigned_codeword_counts",
    "Deprecated Codeword": "deprecated_codeword_counts",
    "Genomic Control": "genomic_control_counts",
}


def read_cells(path):
    src = io.BytesIO(gzip.open(path).read()) if path.endswith(".gz") else path
    cells = pq.read_table(src).to_pandas()
    cells["cell_id"] = cells["cell_id"].astype(str)
    return cells.set_index("cell_id")


def read_graphclust(tar_path):
    with tarfile.open(tar_path) as tf:
        member = next(n for n in tf.getnames() if n.endswith("gene_expression_graphclust/clusters.csv"))
        df = pd.read_csv(tf.extractfile(member))
    return df.set_index(df.columns[0].strip())[df.columns[1]].astype(str).rename(index=str)


def read_author(path):
    # Read as str: a numeric barcode column with blanks would otherwise parse as float ("1.0").
    df = pd.read_csv(path, dtype=str)
    df = df[df.iloc[:, 0].notna()]
    return df.set_index(df.columns[0])[df.columns[1]]


def attach(adata, labels, key, qc):
    common = adata.obs_names.intersection(labels.index)
    adata.obs[key] = labels.reindex(adata.obs_names).astype("category")
    qc[f"{key}_matched_cells"] = int(len(common))
    qc[f"{key}_label_rows"] = int(len(labels))


def main(samples_tsv, row_index, out_root):
    row = list(csv.DictReader(open(samples_tsv), delimiter="\t"))[int(row_index) - 1]
    sid = row["sample_id"]
    out_root = Path(out_root)
    (out_root / "h5ad").mkdir(parents=True, exist_ok=True)
    (out_root / "qc").mkdir(parents=True, exist_ok=True)

    full = sc.read_10x_h5(row["h5"], gex_only=False)
    full.var_names_make_unique()
    ftype = full.var["feature_types"].astype(str)
    adata = full[:, (ftype == "Gene Expression").values].copy()
    for t, col in CONTROL_TYPES.items():
        mask = (ftype == t).values
        if mask.any():
            vals = np.asarray(full[:, mask].X.sum(axis=1)).ravel()
            adata.obs[col] = adata.obs.get(col, 0) + vals
    control_features = {t: int((ftype == t).sum()) for t in CONTROL_TYPES if (ftype == t).any()}
    del full

    cells = read_cells(row["cells"])
    keep = [c for c in ["x_centroid", "y_centroid", "cell_area", "nucleus_area", "nucleus_count",
                        "segmentation_method", "transcript_counts"] if c in cells.columns]
    n_cells_table = len(cells)
    cells = cells.reindex(adata.obs_names)[keep]
    adata.obs = adata.obs.join(cells.rename(columns={"transcript_counts": "tenx_transcript_counts"}))
    missing_xy = int(adata.obs["x_centroid"].isna().sum())
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].to_numpy()

    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None, log1p=False)
    adata.obs["qc_pass"] = adata.obs["total_counts"] >= MIN_TRANSCRIPTS
    # patient/run exist only where the source publishes them (e.g. S18); empty otherwise.
    meta_cols = ["sample_id", "source", "context", "platform", "segmentation", "patient", "run"]
    for c in meta_cols:
        adata.obs[c] = pd.Categorical([row.get(c, "")] * adata.n_obs)
    adata.uns["sample"] = {k: row.get(k, "") for k in ["sample_id", "source", "context", "tissue", "platform",
                                                       "segmentation", "unit_note", "patient", "run"]}

    qc = {"sample_id": sid, "context": row["context"], "patient": row.get("patient", ""), "platform": row["platform"],
          "segmentation": row["segmentation"], "n_cells": int(adata.n_obs), "n_genes_panel": int(adata.n_vars),
          "control_features": control_features, "cells_table_rows": n_cells_table,
          "cells_missing_centroid": missing_xy, "min_transcripts_flag": MIN_TRANSCRIPTS}
    if row["analysis_tar"]:
        attach(adata, read_graphclust(row["analysis_tar"]), "tenx_graphclust", qc)
    for key in ["author_supervised", "author_unsupervised"]:
        if row[key]:
            attach(adata, read_author(row[key]), key, qc)

    obs = adata.obs
    gene_total = float(obs["total_counts"].sum())
    neg_cw = float(obs["neg_codeword_counts"].sum()) if "neg_codeword_counts" in obs else float("nan")
    n_neg_cw = control_features.get("Negative Control Codeword", 0)
    qc.update({
        "median_counts": float(obs["total_counts"].median()),
        "median_genes": float(obs["n_genes_by_counts"].median()),
        "median_cell_area": float(obs["cell_area"].median()) if "cell_area" in obs else None,
        "frac_cells_below_flag": float((~obs["qc_pass"]).mean()),
        "zero_count_cells": int((obs["total_counts"] == 0).sum()),
        # Per-feature background vs signal: mean negative-codeword counts per codeword relative to mean counts per gene.
        "neg_codeword_per_feature_ratio": (neg_cw / n_neg_cw) / (gene_total / adata.n_vars) if n_neg_cw and gene_total else None,
        "x_range_um": [float(obs["x_centroid"].min()), float(obs["x_centroid"].max())],
        "y_range_um": [float(obs["y_centroid"].min()), float(obs["y_centroid"].max())],
    })

    for col in obs.columns:
        if obs[col].dtype == object:
            obs[col] = obs[col].astype(str)
    out_h5ad = out_root / "h5ad" / f"{sid}.counts.h5ad"
    adata.write_h5ad(out_h5ad, compression="gzip")
    qc["output"] = str(out_h5ad.relative_to(out_root))
    (out_root / "qc" / f"{sid}.qc.json").write_text(json.dumps(qc, indent=2))
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:4])
