#!/usr/bin/env python3
"""Build metadata/samples.tsv: one row per directory holding a *cell_feature_matrix.h5.

Stdlib only, so it can run before the analysis env exists.
Usage: META_ROOT=... python3 -I make_samples.py RAW_ROOT OUT_TSV  (config.sh exports META_ROOT)
"""
import csv
import os
import sys
from pathlib import Path

# Context per source prefix (from the download manifest naming; patient/core units are not inferred here).
SOURCE_CONTEXT = {
    "C05": ("lung_cancer", "NSCLC FFPE", "10x Prime 5K", "single section"),
    "C11": ("normal_skin", "non-diseased skin FFPE", "10x Xenium v1 skin", "negative-context control"),
    "S07": ("breast_cancer", "FFPE breast cancer (Janesick 2023)", "10x Xenium v1 breast", "technical replicate sections"),
    "S10": ("ovarian_cancer", "FFPE ovarian cancer", "10x Xenium v1 + IO add-on", "single section"),
    "S11": ("colorectal_cancer", "FFPE colorectal cancer", "10x Xenium v1 + IO add-on", "single section"),
    "S12": ("lung_cancer", "FFPE lung cancer", "10x Xenium v1 + IO add-on", "single section"),
    "S16": ("LUAD_or_mesothelioma", "TMA (GSE300007)", "Xenium uni/multi-modal", "TMA: core-to-patient map pending"),
    "S18": ("HNSCC", "fresh-frozen HNSCC (GSE300147)", "10x Xenium v1 multi-tissue 377 + add-on", "patient unit; run 2 = replicate section"),
}


# GSM -> patient/run from fetch_geo_samples.py tables (META_ROOT/geo/samples_GSE*.tsv);
# replicate runs must collapse to one patient unit downstream.
def load_gsm_meta():
    meta_root = os.environ.get("META_ROOT")
    meta = {}
    if not meta_root:
        return meta
    for tsv in sorted(Path(meta_root, "geo").glob("samples_GSE*.tsv")):
        gse = tsv.stem.replace("samples_", "")
        with open(tsv) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                meta[r["gsm"]] = dict(r, gse=gse)
    return meta


def one(files, pattern):
    hits = sorted(f for f in files if f.name.endswith(pattern) or pattern in f.name)
    return str(hits[0]) if hits else ""


def main(raw_root, out_tsv):
    raw_root = Path(raw_root)
    gsm_meta = load_gsm_meta()
    rows = []
    for h5 in sorted(raw_root.rglob("*cell_feature_matrix.h5")):
        d = h5.parent
        files = [f for f in d.iterdir() if f.is_file()]
        rel = d.relative_to(raw_root)
        source = rel.parts[0]
        prefix = source.split("_")[0]
        context, tissue, platform, unit_note = SOURCE_CONTEXT.get(prefix, ("unknown",) * 4)
        sample_id = "__".join(rel.parts)
        if "multi_modal" in h5.name:
            segmentation = "multimodal"
        elif "unimodal" in h5.name:
            segmentation = "unimodal"
        else:
            segmentation = "default"
        if prefix == "S16":
            context = "mesothelioma" if "Mesothelioma" in h5.name else "LUAD"
        rows.append({
            "sample_id": sample_id,
            "source": source,
            "context": context,
            "tissue": tissue,
            "platform": platform,
            "segmentation": segmentation,
            "unit_note": unit_note,
            "h5": str(h5),
            "cells": one(files, "cells.parquet.gz") or one(files, "cells.parquet"),
            "analysis_tar": one(files, "analysis.tar.gz"),
            "author_supervised": one([f for f in files if "unsupervised" not in f.name], "_supervised.csv.gz"),
            "author_unsupervised": one(files, "unsupervised.csv.gz"),
            "gene_panel": one(files, "gene_panel.json"),
            "patient": "",
            "run": "",
        })
        meta = gsm_meta.get(rel.parts[-1])
        if meta and meta["patient"]:
            rows[-1]["patient"] = meta["gse"] + "_" + meta["patient"].replace("Patient ", "P").replace(" ", "")
            rows[-1]["run"] = meta["run"]
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} samples -> {out_tsv}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
