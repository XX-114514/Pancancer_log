#!/usr/bin/env python3
"""Panel qualification across all ingested samples: which frozen gene panels each Xenium panel covers,
per-gene detection, and a cross-sample QC summary. Descriptive only; no states or effects are fitted.
Usage: python -I scripts/panel_coverage.py OUT_ROOT FROZEN_PLAN_JSON
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

# Lists not in the frozen plan are labelled by provenance and are not themselves frozen definitions.
EXTRA_LISTS = {
    "CM02_core_display_markers": (["CTLA4", "GZMB", "GZMK", "LTB", "IL7R", "CTSW", "LEF1"],
                                  "README: 2026-10-14 immune atlas CM02 v1 core member display markers"),
    "LN_state_RNA_analogs": (["CD4", "CD8A", "FOXP3", "CTLA4", "PDCD1", "HAVCR2", "MKI67", "CCR7", "IL7R"],
                             "RNA analogs of NSCLC-LN protein states; protein-to-RNA bridge not yet frozen"),
}
CD3 = ["CD3D", "CD3E", "CD3G"]


def gene_lists(frozen_plan):
    plan = json.loads(Path(frozen_plan).read_text())
    lists = {f"frozen::{k}": (v, "frozen_plan_before_effects.json /panels") for k, v in plan["panels"].items()}
    lists.update(EXTRA_LISTS)
    return lists


def main(out_root, frozen_plan):
    out_root = Path(out_root)
    lists = gene_lists(frozen_plan)
    wanted = sorted({g for genes, _ in lists.values() for g in genes} | set(CD3))
    gene_rows, set_rows, abund_rows = [], [], []
    qc = pd.DataFrame([json.loads(p.read_text()) for p in sorted((out_root / "qc").glob("*.qc.json"))])

    for sid in qc["sample_id"]:
        a = ad.read_h5ad(out_root / "h5ad" / f"{sid}.counts.h5ad")
        a = a[a.obs["qc_pass"].values]
        present = [g for g in wanted if g in a.var_names]
        sub = a[:, present].X
        sub = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
        det = dict(zip(present, (sub > 0).mean(axis=0)))
        mean = dict(zip(present, sub.mean(axis=0)))
        for g in wanted:
            gene_rows.append({"sample_id": sid, "gene": g, "in_panel": g in det,
                              "frac_cells_detected": det.get(g, np.nan), "mean_counts": mean.get(g, np.nan)})
        for name, (genes, provenance) in lists.items():
            hit = [g for g in genes if g in det]
            set_rows.append({"sample_id": sid, "gene_list": name, "n_present": len(hit), "n_total": len(genes),
                             "missing": ",".join(g for g in genes if g not in det), "provenance": provenance})
        cd3_cols = [present.index(g) for g in CD3 if g in present]
        cd3_pos = (sub[:, cd3_cols] > 0).any(axis=1) if cd3_cols else np.zeros(a.n_obs, bool)
        row = {"sample_id": sid, "qc_pass_cells": int(a.n_obs), "CD3_genes_in_panel": len(cd3_cols),
               "frac_CD3_any_pos": float(cd3_pos.mean()) if cd3_cols else np.nan}
        for g in ["CD8A", "CD4", "FOXP3", "CTLA4", "GZMB", "PDCD1", "HAVCR2", "MKI67"]:
            row[f"frac_{g}_pos"] = float(det[g]) if g in det else np.nan
        abund_rows.append(row)

    summary_dir = out_root / "summary"
    summary_dir.mkdir(exist_ok=True)
    genes = pd.DataFrame(gene_rows)
    sets = pd.DataFrame(set_rows)
    abund = pd.DataFrame(abund_rows)
    qc_cols = ["sample_id", "context", "platform", "segmentation", "n_cells", "n_genes_panel", "median_counts",
               "median_genes", "median_cell_area", "frac_cells_below_flag", "zero_count_cells",
               "neg_codeword_per_feature_ratio", "cells_missing_centroid"]
    qc_out = qc[[c for c in qc_cols if c in qc]]
    genes.to_csv(summary_dir / "panel_gene_detection.tsv", sep="\t", index=False)
    sets.to_csv(summary_dir / "panel_set_coverage.tsv", sep="\t", index=False)
    abund.to_csv(summary_dir / "T_marker_abundance.tsv", sep="\t", index=False)
    qc_out.to_csv(summary_dir / "qc_summary.tsv", sep="\t", index=False)

    cov = sets.assign(cov=sets["n_present"].astype(str) + "/" + sets["n_total"].astype(str)) \
              .pivot(index="sample_id", columns="gene_list", values="cov")
    report = ["# Xenium panel qualification and QC", "",
              f"Frozen plan: `{Path(frozen_plan).name}`; cells flagged qc_pass = total gene counts >= "
              f"{int(qc['min_transcripts_flag'].iloc[0])} (flag only, nothing removed).", "",
              "## QC summary", "", qc_out.round(4).to_markdown(index=False), "",
              "## Gene-list coverage (present / total)", "", cov.to_markdown(), "",
              "## T-cell marker detection among qc_pass cells (descriptive, not a cell-type call)", "",
              abund.round(4).to_markdown(index=False), "",
              "## Missing genes per list", "",
              sets[sets["missing"] != ""][["sample_id", "gene_list", "missing"]].to_markdown(index=False), ""]
    (summary_dir / "REPORT.md").write_text("\n".join(report))
    print("\n".join(report))


if __name__ == "__main__":
    main(*sys.argv[1:3])
