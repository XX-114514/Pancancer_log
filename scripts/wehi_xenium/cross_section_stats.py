#!/usr/bin/env python3
"""Cross-section statistics after section_stats.py: measurement mask, unit structure and replicate concordance.

- Measurement mask: per section x gene list, each gene is not_in_panel / at_background / above_background
  (section_stats.py rule). A list is measurable when every gene is above background, partial when at least one
  is, unmeasurable otherwise. Source level takes the worst section. Missing genes are never dropped.
- Unit structure: S07 technical replicates, S18 same-patient runs, S16 unimodal vs multimodal segmentation of the
  same TMA. Pseudobulk log1p CPM over genes shared by both sections; Pearson and Spearman; same for per-gene
  detection fraction. S18 also compares within-patient pairs with all between-patient pairs.
- S18 patient view: sections, qc_pass cells and pooled medians per patient.
Descriptive only: no states, no effects.
Usage: python -I cross_section_stats.py SAMPLES_TSV OUT_ROOT FROZEN_PLAN_JSON
"""
import importlib.util
import itertools
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ORDER = ["measurable", "partial", "unmeasurable"]


def load_gene_lists(frozen_plan):
    spec = importlib.util.spec_from_file_location("panel_coverage", Path(__file__).with_name("panel_coverage.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.gene_lists(frozen_plan)


def log_cpm(g):
    return np.log1p(g["total_counts"] / g["total_counts"].sum() * 1e6)


def pair_stats(g1, g2):
    shared = g1.index.intersection(g2.index)
    x, y = log_cpm(g1).loc[shared], log_cpm(g2).loc[shared]
    d1, d2 = g1.loc[shared, "frac_cells_detected"], g2.loc[shared, "frac_cells_detected"]
    bg = g1.loc[shared, "above_background"] == g2.loc[shared, "above_background"]
    return {"n_shared_genes": int(len(shared)), "pearson_logcpm": pearsonr(x, y)[0], "spearman_logcpm": spearmanr(x, y)[0],
            "pearson_detection": pearsonr(d1, d2)[0], "background_call_agreement": float(bg.mean())}


def describe(v):
    v = pd.Series(v, dtype=float)
    return f"n={len(v)}, Pearson median {v.median():.4f} (range {v.min():.4f}–{v.max():.4f})" if len(v) else "n=0"


def tma_key(h5):
    return re.sub(r"_Xenium_(unimodal|multi_modal)_cell_feature_matrix\.h5$", "", Path(h5).name.split("_", 1)[1])


def main(samples_tsv, out_root, frozen_plan):
    out_root = Path(out_root)
    samples = pd.read_csv(samples_tsv, sep="\t", dtype=str, keep_default_na=False)
    stats_dir, summary = out_root / "stats", out_root / "summary"
    summary.mkdir(exist_ok=True)
    genes = {s: pd.read_csv(stats_dir / s / "gene_stats.tsv", sep="\t", index_col="gene") for s in samples["sample_id"]}
    sec = pd.DataFrame([json.loads((stats_dir / s / "section_stats.json").read_text()) for s in samples["sample_id"]])
    for key in ["tenx_graphclust", "author_supervised", "author_unsupervised"]:
        col = f"vs_{key}"
        if col in sec:
            sec[f"ARI_{key}"] = sec[col].map(lambda d: d["ARI"] if isinstance(d, dict) else np.nan)
            sec[f"NMI_{key}"] = sec[col].map(lambda d: d["NMI"] if isinstance(d, dict) else np.nan)
            sec = sec.drop(columns=col)
    sec = samples[["sample_id", "source", "context", "segmentation", "patient", "run"]].merge(sec, on="sample_id")
    sec.to_csv(summary / "section_stats.tsv", sep="\t", index=False)

    # Measurement mask
    lists = load_gene_lists(frozen_plan)
    gene_rows, mask_rows = [], []
    for sid, g in genes.items():
        for name, (members, _) in lists.items():
            status = []
            for gene in members:
                s = "not_in_panel" if gene not in g.index else (
                    "above_background" if g.at[gene, "above_background"] else "at_background")
                status.append(s)
                gene_rows.append({"sample_id": sid, "gene_list": name, "gene": gene, "status": s,
                                  "ratio_to_mean_negative": g.at[gene, "ratio_to_mean_negative"] if gene in g.index else np.nan,
                                  "frac_cells_detected": g.at[gene, "frac_cells_detected"] if gene in g.index else np.nan})
            n_ok = status.count("above_background")
            mask_rows.append({"sample_id": sid, "gene_list": name, "n_genes": len(members), "n_above_background": n_ok,
                              "n_at_background": status.count("at_background"), "n_not_in_panel": status.count("not_in_panel"),
                              "status": "measurable" if n_ok == len(members) else ("partial" if n_ok else "unmeasurable"),
                              "not_above_background": ",".join(m for m, s in zip(members, status) if s != "above_background")})
    gene_mask = pd.DataFrame(gene_rows)
    mask = pd.DataFrame(mask_rows).merge(samples[["sample_id", "source"]], on="sample_id")
    gene_mask.to_csv(summary / "measurement_mask_genes.tsv", sep="\t", index=False)
    mask.to_csv(summary / "measurement_mask_sections.tsv", sep="\t", index=False)
    worst = mask.assign(rank=mask["status"].map(ORDER.index)).groupby(["source", "gene_list"])["rank"].max()
    src_mask = worst.map(lambda r: ORDER[r]).unstack("gene_list")
    src_mask.to_csv(summary / "measurement_mask_sources.tsv", sep="\t")

    # Unit-structure pairs
    pairs = []
    s07 = samples[samples["source"] == "S07_GSE243280_breast"]["sample_id"].tolist()
    if len(s07) == 2:
        pairs.append(("S07 technical replicate", "S07", *s07))
    s16 = samples[samples["source"] == "S16_GSE300007_LUAD_meso"].assign(tma=lambda d: d["h5"].map(tma_key))
    for tma, d in s16.groupby("tma"):
        uni, multi = d[d["segmentation"] == "unimodal"]["sample_id"], d[d["segmentation"] == "multimodal"]["sample_id"]
        if len(uni) == 1 and len(multi) == 1:
            pairs.append(("S16 unimodal vs multimodal", tma, uni.iloc[0], multi.iloc[0]))
    s18 = samples[samples["source"] == "S18_GSE300147_HNSCC"]
    for pat, d in s18.groupby("patient"):
        if len(d) == 2:
            pairs.append(("S18 same patient run1 vs run2", pat, *d.sort_values("run")["sample_id"]))
    pair_df = pd.DataFrame([{"comparison": c, "unit": u, "section_a": a, "section_b": b, **pair_stats(genes[a], genes[b])}
                            for c, u, a, b in pairs],
                           columns=["comparison", "unit", "section_a", "section_b", "n_shared_genes", "pearson_logcpm",
                                    "spearman_logcpm", "pearson_detection", "background_call_agreement"])
    between = [pair_stats(genes[a], genes[b])["pearson_logcpm"]
               for a, b in itertools.combinations(s18["sample_id"], 2)
               if s18.set_index("sample_id").at[a, "patient"] != s18.set_index("sample_id").at[b, "patient"]]
    pair_df.to_csv(summary / "unit_concordance_pairs.tsv", sep="\t", index=False)

    # S18 patient view
    s18_sec = sec[sec["source"] == "S18_GSE300147_HNSCC"]
    qc = {s: json.loads((out_root / "qc" / f"{s}.qc.json").read_text()) for s in s18_sec["sample_id"]}
    pat = s18_sec.groupby("patient").agg(n_sections=("sample_id", "size"), qc_pass_cells=("n_cells_qc_pass", "sum"),
                                         all_cells=("n_cells_all", "sum")).reset_index()
    pat["median_counts_range"] = pat["patient"].map(lambda p: "–".join(
        f"{v:g}" for v in sorted({qc[s]["median_counts"] for s in s18_sec[s18_sec["patient"] == p]["sample_id"]})))
    pat.to_csv(summary / "S18_patient_view.tsv", sep="\t", index=False)

    w = pair_df[pair_df["comparison"].str.startswith("S18")]["pearson_logcpm"]
    gl_order = list(lists)
    report = [
        "# Xenium section statistics, measurement mask and unit structure", "",
        "qc_pass cells only (total gene counts >= 10). Descriptive; no cell states or effects. "
        "Leiden clusters are numeric, not identities.", "",
        "## Per-section processing", "",
        "n_features includes S18 TCR clonotype probes; n_genes_above_background counts ordinary genes only.", "",
        sec[["sample_id", "n_cells_qc_pass", "n_features", "n_tcr_clonotype_probes", "n_genes_above_background",
             "n_tcr_probes_above_background", "pca_input", "n_leiden_clusters",
             "largest_cluster_frac", "moran_I_median", "n_genes_moran_bh05_I_gt_0.1"]
            + [c for c in sec if c.startswith("ARI_")]].round(3).to_markdown(index=False), "",
        "S07 Rep1 author_unsupervised (GSM7780153 Fig 3): barcodes run 0–170,429 against 167,780 cells in the "
        "cells table and agree with the supervised labels at chance (ARI ~0.01, also with a +1 shift), so they are "
        "not Xenium cell_ids; its ARI is not interpretable and the labels are unusable until remapped. The supervised "
        "tables write barcode 100000 as '1e+05', so that one cell is unlabelled per file.", "",
        "## Measurement mask by source (worst section)", "",
        "measurable = every gene above the section's largest negative-control feature total; partial = some; "
        "unmeasurable = none.", "",
        src_mask[[c for c in gl_order if c in src_mask]].to_markdown(), "",
        "## Unit-structure concordance (pseudobulk log1p CPM over shared genes)", "",
        pair_df.round(4).to_markdown(index=False), "",
        f"S18 within-patient run pairs: {describe(w)}; between-patient section pairs: {describe(between)}.", "",
        "## S18 patient view", "", pat.to_markdown(index=False), "",
        "## Lists not fully above background, per section", "",
        mask[mask["status"] != "measurable"][["sample_id", "gene_list", "status", "not_above_background"]]
            .to_markdown(index=False), ""]
    (summary / "STATS_REPORT.md").write_text("\n".join(report))
    print("\n".join(report))


if __name__ == "__main__":
    main(*sys.argv[1:4])
