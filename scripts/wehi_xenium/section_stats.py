#!/usr/bin/env python3
"""Per-section descriptive statistics for one ingested Xenium section (qc_pass cells only).

1. Background: per negative-control feature totals re-read from the source h5 (ingest kept only per-cell sums).
   A gene is above_background when its section total exceeds the largest negative-control feature total.
2. Normalisation: library-size to the section median, log1p; raw counts kept in layers['counts'].
3. Numeric clusters: PCA -> kNN graph -> Leiden -> UMAP; Wilcoxon markers; ARI/NMI vs 10x graphclust and
   author labels. Clusters are numbers, not cell identities.
4. Spatial autocorrelation: Moran's I on a 6-NN graph of cell centroids (analytic P, BH within section).
Descriptive only: no states, no effects. Rules are fixed here before results are seen.
Usage: python -I section_stats.py SAMPLES_TSV ROW_INDEX(1-based) OUT_ROOT FROZEN_PLAN_JSON
"""
import csv
import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from pynndescent import PyNNDescentTransformer
from scipy.sparse import csr_matrix
from scipy.stats import false_discovery_control
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

NEG_TYPES = ["Negative Control Probe", "Negative Control Codeword"]
HVG_IF_MORE_THAN = 1000  # targeted panels use every gene; whole-transcriptome-like panels (5K) use HVGs
N_HVG = 2000
N_PCS = 30
N_NEIGHBORS = 15
LEIDEN_RES = 0.5
SPATIAL_NEIGHS = 6
SEED = 0
# S18 add-on panels carry clonotype probes named <run>_<CDR3>_TRA/TRB; counted apart from ordinary genes.
TCR_PROBE = re.compile(r"^[0-9A-Za-z]+_[A-Z]+_TR[ABDG]$")


class NNDescent(PyNNDescentTransformer):
    # scanpy 1.11.5 calls getnnz(), which scipy>=1.15 csr_array lacks; return csr_matrix, same algorithm.
    def fit_transform(self, X, y=None, **kw):
        return csr_matrix(super().fit_transform(X, y, **kw))

    def transform(self, X):
        return csr_matrix(super().transform(X))


def frozen_genes(path):
    plan = json.loads(Path(path).read_text())
    return sorted({g for genes in plan["panels"].values() for g in genes})


def control_totals(h5, keep_barcodes):
    full = sc.read_10x_h5(h5, gex_only=False)
    full.var_names_make_unique()
    full = full[full.obs_names.isin(keep_barcodes)]
    ftype = full.var["feature_types"].astype(str).values
    mask = np.isin(ftype, NEG_TYPES)
    tot = np.asarray(full[:, mask].X.sum(axis=0)).ravel()
    return pd.DataFrame({"feature": full.var_names[mask], "feature_type": ftype[mask], "total_counts": tot})


def concordance(labels, other):
    ok = other.notna().values
    if ok.sum() == 0:
        return None
    a, b = labels[ok].astype(str), other[ok].astype(str)
    return {"matched_cells": int(ok.sum()), "n_ref_labels": int(b.nunique()),
            "ARI": float(adjusted_rand_score(b, a)), "NMI": float(normalized_mutual_info_score(b, a))}


def main(samples_tsv, row_index, out_root, frozen_plan):
    row = list(csv.DictReader(open(samples_tsv), delimiter="\t"))[int(row_index) - 1]
    sid = row["sample_id"]
    out_root = Path(out_root)
    out = out_root / "stats" / sid
    out.mkdir(parents=True, exist_ok=True)
    sc.settings.n_jobs = 4

    a = sc.read_h5ad(out_root / "h5ad" / f"{sid}.counts.h5ad")
    n_all = a.n_obs
    a = a[a.obs["qc_pass"].values].copy()
    a.layers["counts"] = a.X.copy()

    ctrl = control_totals(row["h5"], set(a.obs_names))
    ctrl.to_csv(out / "negative_control_features.tsv", sep="\t", index=False)
    bg_max, bg_mean = float(ctrl["total_counts"].max()), float(ctrl["total_counts"].mean())

    counts = a.layers["counts"]
    gene_total = np.asarray(counts.sum(axis=0)).ravel()
    genes = pd.DataFrame({
        "gene": a.var_names, "total_counts": gene_total,
        "frac_cells_detected": np.asarray((counts > 0).mean(axis=0)).ravel(),
        "mean_counts": gene_total / a.n_obs,
        "ratio_to_mean_negative": gene_total / bg_mean if bg_mean > 0 else np.inf,
        "above_background": gene_total > bg_max,
        "feature_class": ["tcr_clonotype_probe" if TCR_PROBE.match(g) else "gene" for g in a.var_names],
    }).set_index("gene")
    is_gene = genes["feature_class"] == "gene"

    sc.pp.normalize_total(a)
    sc.pp.log1p(a)
    use_hvg = a.n_vars > HVG_IF_MORE_THAN
    if use_hvg:
        sc.pp.highly_variable_genes(a, flavor="seurat", n_top_genes=N_HVG)  # on log1p data; seurat_v3 needs scikit-misc, not in the locked env
    else:
        a.var["highly_variable"] = True
    a.var["highly_variable"] &= (genes["feature_class"] == "gene").reindex(a.var_names).values  # PCA on ordinary genes
    p = a[:, a.var["highly_variable"].values].copy()
    sc.pp.scale(p, max_value=10)
    n_pcs = min(N_PCS, p.n_vars - 1)
    sc.pp.pca(p, n_comps=n_pcs, random_state=SEED)
    a.obsm["X_pca"] = p.obsm["X_pca"]
    a.uns["pca"] = p.uns["pca"]
    del p
    nn = NNDescent(n_neighbors=N_NEIGHBORS, metric="euclidean", random_state=SEED, n_jobs=4)
    sc.pp.neighbors(a, n_neighbors=N_NEIGHBORS, n_pcs=n_pcs, transformer=nn, random_state=SEED)
    sc.tl.leiden(a, resolution=LEIDEN_RES, flavor="igraph", n_iterations=2, random_state=SEED)
    sc.tl.umap(a, random_state=SEED)

    sc.tl.rank_genes_groups(a, "leiden", method="wilcoxon", pts=True)
    mk = sc.get.rank_genes_groups_df(a, None)
    mk = mk.sort_values(["group", "scores"], ascending=[True, False]).groupby("group", observed=True).head(10)
    mk.to_csv(out / "cluster_markers_top10.tsv", sep="\t", index=False)
    sizes = a.obs["leiden"].value_counts().sort_index()
    sizes.rename("n_cells").to_frame().assign(frac=sizes / a.n_obs).to_csv(out / "cluster_sizes.tsv", sep="\t")

    sq.gr.spatial_neighbors(a, coord_type="generic", n_neighs=SPATIAL_NEIGHS, key_added="spatial")
    moran_genes = list(a.var_names) if not use_hvg else sorted(
        set(a.var_names[a.var["highly_variable"].values]) | (set(frozen_genes(frozen_plan)) & set(a.var_names)))
    sq.gr.spatial_autocorr(a, mode="moran", genes=moran_genes, n_perms=None, n_jobs=4)
    mor = a.uns["moranI"][["I", "pval_norm"]].rename(columns={"I": "moran_I", "pval_norm": "moran_p"})
    # BH over genes with a defined P; squidpy's own BH turns the whole column NaN if any P is NaN (zero-count gene).
    ok = mor["moran_p"].notna() & mor["moran_I"].notna()
    mor["moran_p_bh"] = np.nan
    mor.loc[ok, "moran_p_bh"] = false_discovery_control(mor.loc[ok, "moran_p"].to_numpy(), method="bh")
    genes = genes.join(mor)
    genes.to_csv(out / "gene_stats.tsv", sep="\t")

    stats = {"sample_id": sid, "n_cells_all": int(n_all), "n_cells_qc_pass": int(a.n_obs), "n_features": int(a.n_vars),
             "neg_control_features": int(len(ctrl)), "neg_feature_total_max": bg_max, "neg_feature_total_mean": bg_mean,
             "n_tcr_clonotype_probes": int((~is_gene).sum()),
             "n_genes_above_background": int(genes.loc[is_gene, "above_background"].sum()),
             "n_tcr_probes_above_background": int(genes.loc[~is_gene, "above_background"].sum()),
             "pca_input": f"{N_HVG} seurat HVGs" if use_hvg else "all panel genes (no TCR probes)",
             "n_pca_genes": int(a.var["highly_variable"].sum()), "n_pcs": n_pcs,
             "n_neighbors": N_NEIGHBORS, "leiden_resolution": LEIDEN_RES, "n_leiden_clusters": int(sizes.size),
             "largest_cluster_frac": float(sizes.max() / a.n_obs), "moran_genes": len(moran_genes),
             "moran_I_median": float(mor["moran_I"].median()),
             "n_genes_moran_bh05_I_gt_0.1": int(((mor["moran_p_bh"] < 0.05) & (mor["moran_I"] > 0.1)).sum())}
    for key in ["tenx_graphclust", "author_supervised", "author_unsupervised"]:
        if key in a.obs:
            stats[f"vs_{key}"] = concordance(a.obs["leiden"], a.obs[key])
    (out / "section_stats.json").write_text(json.dumps(stats, indent=2))

    fig, ax = plt.subplots(1, 2, figsize=(14, 6.5))
    sc.pl.umap(a, color="leiden", ax=ax[0], show=False, legend_loc="on data", size=1, title=f"{sid} UMAP")
    sc.pl.embedding(a, "spatial", color="leiden", ax=ax[1], show=False, legend_loc=None, size=1, title="centroids (µm)")
    ax[1].set_aspect("equal")
    fig.savefig(out / "leiden_umap_spatial.png", dpi=110, bbox_inches="tight")
    for c in a.obs.columns:
        if a.obs[c].dtype == object:
            a.obs[c] = a.obs[c].astype(str)
    a.write_h5ad(out_root / "h5ad" / f"{sid}.norm.h5ad", compression="gzip")
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main(*sys.argv[1:5])
