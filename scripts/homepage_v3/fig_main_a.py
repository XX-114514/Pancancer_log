"""Main figures 1-3: concept, data organisation, processing and annotation."""
from __future__ import annotations

import numpy as np
import pandas as pd
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch

from style import (REPO, SRC2, INK, MUTED, GRID, CM_COL, LINEAGE_COL, SEQ, DIV, panel, fig_letter, box, arrow,
                   export, save)


# ---------------------------------------------------------------------------
# Figure 1  concept: hubs -> ecotypes -> cellular modules; CoVarNet; questions; design
# ---------------------------------------------------------------------------
def fig1() -> None:
    rng = np.random.default_rng(7)
    fig = plt.figure(figsize=(7.2, 6.4))

    # ---- a: three frameworks -------------------------------------------------
    titles = [("Spatial immune hub", "Pelka et al., Cell 2021"),
              ("Cell-state ecotype", "Luca et al., Cell 2021; Spatial EcoTyper 2026"),
              ("Cellular module (CM)", "CoVarNet, Nature 2025")]
    for i, (t, ref) in enumerate(titles):
        ax = fig.add_axes([0.02 + i * 0.33, 0.66, 0.30, 0.28])
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
        ax.add_patch(FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0,rounding_size=0.03",
                                    fc="#fafbfc", ec="#d5dbe3", lw=0.6, zorder=0))
        ax.text(0.04, 0.93, t, fontsize=7, fontweight="bold", color=INK, va="top")
        ax.text(0.04, 0.84, ref, fontsize=5.2, color=MUTED, va="top", style="italic")
        if i == 0:
            # tissue with a hub of CXCL13+ T, myeloid and malignant cells
            for _ in range(70):
                x, y = rng.uniform(0.08, 0.92), rng.uniform(0.08, 0.72)
                ax.add_patch(Circle((x, y), 0.016, fc="#ead3dc", ec="none"))
            cx, cy = 0.52, 0.38
            ax.add_patch(Circle((cx, cy), 0.2, fc="none", ec="#c4502a", lw=0.8, ls=(0, (2, 1.5))))
            pts = rng.normal(0, 0.08, (26, 2))
            cols = ["#3b6ea8"] * 10 + ["#c2532b"] * 9 + ["#d16d8e"] * 7
            for (dx, dy), c in zip(pts, cols):
                ax.add_patch(Circle((cx + dx, cy + dy * 0.9), 0.018, fc=c, ec="white", lw=0.3))
            for k, (lab, c) in enumerate([("CXCL13⁺ T", "#3b6ea8"), ("myeloid", "#c2532b"), ("malignant", "#d16d8e")]):
                ax.add_patch(Circle((0.08, 0.17 - k * 0.065), 0.016, fc=c, ec="none"))
                ax.text(0.12, 0.17 - k * 0.065, lab, fontsize=5.2, va="center", color=INK)
            ax.text(0.76, 0.12, "co-localised\nmulticellular\nprogramme", fontsize=5, color="#c4502a", ha="center")
        elif i == 1:
            # cell states per type -> ecotypes (block structure)
            M = rng.uniform(0, 0.25, (9, 12))
            for b, (r0, c0) in enumerate([(0, 0), (3, 4), (6, 8)]):
                M[r0:r0 + 3, c0:c0 + 4] += rng.uniform(0.6, 1.0, (3, 4))
            ax.imshow(M, extent=(0.1, 0.72, 0.1, 0.72), cmap=SEQ, aspect="auto", vmin=0, vmax=1.1, zorder=3)
            ax.text(0.41, 0.04, "samples", fontsize=5.2, ha="center", color=MUTED)
            ax.text(0.06, 0.41, "cell states", fontsize=5.2, rotation=90, va="center", color=MUTED)
            for k, (lab, y) in enumerate([("E1", 0.62), ("E2", 0.41), ("E3", 0.2)]):
                ax.text(0.8, y, lab, fontsize=6, fontweight="bold", color=["#1b7f5f", "#c4502a", "#2f63a3"][k], va="center")
            ax.text(0.86, 0.73, "states that\nco-occur", fontsize=5, color=MUTED, ha="center", va="bottom")
        else:
            # sample x subcluster frequency -> W (subclusters x K) network
            M = rng.uniform(0, 0.3, (8, 10))
            M[:3, :5] += 0.6; M[3:6, 4:] += 0.6; M[6:, ::2] += 0.5
            ax.imshow(M, extent=(0.06, 0.40, 0.1, 0.70), cmap=SEQ, aspect="auto", vmin=0, vmax=1, zorder=3)
            ax.text(0.23, 0.04, "samples", fontsize=5.2, ha="center", color=MUTED)
            ax.text(0.025, 0.40, "subclusters", fontsize=5.2, rotation=90, va="center", color=MUTED)
            ax.annotate("", xy=(0.52, 0.40), xytext=(0.43, 0.40),
                        arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.8))
            ax.text(0.475, 0.47, "nsNMF", fontsize=5.2, ha="center", color=MUTED)
            # module networks
            centers = [(0.66, 0.56, "#1b7f5f"), (0.84, 0.56, "#c4502a"), (0.75, 0.24, "#2f63a3")]
            for (x, y, c) in centers:
                ang = np.linspace(0, 2 * np.pi, 5)[:-1] + 0.4
                pts = [(x + 0.07 * np.cos(a), y + 0.08 * np.sin(a)) for a in ang]
                for p in range(4):
                    for q in range(p + 1, 4):
                        if (p + q) % 2 or p == 0:
                            ax.plot([pts[p][0], pts[q][0]], [pts[p][1], pts[q][1]], color=c, lw=0.6, alpha=0.7)
                for (px, py) in pts:
                    ax.add_patch(Circle((px, py), 0.022, fc=c, ec="white", lw=0.4))
            ax.text(0.75, 0.72, "co-varying cell\nmodules + networks", fontsize=5, ha="center", color=MUTED)
    fig_letter(fig, 0.005, 0.95, "a")
    for i in range(2):
        fig.text(0.325 + i * 0.33, 0.80, "→", fontsize=12, color="#9aa5b1", ha="center", va="center")
    fig.text(0.02, 0.635, "From co-localised hubs (space) to co-occurring states (ecotypes) to co-varying composition (cellular modules)",
             fontsize=5.8, color=MUTED)

    # ---- b: CoVarNet algorithm ----------------------------------------------
    ax = fig.add_axes([0.02, 0.36, 0.55, 0.25]); ax.axis("off")
    fig_letter(fig, 0.005, 0.605, "b")
    ax.text(0.0, 0.98, "CoVarNet in four steps (as implemented for K = 9, 0906 version)", fontsize=7, fontweight="medium", va="top")
    steps = [("1  Frequency", "subcluster ÷ parent\nlineage per sample;\nmissing = 0;\nmin–max per row", "#eef3f8"),
             ("2  nsNMF", "X ≈ W·S·H\nS = (1−θ)I + (θ/K)J\nθ = 0.5; 30 starts\nper K", "#eef3f8"),
             ("3  Rank, members", "K by cophenetic\ncurve; top-N rows\nof W = members", "#eef3f8"),
             ("4  Network", "member–member\ncross-sample\ncorrelation; H =\nsample activity", "#eef3f8")]
    for k, (h, t, c) in enumerate(steps):
        x = 0.0 + k * 0.25
        box(ax, x, 0.14, 0.225, 0.68, "", fc=c, ec="#c5d3e3")
        ax.text(x + 0.015, 0.76, h, fontsize=6.0, fontweight="bold", color="#2f4b6b", va="top", transform=ax.transAxes)
        ax.text(x + 0.015, 0.60, t, fontsize=5.4, color=INK, va="top", transform=ax.transAxes, linespacing=1.3)
        if k < 3:
            arrow(ax, x + 0.227, 0.48, x + 0.249, 0.48, color="#7d8a99")
    ax.text(0.0, 0.06, "Input semantics are composition (not expression); scaling, rank heuristic and top-N choice shape what a CM can mean (Subpage 4).",
            fontsize=5.4, color=MUTED, transform=ax.transAxes, va="top")

    # ---- c: two questions ----------------------------------------------------
    ax = fig.add_axes([0.60, 0.36, 0.39, 0.25]); ax.axis("off")
    fig_letter(fig, 0.59, 0.605, "c")
    ax.text(0.0, 0.98, "Questions asked in this project", fontsize=7, fontweight="medium", va="top")
    box(ax, 0.0, 0.50, 1.0, 0.38,
        "Q1  Do CMs matter in tumour development and therapy?\n"
        "Is any CM higher in immunotherapy cohorts or associated with\nfavourable / unfavourable outcome, and in which cancer types?",
        fc="#fdf1ec", ec="#e8b9a5", fs=5.7, ha="left")
    box(ax, 0.0, 0.06, 1.0, 0.38,
        "Q2  Is intra-CM communication linked to immune infiltration?\n"
        "Tumour–immune, myeloid–stromal ligand–receptor axes inside a CM:\ndo they align with T-cell infiltration or exclusion in tissue?",
        fc="#eef3fa", ec="#b5c8e0", fs=5.7, ha="left")

    # ---- d: study design -----------------------------------------------------
    ax = fig.add_axes([0.02, 0.02, 0.97, 0.30]); ax.axis("off")
    fig_letter(fig, 0.005, 0.315, "d")
    ax.text(0.0, 0.98, "Study design shown on this page (0906 version) and the follow-up versions", fontsize=7, fontweight="medium", va="top")
    row1 = [("Public scRNA-seq\n54 GSE · 43 labels\n1,322 samples · 4.68 M cells", "#f3f5f8"),
            ("QC · merge · BBKNN\nlayered annotation (V5)\nper-sample CNV (V7)", "#f3f5f8"),
            ("CoVarNet K = 9\n1,303 samples\nCM01 / CM07 / CM08", "#eaf2ee"),
            ("Within-CM LIANA\n1,301 samples\n12 focus axes", "#eaf2ee")]
    row2 = [("Visium\n5 datasets · 61 sections", "#fbf1e9"),
            ("CosMx ccRCC\n3 tissues · 199,112 cells", "#fbf1e9"),
            ("TCGA\n9,034 tumours · 27 types", "#fbf1e9"),
            ("ICI bulk cohorts\nRiaz · Hugo · IMvigor210\nRose · Jung", "#fbf1e9")]
    w, h = 0.215, 0.30
    for k, (t, c) in enumerate(row1):
        x = k * 0.245
        box(ax, x, 0.52, w, h, t, fc=c, fs=5.7)
        if k < 3:
            arrow(ax, x + w + 0.003, 0.67, x + 0.243, 0.67)
    for k, (t, c) in enumerate(row2):
        x = k * 0.245
        box(ax, x, 0.08, w, 0.28, t, fc=c, fs=5.7)
    ax.text(0.0, 0.84, "discovery (sample = unit of replication)", fontsize=5.2, color=MUTED)
    ax.text(0.0, -0.01, "independent readouts (associations; not proof of communication or causality)", fontsize=5.2, color=MUTED)
    ax.plot([0.1075, 0.8425], [0.445, 0.445], color="#9aa5b1", lw=0.7, transform=ax.transAxes, clip_on=False)
    ax.plot([0.6975, 0.6975], [0.52, 0.445], color="#9aa5b1", lw=0.7, transform=ax.transAxes, clip_on=False)
    for k in range(4):
        arrow(ax, 0.1075 + k * 0.245, 0.445, 0.1075 + k * 0.245, 0.37, color="#9aa5b1")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.text(0.985, 0.012, "Follow-up discovery versions (Section 9): V4 K14 (2026-11-06) · WEHI R4 K15 (2026-11-15)",
             fontsize=5.2, color=MUTED, ha="right")
    save(fig, "fig1_concept_design", ["schematic; numbers from README text and releases/communication_evidence_chain_20260906"])


# ---------------------------------------------------------------------------
# Figure 2  data organisation
# ---------------------------------------------------------------------------
GENERATIONS = pd.DataFrame([
    ("2026-07", "Project_v3 discovery", 1963745, 773, "31 GSE", "first integration + pan-cancer CoVarNet"),
    ("2026-08/09", "5M merged (0906 basis)", 4676787, 1322, "54 GSE · 43 labels", "K=9 CMs, LIANA, V5→V7 identity"),
    ("2026-10-14", "Immune atlas v1", 2670057, 400, "9 sources (libraries)", "ICI / immune-enriched; T-cell states"),
    ("2026-10-23", "Full-context archive", 9439479, None, "65 source keys", "unified identity, crosswalk"),
    ("2026-10-30", "Eligible TME scope", 2968368, 770, "43 sources (units)", "non-enriched TME composition + states"),
    ("2026-11-06", "V4 non-epithelial K14", 560871, 164, "13 sources", "author-configured CoVarNet replication"),
    ("2026-11-15", "WEHI R4 K15", 1974301, 1194, "R4 units", "vascular CM12 validation"),
], columns=["date", "version", "cells", "units", "sources", "role"])


def fig2() -> None:
    cc = pd.read_csv(SRC2 / "fig2_cancer_composition.tsv", sep="\t").sort_values("n_cells")
    sc = pd.read_csv(SRC2 / "fig2_sample_class_composition.tsv", sep="\t")
    export(GENERATIONS, "fig2_dataset_generations.tsv")

    fig = plt.figure(figsize=(7.2, 6.2))
    gs = fig.add_gridspec(3, 3, width_ratios=[1.35, 1, 1], height_ratios=[1, 1, 0.95], hspace=0.95, wspace=0.65)

    ax = fig.add_subplot(gs[:2, 0])
    cols = {1: "#c9d6e4", 2: "#8fb0cf", 3: "#5f8fbf"}
    c = [cols.get(min(g, 3), "#1d3f66") if g < 4 else "#1d3f66" for g in cc.n_gse]
    y = np.arange(len(cc))
    ax.barh(y, cc.n_cells / 1e3, color=c, height=0.72, lw=0)
    for yi, (n, s) in enumerate(zip(cc.n_cells, cc.n_samples)):
        ax.text(n / 1e3 + 8, yi, f"{s}", va="center", fontsize=4.8, color=MUTED)
    ax.set_yticks(y, cc.cancer_group, fontsize=5.0)
    ax.set_xlabel("Cells (×10³)  · number = samples")
    ax.set_ylim(-0.6, len(cc) - 0.4)
    ax.set_title("Cancer groups in the 0906 discovery set", fontsize=6.8)
    ax.legend([Rectangle((0, 0), 1, 1, fc=c_) for c_ in [cols[1], cols[2], cols[3], "#1d3f66"]],
              ["1", "2", "3", "≥4"], title="GSEs per group", title_fontsize=5.2, ncol=4, loc="lower right",
              handlelength=0.9, columnspacing=0.8, handletextpad=0.3, bbox_to_anchor=(1.0, 0.03))
    ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)
    panel(ax, "a", x=-0.55)

    # b: two batches
    ax = fig.add_subplot(gs[0, 1])
    batches = pd.DataFrame({"batch": ["July cohort\n(31 GSE)", "August cohort\n(23 GSE)"],
                            "cells": [1963745, 2713042], "samples": [773, 549]})
    export(batches, "fig2_batches.tsv")
    ax.bar([0, 1], batches.cells / 1e6, color=["#8fb0cf", "#e3a38c"], width=0.6)
    for i, r in batches.iterrows():
        ax.text(i, r.cells / 1e6 + 0.08, f"{r.cells/1e6:.2f} M\n{r.samples} samples", ha="center", fontsize=5.2)
    ax.set_xticks([0, 1], batches.batch, fontsize=5.4)
    ax.set_ylabel("Cells (M)"); ax.set_ylim(0, 3.4)
    ax.set_title("Two merged batches\n(10,432 shared genes)", fontsize=6.6)
    panel(ax, "b", x=-0.35)

    # c: QC of the August batch
    ax = fig.add_subplot(gs[0, 2])
    vals = [3282991, 2713042]
    ax.bar([0, 1], np.array(vals) / 1e6, color=["#c9ced6", "#e3a38c"], width=0.6)
    ax.text(0, vals[0] / 1e6 + 0.08, "3.28 M", ha="center", fontsize=5.4)
    ax.text(1, vals[1] / 1e6 + 0.08, "2.71 M\n(82.6%)", ha="center", fontsize=5.4)
    ax.set_xticks([0, 1], ["before QC", "after QC +\nScrublet"], fontsize=5.4)
    ax.set_ylim(0, 4.9); ax.set_ylabel("Cells (M)")
    ax.set_title("Per-sample adaptive QC\n(August batch)", fontsize=6.6)
    ax.text(0.5, 0.99, "13,821 doublets removed\n1 GSE excluded (log1p .X)", transform=ax.transAxes, ha="center",
            fontsize=4.9, color=MUTED, va="top")
    panel(ax, "c", x=-0.35)

    # d: sample material
    ax = fig.add_subplot(gs[1, 1:])
    order = sc.sort_values("n_samples", ascending=True)
    colmap = {"Primary/unspecified tumor": "#c4502a", "Metastatic tumor": "#8c2f18", "Normal/adjacent tissue": "#5fa58a",
              "Blood": "#b23a24", "Bone marrow": "#795548", "Lymph node": "#7b5ea7", "Unknown": "#c9ced6"}
    left = 0
    for _, r in sc.iterrows():
        ax.barh(1, r.sample_percent, left=left, color=colmap[r.sample_class_plot], height=0.55, lw=0.4, ec="white")
        ax.barh(0, r.cell_percent, left=left if False else None, color="none")
        left += r.sample_percent
    left = 0
    for _, r in sc.iterrows():
        ax.barh(0, r.cell_percent, left=left, color=colmap[r.sample_class_plot], height=0.55, lw=0.4, ec="white")
        if r.cell_percent > 5:
            ax.text(left + r.cell_percent / 2, 0, f"{r.cell_percent:.0f}%", ha="center", va="center", fontsize=5, color="white")
        left += r.cell_percent
    left = 0
    for _, r in sc.iterrows():
        if r.sample_percent > 5:
            ax.text(left + r.sample_percent / 2, 1, f"{r.n_samples}", ha="center", va="center", fontsize=5, color="white")
        left += r.sample_percent
    ax.set_yticks([0, 1], ["cells", "samples"]); ax.set_xlim(0, 100); ax.set_xlabel("%")
    ax.set_title("Sample material (numbers on sample bar = samples)", fontsize=6.6, pad=22)
    handles = [Rectangle((0, 0), 1, 1, fc=v) for v in colmap.values()]
    ax.legend(handles, colmap.keys(), ncol=4, bbox_to_anchor=(0.0, 1.0), loc="lower left", fontsize=5.2,
              handlelength=1, columnspacing=1)
    panel(ax, "d", x=-0.1)

    # e: metadata availability
    ax = fig.add_subplot(gs[2, 0])
    meta = pd.DataFrame({"field": ["Immunotherapy status", "Metastasis status", "Material type", "Assay (GEO text)"],
                         "known": [127, 1322 - 968, 1322 - 284, 1322 - 850]})
    meta["unknown"] = 1322 - meta.known
    export(meta, "fig2_metadata_availability.tsv")
    y = np.arange(len(meta))[::-1]
    ax.barh(y, meta.known / 1322 * 100, color="#5f8fbf", height=0.6)
    ax.barh(y, meta.unknown / 1322 * 100, left=meta.known / 1322 * 100, color="#e3e7ec", height=0.6)
    for yi, k in zip(y, meta.known):
        ax.text(1.5, yi, f"{k} known", va="center", fontsize=5, color="white" if k > 200 else INK)
    ax.set_yticks(y, meta.field, fontsize=5.6); ax.set_xlim(0, 100); ax.set_xlabel("% of 1,322 samples")
    ax.set_title("Clinical metadata readable from GEO", fontsize=6.6)
    panel(ax, "e", x=-0.55)

    # f: dataset generations
    ax = fig.add_subplot(gs[2, 1:])
    g = GENERATIONS.iloc[::-1].reset_index(drop=True)
    colr = {"5M merged (0906 basis)": "#c4502a", "Immune atlas v1": "#7b5ea7", "Eligible TME scope": "#2a7f62",
            "V4 non-epithelial K14": "#2f63a3", "WEHI R4 K15": "#2f63a3"}
    cl = [colr.get(v, "#9aa5b1") for v in g.version]
    ax.barh(range(len(g)), g.cells / 1e6, color=cl, height=0.62)
    for i, r in g.iterrows():
        ax.text(r.cells / 1e6 + 0.15, i, f"{r.date} · {r.cells/1e6:.2f} M · {r.sources} — {r.role}", va="center", fontsize=4.7, color=INK)
    ax.set_yticks(range(len(g)), list(g.version), fontsize=4.9)
    ax.set_xlim(0, 24); ax.set_xticks([0, 2.5, 5, 7.5, 10]); ax.set_xlabel("Cells / observations (M)")
    ax.set_title("Dataset generations (orange = this page's main line; blue = follow-up CM versions)", fontsize=6.4)
    for t, c_ in zip(ax.get_yticklabels(), cl):
        t.set_color(c_ if c_ != "#9aa5b1" else INK)
    panel(ax, "f", x=-0.03, y=1.06)
    save(fig, "fig2_data_organisation", ["inventories/publication/2026-11-02/homepage_v2/fig2_*.tsv", "docs/homepage/01_datasets.md"])


# ---------------------------------------------------------------------------
# Figure 3  processing pipeline, integration choice, UMAP and annotation
# ---------------------------------------------------------------------------
def _crop(png, box_frac):
    im = Image.open(png)
    W, H = im.size
    l, t, r, b = box_frac
    return np.asarray(im.crop((int(l * W), int(t * H), int(r * W), int(b * H))))


def fig3() -> None:
    old = REPO / "assets/publication/20261102/fig3_atlas_umap_annotation.png"
    vc = pd.read_csv(SRC2 / "fig3_v7_major_counts.tsv", sep="\t")
    tr = pd.read_csv(SRC2 / "fig3_v5_to_v7_unknown_resolution.tsv", sep="\t")

    fig = plt.figure(figsize=(7.2, 7.4))
    # a pipeline
    ax = fig.add_axes([0.02, 0.80, 0.96, 0.17]); ax.axis("off")
    fig_letter(fig, 0.005, 0.965, "a")
    ax.text(0, 1.0, "Processing pipeline of the 0906 version", fontsize=7, fontweight="medium", va="top")
    steps = [("Per-sample QC", "strict / floor MAD\nthresholds; Scrublet"),
             ("Counts check", "non-negative integer\n.X → layers['counts']"),
             ("Merge", "gene intersection\n10,432 genes"),
             ("HVG / PCA", "4,000 HVG (seurat_v3)\n50 PCs; batch = GSE"),
             ("Integration", "Harmony · BBKNN · scVI\nscIB comparison"),
             ("Major identity", "author labels + marker\nscores + CellTypist"),
             ("Fine identity", "lineage Leiden +\nextended panels"),
             ("Malignancy", "per-sample inferCNVpy\n+ CopyKAT tiers")]
    w = 0.112
    for k, (h, t) in enumerate(steps):
        x = k * 0.125
        fc = "#eaf2ee" if k >= 5 else "#f3f5f8"
        box(ax, x, 0.08, w, 0.66, "", fc=fc, ec="#c5ccd6")
        ax.text(x + w / 2, 0.62, h, ha="center", fontsize=5.9, fontweight="bold", color="#2f4b6b", transform=ax.transAxes)
        ax.text(x + w / 2, 0.33, t, ha="center", va="center", fontsize=4.9, color=INK, transform=ax.transAxes)
        if k < 7:
            arrow(ax, x + w + 0.002, 0.41, x + 0.124, 0.41)

    # b scIB scores
    ax = fig.add_axes([0.07, 0.56, 0.22, 0.17])
    sc = pd.DataFrame({"method": ["BBKNN", "scVI", "Harmony"], "overall": [0.4917, 0.4862, 0.4807]})
    export(sc, "fig3_scib_overall.tsv")
    ax.barh([2, 1, 0], sc.overall, color=["#c4502a", "#9aa5b1", "#9aa5b1"], height=0.55)
    for yi, v in zip([2, 1, 0], sc.overall):
        ax.text(v + 0.0006, yi, f"{v:.4f}", va="center", fontsize=5.4)
    ax.set_yticks([2, 1, 0], sc.method); ax.set_xlim(0.475, 0.497)
    ax.set_xlabel("scIB overall\n(0.4 batch + 0.6 bio)")
    ax.set_title("Integration comparison", fontsize=6.6)
    ax.text(0.0, -0.62, "ranks swap at batch weight ≈ 0.317;\nBBKNN graph used downstream", transform=ax.transAxes,
            fontsize=5, color=MUTED, va="top")
    panel(ax, "b", x=-0.42)

    # c UMAP (re-used raster)
    ax = fig.add_axes([0.36, 0.43, 0.62, 0.33])
    ax.imshow(_crop(old, (0.005, 0.058, 0.60, 0.485))); ax.axis("off")
    ax.set_title("Major lineage on the BBKNN UMAP (V5 identity = CoVarNet input; 200,000 of 4,676,787 cells)", fontsize=6.4)
    panel(ax, "c", x=-0.05, y=1.0)

    # d, e: CNV tier & confidence UMAP
    ax = fig.add_axes([0.02, 0.06, 0.30, 0.34])
    ax.imshow(_crop(old, (0.035, 0.555, 0.335, 0.955))); ax.axis("off")
    ax.set_title("Per-sample CNV malignancy tier (V7)", fontsize=6.4)
    panel(ax, "d", x=-0.1, y=1.0)
    ax = fig.add_axes([0.33, 0.06, 0.30, 0.34])
    ax.imshow(_crop(old, (0.37, 0.555, 0.61, 0.955))); ax.axis("off")
    ax.set_title("Fine-identity confidence (V7)", fontsize=6.4)
    panel(ax, "e", x=-0.1, y=1.0)

    # f V7 major counts
    ax = fig.add_axes([0.73, 0.235, 0.25, 0.15])
    top = vc.head(14).iloc[::-1]
    cols = [LINEAGE_COL.get(m, "#b9c2cc") for m in top.major_celltype_v7]
    ax.barh(range(len(top)), top.n_cells / 1e6, color=cols, height=0.7)
    ax.set_yticks(range(len(top)), [m.replace("_", " ") for m in top.major_celltype_v7], fontsize=4.9)
    ax.set_xlabel("Cells (M)"); ax.set_title("V7 major identities (top 14 of 26)", fontsize=6.4)
    panel(ax, "f", x=-0.38)
    # g V5 -> V7 fate
    ax = fig.add_axes([0.73, 0.075, 0.25, 0.085])
    cats = ["Resolved lineage", "Boundary class", "Still Unknown/Ambig."]
    colz = ["#2f63a3", "#9cbad6", "#a7afba"]
    for i, r in enumerate(tr.itertuples()):
        left = 0
        for c, cl in zip(cats, colz):
            v = getattr(r, c.replace(" ", "_").replace("/", "_").replace(".", "_"), None) if False else tr.iloc[i][c]
            ax.barh(i, v * 100, left=left, color=cl, height=0.6)
            if v > 0.12:
                ax.text(left + v * 50, i, f"{v*100:.0f}%", ha="center", va="center", fontsize=5, color="white" if cl != "#a7afba" else INK)
            left += v * 100
    ax.set_yticks([0, 1], [f"V5 {v}" for v in tr.v5_major], fontsize=5.4)
    ax.set_xlim(0, 100); ax.set_xlabel("% re-assigned in V7"); ax.set_title("V5 uncertain cells in V7", fontsize=6.4)
    ax.legend([Rectangle((0, 0), 1, 1, fc=c) for c in colz], ["resolved", "boundary", "still unresolved"],
              ncol=3, bbox_to_anchor=(-0.5, -0.62), loc="upper left", fontsize=5, handlelength=1, columnspacing=0.8)
    panel(ax, "g", x=-0.38, y=1.15)
    save(fig, "fig3_processing_annotation", ["assets/publication/20261102/fig3_atlas_umap_annotation.png (UMAP raster)",
                                              "inventories/publication/2026-11-02/homepage_v2/fig3_*.tsv"])
