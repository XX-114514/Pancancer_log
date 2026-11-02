#!/usr/bin/env python3
"""Render the 2026-11-02 paper-style homepage figures (homepage v2).

Inputs are read-only copies of aggregate tables exported from the research
server (V5 CoVarNet K=9 / LIANA run 20260825_v5_final_v1, the 2026-09-06
communication evidence chain, the 2026-09-03 Visium expansion and the
2026-09-04 bulk cohorts) plus the V7 release tables already in this repository.

Usage
-----
    RAW=/path/to/pulled/server/tables python3 scripts/render_homepage_v2_figures.py

``RAW`` must contain the sub-directories ``v5final/``, ``pv5/``, ``cmp0829/`` and
``p5m/`` with the same relative layout as on the server. Cell-level UMAP
coordinates are only read to draw Fig. 3 and are never written to the
repository. Every other figure is drawn from sample- or cohort-level
aggregates, which are re-exported to ``inventories/publication/2026-11-02/
homepage_v2/`` so that the plotted numbers can be audited.
"""
from __future__ import annotations

import json
import os
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import patches  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
RAW = Path(os.environ.get("RAW", ""))
OUT = REPO / "assets" / "publication" / "20261102"
AGG = REPO / "inventories" / "publication" / "2026-11-02" / "homepage_v2"
V7 = REPO / "releases" / "annotation_v7_20260905" / "tables"
CHAIN = REPO / "releases" / "communication_evidence_chain_20260906" / "tables"

V5F = RAW / "v5final"
PV5 = RAW / "pv5"
CMP = RAW / "cmp0829" / "tables"
P5M = RAW / "p5m"
SPAT = PV5 / "spatial_expansion_20260903" / "runs" / "20260904_all_visium_cm_region_communication_v1" / "tables"
BULK = PV5 / "bulk_cohorts_20260904"
ECH = PV5 / "communication_evidence_chain_20260906" / "tables"

# --------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------
plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 7,
        "axes.titlesize": 7.5,
        "axes.labelsize": 7,
        "xtick.labelsize": 6.3,
        "ytick.labelsize": 6.3,
        "legend.fontsize": 6.2,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "figure.dpi": 150,
        "savefig.dpi": 300,
    }
)

INK = "#1f2933"
GREY = "#8a94a3"
LIGHT = "#e7eaee"
CM_COL = {
    "CM01": "#2a7f62",
    "CM07": "#c2532b",
    "CM08": "#3b6ea8",
}
OTHER_CM = "#9aa5b1"
DIV = LinearSegmentedColormap.from_list("div", ["#2b5c8a", "#f7f7f7", "#b2402f"])
SEQ = LinearSegmentedColormap.from_list("seq", ["#f7f8fa", "#9cc3d5", "#2c6e8f", "#14324a"])
SEQ_WARM = LinearSegmentedColormap.from_list("seqw", ["#fbf7f2", "#efc5a3", "#c76b3c", "#6d2a12"])

LINEAGE_COL = {
    "Tcell": "#3b6ea8",
    "NK": "#6a9fd4",
    "ILC": "#a3c4e8",
    "Bcell": "#7b5ea7",
    "Plasma": "#b39ddb",
    "Macrophage": "#c2532b",
    "Monocyte": "#e08a5c",
    "Dendritic": "#f0b48f",
    "Mast": "#d4a017",
    "Granulocyte": "#b8860b",
    "Platelet": "#8d6e63",
    "Erythroid": "#a1887f",
    "HSPC": "#795548",
    "Fibroblast": "#2a7f62",
    "Endothelial": "#5fa58a",
    "Epithelial": "#d16d8e",
    "Neural": "#9e7b9b",
    "Unknown": "#c9ced6",
    "Ambiguous": "#a7afba",
}


def lineage_of(label: str) -> str:
    for key in LINEAGE_COL:
        if label.startswith(key) or f"::{key}" in label:
            return key
    if label.startswith("Malignan"):
        return "Epithelial"
    return "Unknown"


def panel(ax, letter: str, x: float = -0.12, y: float = 1.06) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left", color=INK)


def save(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", name)


def export(df: pd.DataFrame, name: str) -> None:
    AGG.mkdir(parents=True, exist_ok=True)
    df.to_csv(AGG / name, sep="\t", index=False)


def box(ax, x, y, w, h, text, fc="#ffffff", ec=INK, fs=6.4, bold=False, color=INK, lw=0.7, r=0.02):
    ax.add_patch(
        patches.FancyBboxPatch(
            (x, y), w, h, boxstyle=f"round,pad=0.004,rounding_size={r}", fc=fc, ec=ec, lw=lw, transform=ax.transAxes
        )
    )
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fs,
        fontweight="bold" if bold else "normal",
        color=color,
        transform=ax.transAxes,
        linespacing=1.25,
    )


def arrow(ax, x0, y0, x1, y1, color=INK, lw=0.8, style="-|>", ms=6, ls="-", rad=0.0):
    ax.annotate(
        "",
        xy=(x1, y1),
        xytext=(x0, y0),
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle=style, color=color, lw=lw, mutation_scale=ms, linestyle=ls,
                        connectionstyle=f"arc3,rad={rad}"),
    )


def short_state(label: str) -> str:
    rep = {
        "Monocyte_Classical_FCN1_FCN1": "FCN1 Mono",
        "Macrophage_SPP1_TREM2_SPP1": "SPP1 Mac",
        "Endothelial_Blood_PECAM1": "Endo (blood)",
        "Fibroblast_myCAF_ACTA2": "myCAF",
        "Malignant::Epithelial": "Malignant epi",
        "Monocyte_Inflammatory_IL1B_IL1B": "Infl. Mono",
        "Macrophage_Inflammatory_IL1B_IL1B": "Infl. Mac",
        "Tcell_CD8_cytotoxic_CCL5": "CD8 cytotoxic",
    }
    return rep.get(label, label.replace("_unresolved", " (unres.)").replace("_", " "))


# --------------------------------------------------------------------------
# Fig 1  concept and workflow
# --------------------------------------------------------------------------
def fig1() -> None:
    fig = plt.figure(figsize=(7.2, 6.0))
    ax = fig.add_axes([0, 0.52, 1, 0.48])
    ax.axis("off")
    panel(ax, "a", x=0.0, y=0.95)
    ax.text(0.03, 0.965, "From cellular hubs to cellular modules", fontsize=7.5, fontweight="bold", transform=ax.transAxes)

    # three conceptual tiles
    tiles = [
        (0.02, "Co-occurring cell hubs",
         "Pelka et al., Cell 2021\nprogramme co-variation across\npatients; spatially co-localised\nimmune hubs (CXCL13+ T / myeloid)"),
        (0.35, "Cellular ecosystems",
         "EcoTyper, Cell 2021;\nSpatial EcoTyper, Nature 2026\ncell states → recurrent\nmulticellular ecotypes"),
        (0.68, "Cellular modules (CMs)",
         "CoVarNet, Nature 2025\nsubtype frequencies across samples\n→ nsNMF → co-varying modules\n+ member networks"),
    ]
    rng = np.random.default_rng(3)
    for i, (x0, title, sub) in enumerate(tiles):
        ax.add_patch(patches.FancyBboxPatch((x0, 0.06), 0.30, 0.82, boxstyle="round,pad=0.006,rounding_size=0.02",
                                            fc="#fafbfc", ec=LIGHT, lw=0.8, transform=ax.transAxes))
        ax.text(x0 + 0.15, 0.82, title, ha="center", fontsize=7, fontweight="bold", transform=ax.transAxes)
        ax.text(x0 + 0.15, 0.13, sub, ha="center", va="bottom", fontsize=5.8, color="#475260", transform=ax.transAxes,
                linespacing=1.3)
        # mini illustration
        ins = ax.inset_axes([x0 + 0.04, 0.38, 0.22, 0.38])
        ins.set_xticks([]); ins.set_yticks([])
        for s in ins.spines.values():
            s.set_visible(False)
        ins.set_xlim(0, 1); ins.set_ylim(0, 1)
        if i == 0:
            cols = ["#3b6ea8", "#c2532b", "#2a7f62", "#d16d8e"]
            for k in range(70):
                c = cols[k % 4]
                cx, cy = (0.32, 0.55) if k % 4 in (0, 1) and k < 40 else (rng.uniform(0.05, 0.95), rng.uniform(0.05, 0.95))
                if k < 40 and k % 4 in (0, 1):
                    cx += rng.normal(0, 0.1); cy += rng.normal(0, 0.1)
                ins.scatter(cx, cy, s=9, color=c, lw=0)
            ins.add_patch(patches.Circle((0.32, 0.55), 0.24, fill=False, ls="--", lw=0.7, ec=INK))
            ins.text(0.62, 0.86, "hub", fontsize=6, color=INK)
        elif i == 1:
            for j, (cx, cy, col) in enumerate([(0.25, 0.7, "#3b6ea8"), (0.7, 0.72, "#c2532b"), (0.5, 0.3, "#2a7f62")]):
                for _ in range(16):
                    ins.scatter(cx + rng.normal(0, 0.07), cy + rng.normal(0, 0.07), s=8, color=col, lw=0)
            ins.plot([0.25, 0.7, 0.5, 0.25], [0.7, 0.72, 0.3, 0.7], color=GREY, lw=0.7, zorder=0)
            ins.text(0.5, 0.02, "ecotype = state co-occurrence", fontsize=5.3, ha="center", color="#475260")
        else:
            m = rng.random((6, 8)) ** 2
            m[0:3, 0:4] += 0.9
            m[3:6, 4:8] += 0.9
            ins.imshow(m, cmap=SEQ, aspect="auto", extent=(0, 0.62, 0.05, 0.95))
            ins.text(0.31, -0.06, "samples", ha="center", fontsize=5.3, color="#475260")
            ins.text(-0.06, 0.5, "subtypes", rotation=90, va="center", fontsize=5.3, color="#475260")
            for k, (yy, col) in enumerate([(0.72, "#c2532b"), (0.28, "#3b6ea8")]):
                ins.add_patch(patches.Circle((0.82, yy), 0.11, fc=col, alpha=0.18, ec=col, lw=0.8))
                ins.text(0.82, yy, f"CM{k + 1}", ha="center", va="center", fontsize=5.6, color=col, fontweight="bold")
            arrow(ins, 0.64, 0.5, 0.70, 0.5, lw=0.7)
    for x in (0.325, 0.655):
        arrow(ax, x, 0.47, x + 0.022, 0.47, color=GREY)

    # panel b workflow
    bx = fig.add_axes([0, 0.0, 1, 0.5])
    bx.axis("off")
    panel(bx, "b", x=0.0, y=0.94)
    bx.text(0.03, 0.955, "Data and analysis chain shown on this homepage (version 0906)", fontsize=7.5, fontweight="bold",
            transform=bx.transAxes)
    steps = [
        ("Public scRNA-seq\n54 GSE · 43 cancers\n1,322 samples", "#eef3f8"),
        ("Per-sample QC\nraw counts\n4,676,787 cells", "#eef3f8"),
        ("HVG4000 · PCA\nHarmony / BBKNN / scVI\nscIB comparison", "#eef3f8"),
        ("Identity V5→V7\nmarker + source labels\n+ CellTypist evidence", "#f6efe7"),
        ("Per-sample CNV\ninferCNVpy + CopyKAT\nmalignancy tiers", "#f6efe7"),
        ("CoVarNet nsNMF\nK = 9 cellular modules\n(sample × subtype freq.)", "#eaf4ef"),
    ]
    xs = np.linspace(0.025, 0.835, len(steps))
    for (t, fc), x in zip(steps, xs):
        box(bx, x, 0.58, 0.14, 0.26, t, fc=fc, ec="#b7c0cb", fs=6.0)
    for a, b in zip(xs[:-1], xs[1:]):
        arrow(bx, a + 0.142, 0.71, b - 0.002, 0.71, color=GREY)
    steps2 = [
        ("LIANA per sample\nwithin-CM sender→receiver\n1,301 / 1,322 samples", "#eaf4ef"),
        ("Focus axes\nCM01 · CM07 · CM08\nsample-level replication", "#eaf4ef"),
        ("Spatial readouts\nVisium 61 sections\nccRCC CosMx 199,112 cells", "#fbefe9"),
        ("Bulk cohorts\nTCGA 9,034 patients\nICI: Riaz, Rose, IMvigor210…", "#fbefe9"),
    ]
    xs2 = np.linspace(0.12, 0.74, len(steps2))
    for (t, fc), x in zip(steps2, xs2):
        box(bx, x, 0.14, 0.16, 0.26, t, fc=fc, ec="#b7c0cb", fs=6.0)
    for a, b in zip(xs2[:-1], xs2[1:]):
        arrow(bx, a + 0.162, 0.27, b - 0.002, 0.27, color=GREY)
    bx.plot([xs[-1] + 0.07, xs[-1] + 0.07, xs2[0] + 0.08], [0.58, 0.49, 0.49], color=GREY, lw=0.8,
            transform=bx.transAxes)
    arrow(bx, xs2[0] + 0.08, 0.49, xs2[0] + 0.08, 0.405, color=GREY)
    bx.text(0.5, 0.03,
            "discovery = pooled public scRNA-seq (sample as replicate)   ·   readouts = independent spatial and bulk "
            "cohorts (association, not causation)",
            ha="center", fontsize=5.8, color="#5b6573", transform=bx.transAxes)
    save(fig, "fig1_concept_workflow")


# --------------------------------------------------------------------------
# Fig 2  data composition
# --------------------------------------------------------------------------
def fig2() -> None:
    md = V5F / "metadata_audit_20260829" / "visualization_r_20260829" / "tables"
    cc = pd.read_csv(md / "cancer_composition.tsv", sep="\t").sort_values("n_cells")
    sc = pd.read_csv(md / "sample_class_composition.tsv", sep="\t")
    it = pd.read_csv(md / "immunotherapy_composition.tsv", sep="\t")
    mt = pd.read_csv(md / "metastasis_composition.tsv", sep="\t")
    co = pd.read_csv(md / "cohort_composition.tsv", sep="\t")
    pl = pd.read_csv(md / "assay_composition.tsv", sep="\t")
    export(cc, "fig2_cancer_composition.tsv")
    export(sc, "fig2_sample_class_composition.tsv")

    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(3, 3, width_ratios=[1.25, 1, 1], height_ratios=[1, 1, 1], wspace=0.95, hspace=0.95)

    ax = fig.add_subplot(gs[:, 0])
    y = np.arange(len(cc))
    cmap = plt.get_cmap("Blues")
    gmax = cc.n_gse.max()
    ax.barh(y, cc.n_cells / 1e3, color=[cmap(0.35 + 0.6 * g / gmax) for g in cc.n_gse], height=0.72, lw=0)
    for yi, (n, s) in enumerate(zip(cc.n_cells, cc.n_samples)):
        ax.text(n / 1e3 + 8, yi, f"{s}", va="center", fontsize=5.2, color="#5b6573")
    ax.set_yticks(y, cc.cancer_group, fontsize=5.6)
    ax.set_xlabel("Cells (thousands)")
    ax.set_ylim(-0.7, len(cc) - 0.3)
    ax.set_title("Cells per cancer group (n samples at bar end)", loc="left")
    sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.ListedColormap([cmap(0.35 + 0.6 * g / gmax) for g in range(1, gmax + 1)]),
                               norm=matplotlib.colors.Normalize(0.5, gmax + 0.5))
    cax = ax.inset_axes([0.62, 0.08, 0.3, 0.025])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=[1, 4, 8])
    cb.set_label("GSE per group", fontsize=5.6)
    cb.ax.tick_params(labelsize=5.4, length=1.5)
    cb.outline.set_linewidth(0.4)
    panel(ax, "a", x=-0.55)

    # b cohorts
    ax = fig.add_subplot(gs[0, 1])
    cells = co[co.metric == "Cells"].set_index("cohort_plot")["value"]
    samples = co[co.metric == "Samples"].set_index("cohort_plot")["value"]
    labs = ["Existing 2M cohort", "New 2026 cohort"]
    colr = ["#7c93ad", "#c2532b"]
    left = 0
    for lab, c in zip(labs, colr):
        ax.barh(1, cells[lab] / 1e6, left=left, color=c, height=0.55)
        ax.text(left + cells[lab] / 2e6, 1, f"{cells[lab] / 1e6:.2f}M", ha="center", va="center", color="white", fontsize=6)
        left += cells[lab] / 1e6
    left = 0
    for lab, c in zip(labs, colr):
        frac = samples[lab] / samples.sum() * left_total(cells)
        ax.barh(0, frac, left=left, color=c, height=0.55)
        ax.text(left + frac / 2, 0, f"{int(samples[lab])}", ha="center", va="center", color="white", fontsize=6)
        left += frac
    ax.set_yticks([0, 1], ["Samples\n(scaled)", "Cells"])
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)
    ax.set_title("Two cohorts merged (Aug 2026)", loc="left")
    ax.legend(handles=[patches.Patch(color=c, label=l.replace(" cohort", "")) for l, c in zip(labs, colr)],
              frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2, handlelength=1)
    panel(ax, "b")

    # c sample class
    ax = fig.add_subplot(gs[0, 2])
    order = ["Primary/unspecified tumor", "Metastatic tumor", "Normal/adjacent tissue", "Lymph node", "Blood", "Bone marrow", "Unknown"]
    sc = sc.set_index("sample_class_plot").loc[order]
    colc = ["#c2532b", "#e08a5c", "#5fa58a", "#7b5ea7", "#d16d8e", "#b39ddb", "#c9ced6"]
    ax.barh(range(len(order))[::-1], sc.n_samples, color=colc, height=0.7)
    for i, v in zip(range(len(order))[::-1], sc.n_samples):
        ax.text(v + 8, i, str(v), va="center", fontsize=5.4)
    ax.set_yticks(range(len(order))[::-1], [o.replace("/unspecified", "/unspec.") for o in order], fontsize=5.6)
    ax.set_xlabel("Samples")
    ax.set_title("Sample material", loc="left")
    panel(ax, "c", x=-0.75)

    # d immunotherapy + metastasis known
    ax = fig.add_subplot(gs[1, 1:])
    rows = [
        ("Immunotherapy\nannotated", it.set_index("immunotherapy_plot")["n_samples"].reindex(["Yes", "No", "Unknown"]).fillna(0)),
        ("Metastasis\nannotated", mt.set_index("metastasis_plot")["n_samples"].reindex(
            ["Yes", "No", "Primary sample only", "Normal sample (N/A)", "Unknown", "Mixed/not mappable"]).fillna(0)),
    ]
    palette = {"Yes": "#c2532b", "No": "#3b6ea8", "Unknown": "#d5d9df", "Primary sample only": "#7fa7cf",
               "Normal sample (N/A)": "#5fa58a", "Mixed/not mappable": "#a7afba"}
    for i, (lab, s) in enumerate(rows):
        left = 0
        tot = s.sum()
        for k, v in s.items():
            ax.barh(1 - i, v / tot * 100, left=left, color=palette[k], height=0.55, lw=0)
            if v / tot > 0.045:
                ax.text(left + v / tot * 50, 1 - i, f"{int(v)}", ha="center", va="center", fontsize=5.4,
                        color="white" if k not in ("Unknown",) else INK)
            left += v / tot * 100
    ax.set_yticks([1, 0], [r[0] for r in rows])
    ax.legend(handles=[patches.Patch(color=c, label=k) for k, c in palette.items()], frameon=False, fontsize=5.2,
              ncol=6, loc="upper left", bbox_to_anchor=(-0.02, -0.32), handlelength=0.8, columnspacing=0.8)
    ax.set_xlabel("% of 1,322 samples")
    ax.set_xlim(0, 100)
    ax.set_title("Clinical context is sparse in public metadata", loc="left")
    panel(ax, "d", x=-0.17)

    # e assay
    ax = fig.add_subplot(gs[2, 1:])
    pl = pl.sort_values("n_cells")
    lab = [textwrap.shorten(x.replace("single-cell RNA-seq", "scRNA-seq"), 46) for x in pl.library_prep_or_assay]
    ax.barh(range(len(pl)), pl.n_cells / 1e6, color="#7c93ad", height=0.6)
    for i, (n, s, g) in enumerate(zip(pl.n_cells, pl.n_samples, pl.n_gse)):
        ax.text(n / 1e6 + 0.03, i, f"{s} samples · {g} GSE", va="center", fontsize=5.4)
    ax.set_yticks(range(len(pl)), lab, fontsize=5.6)
    ax.set_xlabel("Cells (millions)")
    ax.set_xlim(0, 3.9)
    ax.set_title("Library preparation recorded in GEO text", loc="left")
    panel(ax, "e", x=-0.17)
    save(fig, "fig2_data_composition")


def left_total(cells: pd.Series) -> float:
    return float(cells.sum() / 1e6)


# --------------------------------------------------------------------------
# Fig 3  UMAP atlas (cell coordinates read only, never exported)
# --------------------------------------------------------------------------
def fig3() -> None:
    u = pd.read_csv(V5F / "tables" / "final_annotation_umap_review_200k.tsv.gz", sep="\t",
                    usecols=["major_celltype_final", "malignancy_tier_v3", "annotation_confidence", "umap1", "umap2"])
    u = u.sample(frac=1.0, random_state=7)
    fig = plt.figure(figsize=(7.2, 5.9))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], wspace=0.35, hspace=0.32)

    ax = fig.add_subplot(gs[0, :2])
    for lab, g in u.groupby("major_celltype_final"):
        ax.scatter(g.umap1, g.umap2, s=0.05, color=LINEAGE_COL.get(lab, "#cccccc"), lw=0, rasterized=True)
    cent = u.groupby("major_celltype_final")[["umap1", "umap2"]].median()
    nudge = {"Dendritic": (0.0, -2.2), "Monocyte": (-1.0, 1.2), "Endothelial": (-2.5, -0.6), "Fibroblast": (1.5, 0.4),
             "HSPC": (-1.2, -0.8), "Erythroid": (1.6, 0.9)}
    for lab, r in cent.iterrows():
        if lab in ("Unknown", "Ambiguous", "Platelet"):
            continue
        dx, dy = nudge.get(lab, (0, 0))
        ax.text(r.umap1 + dx, r.umap2 + dy, lab, fontsize=5.8, ha="center", va="center", color=INK,
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
    ax.set_title("Major lineage (V5 identity used by CoVarNet K=9); 200,000 of 4,676,787 cells", loc="left")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel("UMAP1 (BBKNN graph)"); ax.set_ylabel("UMAP2")
    handles = [patches.Patch(color=LINEAGE_COL[k], label=k) for k in ["Unknown", "Ambiguous"]]
    ax.legend(handles=handles, frameon=False, loc="lower left", fontsize=5.6, handlelength=0.9)
    panel(ax, "a", x=-0.04)

    ax = fig.add_subplot(gs[0, 2])
    vc = u.major_celltype_final.value_counts()
    full = pd.read_csv(REPO / "releases" / "annotation_v7_20260905" / "tables" / "major_celltype_v7_counts.tsv", sep="\t")
    strict_major = 4057641
    total = 4676787
    ax.barh([2, 1, 0], [strict_major / total * 100, 13.24, 2.78], color=["#2c6e8f", "#9cc3d5", "#a7afba"], height=0.6)
    for yy, v, t in zip([2, 1, 0], [86.76, 13.24, 2.78], ["strict major", "boundary / provisional", "exact Unknown+Ambiguous"]):
        ax.text(v + 1.5, yy, f"{v:.1f}%", va="center", fontsize=6)
    ax.set_yticks([2, 1, 0], ["Strict major\n(V7)", "Non-strict major\n(V7)", "Unknown +\nAmbiguous (V7)"], fontsize=5.8)
    ax.set_xlim(0, 112)
    ax.set_xlabel("% of 4,676,787 cells")
    ax.set_title("V7 (frozen 2026-09-05) uncertainty", loc="left")
    v5_unres = (vc.get("Unknown", 0) + vc.get("Ambiguous", 0)) / vc.sum() * 100
    ax.text(0.98, 0.08, f"V5 UMAP sample:\nUnknown+Ambiguous = {v5_unres:.1f}%", transform=ax.transAxes, ha="right",
            fontsize=5.6, color="#5b6573")
    panel(ax, "b", x=-0.42)
    export(full, "fig3_v7_major_counts.tsv")

    ax = fig.add_subplot(gs[1, 0])
    tier_col = {"high_confidence_malignant": "#8c1d18", "cluster_supported_malignant": "#d9603b",
                "sensitive_union_only": "#f0a35e", "uncertain": "#e9d8a6", "non_malignant": "#5fa58a",
                "not_evaluable": "#b0b7c1", "not_candidate": "#e4e7eb"}
    for lab in ["not_candidate", "not_evaluable", "non_malignant", "uncertain", "sensitive_union_only",
                "cluster_supported_malignant", "high_confidence_malignant"]:
        g = u[u.malignancy_tier_v3 == lab]
        ax.scatter(g.umap1, g.umap2, s=0.05, color=tier_col[lab], lw=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("Per-sample CNV malignancy tier", loc="left")
    ax.legend(handles=[patches.Patch(color=c, label=k.replace("_", " ")) for k, c in tier_col.items()], frameon=False,
              fontsize=5.0, loc="upper left", bbox_to_anchor=(-0.02, -0.0), ncol=2, handlelength=0.8, columnspacing=0.6)
    panel(ax, "c", x=-0.08)

    ax = fig.add_subplot(gs[1, 1])
    conf_col = {"high": "#14324a", "medium": "#2c6e8f", "provisional": "#9cc3d5", "coarse": "#e0b04a",
                "unresolved": "#e4e7eb"}
    for lab in ["unresolved", "coarse", "provisional", "medium", "high"]:
        g = u[u.annotation_confidence == lab]
        ax.scatter(g.umap1, g.umap2, s=0.05, color=conf_col[lab], lw=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("Fine-identity confidence", loc="left")
    ax.legend(handles=[patches.Patch(color=c, label=k) for k, c in conf_col.items()], frameon=False, fontsize=5.0,
              loc="upper left", bbox_to_anchor=(-0.02, 0.0), ncol=3, handlelength=0.8, columnspacing=0.6)
    panel(ax, "d", x=-0.08)

    ax = fig.add_subplot(gs[1, 2])
    tr = pd.read_csv(V7 / "major_v5_to_v7_transition.tsv", sep="\t").set_index("major_celltype_final")
    rows = {}
    for src in ["Unknown", "Ambiguous"]:
        r = tr.loc[src]
        tot = r.sum()
        rows[src] = {
            "Resolved lineage": r[[c for c in r.index if c in LINEAGE_COL and c not in ("Unknown", "Ambiguous")]].sum() / tot,
            "Boundary class": r[[c for c in r.index if "boundary" in c]].sum() / tot,
            "Still Unknown/Ambig.": r[["Unknown", "Ambiguous"]].sum() / tot,
        }
    rr = pd.DataFrame(rows).T
    export(rr.reset_index().rename(columns={"index": "v5_major"}), "fig3_v5_to_v7_unknown_resolution.tsv")
    cols3 = ["#2c6e8f", "#9cc3d5", "#a7afba"]
    left = np.zeros(len(rr))
    for c, col in zip(rr.columns, cols3):
        ax.barh(range(len(rr)), rr[c] * 100, left=left, color=col, height=0.55, label=c)
        for i, (l0, v) in enumerate(zip(left, rr[c] * 100)):
            if v > 8:
                ax.text(l0 + v / 2, i, f"{v:.0f}%", ha="center", va="center", fontsize=5.6, color="white" if col != "#a7afba" else INK)
        left += rr[c].values * 100
    v5n = tr.loc[["Unknown", "Ambiguous"]].sum(axis=1)
    ax.set_yticks(range(len(rr)), [f"V5 {k}\n(n={v5n[k] / 1e3:.0f}k)" for k in rr.index], fontsize=5.8)
    ax.set_xlabel("% re-assigned in V7")
    ax.set_xlim(0, 100)
    ax.legend(frameon=False, fontsize=5.2, loc="upper left", bbox_to_anchor=(-0.05, -0.2), ncol=2, handlelength=0.8)
    ax.set_title("V5 → V7 fate of uncertain cells", loc="left")
    panel(ax, "e", x=-0.42)
    save(fig, "fig3_atlas_umap_annotation")


# --------------------------------------------------------------------------
# Fig 4  CM landscape (K = 9)
# --------------------------------------------------------------------------
def fig4() -> None:
    w = pd.read_csv(V5F / "tables" / "cm_membership_weights.tsv", sep="\t")
    act = pd.read_csv(V5F / "tables" / "cm_cancer_activity.tsv", sep="\t")
    inf = pd.read_csv(CMP / "cm_immune_infiltration_associations.tsv", sep="\t")
    cms = [f"CM{i:02d}" for i in range(1, 10)]

    top = w[w["rank"] <= 6]
    members = []
    for cm in cms:
        for s in top[top.cm == cm].sort_values("rank").subCluster:
            if s not in members:
                members.append(s)
    mat = w.pivot_table(index="subCluster", columns="cm", values="weight_fraction").reindex(index=members, columns=cms).fillna(0)
    export(w[w.top10], "fig4_cm_top10_membership.tsv")

    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.05, 1], height_ratios=[1.35, 1], wspace=0.75, hspace=0.5)

    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(mat.values, cmap=SEQ, aspect="auto", vmin=0, vmax=0.2)
    ax.set_xticks(range(9), cms, rotation=90)
    labs = []
    for m in members:
        unres = "unresolved" in m or m in ("Unknown", "Ambiguous") or "uncertain" in m or "not_evaluable" in m
        labs.append(("◦ " if unres else "") + short_state(m))
    ax.set_yticks(range(len(members)), labs, fontsize=5.0)
    for i, m in enumerate(members):
        ax.add_patch(patches.Rectangle((-0.95, i - 0.5), 0.4, 1, color=LINEAGE_COL[lineage_of(m)], clip_on=False,
                                       transform=ax.transData))
    ax.tick_params(axis="y", pad=7)
    for j, cm in enumerate(cms):
        if cm in CM_COL:
            ax.get_xticklabels()[j].set_color(CM_COL[cm]); ax.get_xticklabels()[j].set_fontweight("bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, shrink=0.4)
    cb.set_label("Member weight fraction (W)", fontsize=5.8); cb.ax.tick_params(labelsize=5.4)
    ax.set_title("Top-6 members per module  (◦ = unresolved / uncertain node)", loc="left", fontsize=6.6)
    panel(ax, "a", x=-0.75, y=1.02)

    # b lineage composition of top10 weight
    ax = fig.add_subplot(gs[0, 1])
    t10 = w[w.top10].copy()
    t10["lin"] = t10.subCluster.map(lineage_of)
    comp = t10.groupby(["cm", "lin"]).weight.sum().unstack(fill_value=0)
    comp = comp.div(comp.sum(axis=1), axis=0).reindex(cms)
    order = [k for k in LINEAGE_COL if k in comp.columns]
    left = np.zeros(9)
    for lin in order:
        ax.barh(range(9)[::-1], comp[lin] * 100, left=left, color=LINEAGE_COL[lin], height=0.72, lw=0, label=lin)
        left += comp[lin].values * 100
    ax.set_yticks(range(9)[::-1], cms)
    for t in ax.get_yticklabels():
        if t.get_text() in CM_COL:
            t.set_color(CM_COL[t.get_text()]); t.set_fontweight("bold")
    unres = t10.assign(u=t10.review_flag).groupby("cm").apply(lambda g: (g.weight * g.u).sum() / g.weight.sum()).reindex(cms)
    for i, cm in enumerate(cms):
        ax.text(102, 8 - i, f"{unres[cm] * 100:.0f}%", va="center", fontsize=5.6, color="#5b6573")
    ax.text(102, 9.0, "review-\nflagged", fontsize=5.2, color="#5b6573", va="bottom")
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of top-10 member weight by lineage")
    ax.legend(frameon=False, fontsize=4.8, ncol=5, loc="upper left", bbox_to_anchor=(-0.12, -0.11), handlelength=0.7,
              columnspacing=0.4, handletextpad=0.3)
    ax.set_title("Lineage make-up of each module", loc="left")
    panel(ax, "b", x=-0.2, y=1.02)
    export(comp.reset_index(), "fig4_cm_lineage_composition.tsv")

    # c cancer activity
    ax = fig.add_subplot(gs[1, 0])
    a2 = act[act.n_samples >= 10]
    piv = a2.pivot_table(index="cancer_type", columns="cm", values="median_activity").reindex(columns=cms)
    z = (piv - piv.mean()) / piv.std()
    order_c = z[["CM01", "CM07", "CM08"]].mean(axis=1).sort_values().index
    z = z.loc[order_c]
    im = ax.imshow(z.values, cmap=DIV, norm=TwoSlopeNorm(0, -2.5, 2.5), aspect="auto")
    ax.set_xticks(range(9), cms, rotation=90)
    for j, cm in enumerate(cms):
        if cm in CM_COL:
            ax.get_xticklabels()[j].set_color(CM_COL[cm]); ax.get_xticklabels()[j].set_fontweight("bold")
    ax.set_yticks(range(len(z)), [textwrap.shorten(c, 34) for c in z.index], fontsize=4.9)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, shrink=0.45)
    cb.set_label("Median activity\n(z across cancers)", fontsize=5.6); cb.ax.tick_params(labelsize=5.4)
    ax.set_title(f"Module activity by cancer label ({len(z)} labels with ≥10 samples)", loc="left", fontsize=6.6)
    panel(ax, "c", x=-0.75, y=1.02)
    export(a2, "fig4_cm_cancer_activity_ge10.tsv")

    # d immune infiltration within cancer
    ax = fig.add_subplot(gs[1, 1])
    feats = ["spp1_macrophage_fraction", "inflammatory_myeloid_fraction", "myeloid_fraction", "caf_fraction",
             "stromal_activation_proxy", "cytotoxic_fraction", "lymphoid_fraction", "immune_hot_proxy",
             "malignant_strict_fraction"]
    flab = ["SPP1 macrophage", "Inflammatory myeloid", "Myeloid", "CAF", "Stromal activation", "Cytotoxic",
            "Lymphoid", "Immune-hot proxy", "Strict malignant"]
    p = inf.pivot_table(index="cm", columns="feature", values="within_cancer_spearman").reindex(index=cms, columns=feats)
    q = inf.pivot_table(index="cm", columns="feature", values="within_cancer_fdr").reindex(index=cms, columns=feats)
    im = ax.imshow(p.values.T, cmap=DIV, norm=TwoSlopeNorm(0, -0.6, 0.6), aspect="auto")
    for i in range(len(feats)):
        for j in range(9):
            if q.values[j, i] < 0.05 and abs(p.values[j, i]) >= 0.2:
                ax.text(j, i, "•", ha="center", va="center", fontsize=6, color="white" if abs(p.values[j, i]) > 0.4 else INK)
    ax.set_xticks(range(9), cms, rotation=90)
    for j, cm in enumerate(cms):
        if cm in CM_COL:
            ax.get_xticklabels()[j].set_color(CM_COL[cm]); ax.get_xticklabels()[j].set_fontweight("bold")
    ax.set_yticks(range(len(feats)), flab, fontsize=5.6)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03, shrink=0.6)
    cb.set_label("Within-cancer Spearman ρ\n(n = 1,303 samples)", fontsize=5.6); cb.ax.tick_params(labelsize=5.4)
    ax.set_title("Activity vs. infiltration (• FDR<0.05, |ρ|≥0.2)", loc="left", fontsize=6.4)
    panel(ax, "d", x=-0.45, y=1.02)
    export(inf, "fig4_cm_infiltration_associations.tsv")
    save(fig, "fig4_cm_landscape")


# --------------------------------------------------------------------------
# Fig 5  communication focus axes
# --------------------------------------------------------------------------
def fig5() -> None:
    rep = pd.read_csv(CHAIN / "focus_axes_replication_summary.tsv", sep="\t")
    byc = pd.read_csv(ECH / "focus_axes_by_cancer.tsv", sep="\t")
    rep["label"] = rep.apply(lambda r: f"{r.ligand}→{r.receptor}  {short_state(r.source)}→{short_state(r.target)}", axis=1)
    order_cm = {"CM01": 0, "CM07": 1, "CM08": 2}
    rep = rep.sort_values(["cm", "support_frac_standard_0.05"], key=lambda s: s.map(order_cm) if s.name == "cm" else s)
    rep = rep.reset_index(drop=True)
    export(rep, "fig5_focus_axes_replication.tsv")

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(3, 3, height_ratios=[0.8, 1.25, 1.05], hspace=0.62, wspace=0.55)

    # a: sender -> receiver ladders
    for k, cm in enumerate(["CM01", "CM07", "CM08"]):
        ax = fig.add_subplot(gs[0, k])
        ax.set_xlim(0, 1); ax.set_ylim(-0.6, 5.6); ax.axis("off")
        sub = rep[rep.cm == cm].sort_values("support_frac_standard_0.05", ascending=False).reset_index(drop=True)
        col = CM_COL[cm]
        for i, r in sub.iterrows():
            yy = 5 - i * (5.0 / max(1, len(sub) - 1)) if len(sub) > 1 else 2.5
            for xx, lab, ha in ((0.0, r.source, "left"), (1.0, r.target, "right")):
                ax.add_patch(patches.Circle((0.035 if ha == "left" else 0.965, yy - 0.42), 0.022,
                                            color=LINEAGE_COL[lineage_of(lab)], transform=ax.transData, clip_on=False))
                ax.text(0.07 if ha == "left" else 0.93, yy - 0.42, short_state(lab), ha=ha, va="center", fontsize=5.3,
                        color=INK)
            ax.annotate("", xy=(0.62, yy - 0.42), xytext=(0.38, yy - 0.42),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=0.5 + 2.2 * r.robust_strict_vs_std,
                                        mutation_scale=6))
            ax.text(0.5, yy + 0.05, f"{r.ligand}→{r.receptor}", ha="center", va="center", fontsize=5.4, color=col,
                    fontweight="bold" if r.hypothesis == "primary" else "normal")
        title = {"CM01": "CM01 · perivascular SPP1-TAM axis",
                 "CM07": "CM07 · myCAF–epithelial desmoplasia",
                 "CM08": "CM08 · inflammatory myeloid / myCAF–CD8"}[cm]
        ax.set_title(title, fontsize=6.4, color=col, fontweight="bold", loc="left")
        if k == 0:
            panel(ax, "a", x=-0.1, y=1.02)
    fig.text(0.125, 0.69, "sender (left) → receiver (right); arrow width proportional to strict/standard retention; "
             "bold = primary hypothesis axis", fontsize=5.6, color="#5b6573")

    # b: support at three thresholds
    ax = fig.add_subplot(gs[1, :2])
    y = np.arange(len(rep))[::-1]
    for th, mk, col in [("strict_0.01", "o", "#14324a"), ("standard_0.05", "o", "#2c6e8f"), ("loose_0.10", "o", "#9cc3d5")]:
        ax.scatter(rep[f"support_frac_{th}"], y, s=14, marker=mk, color=col, zorder=3, lw=0, label=f"rank ≤ {th.split('_')[1]}")
    for yi, (lo, hi) in zip(y, zip(rep["support_frac_strict_0.01"], rep["support_frac_loose_0.10"])):
        ax.plot([lo, hi], [yi, yi], color=GREY, lw=0.6, zorder=1)
    ax.set_yticks(y, rep.label, fontsize=5.4)
    for t, cm in zip(ax.get_yticklabels(), rep.cm):
        t.set_color(CM_COL[cm])
    ax.set_xlim(0, 1.05)
    ax.axvline(0.5, color=LIGHT, lw=0.6, zorder=0)
    ax.set_xlabel("Support fraction = samples passing / samples where both cell states are present (≥20 cells)")
    ax.legend(frameon=False, loc="lower left", fontsize=5.4, handletextpad=0.2)
    ax.set_title("Sample-level replication of 12 focus axes at three LIANA rank thresholds", loc="left")
    panel(ax, "b", x=-0.62, y=1.02)

    ax = fig.add_subplot(gs[1, 2])
    ax.scatter(rep["n_gse_standard_0.05"], rep["n_cancer_standard_0.05"], s=rep["n_samples_standard_0.05"] / 4,
               color=[CM_COL[c] for c in rep.cm], alpha=0.65, lw=0.4, ec="white")
    offs = {"S100A9": (0.6, 0.9), "COL1A1": (1.0, -0.2), "LGALS1": (-9.5, 0.9)}
    for r in rep.itertuples():
        if r.hypothesis == "primary":
            dx, dy = offs.get(r.ligand, (0.6, 0))
            ax.text(rep.loc[r.Index, "n_gse_standard_0.05"] + dx, rep.loc[r.Index, "n_cancer_standard_0.05"] + dy,
                    f"{r.ligand}–{r.receptor}", fontsize=5.3, color=CM_COL[r.cm])
    ax.set_ylim(20.5, 35.5)
    ax.set_xlabel("GSE supporting (rank ≤ 0.05)")
    ax.set_ylabel("Cancer labels supporting")
    for s in (100, 300):
        ax.scatter([], [], s=s / 4, color=GREY, alpha=0.5, label=f"{s} samples")
    ax.legend(frameon=False, fontsize=5.2, loc="upper left", labelspacing=1.1)
    ax.set_title("Breadth of support", loc="left")
    panel(ax, "c", x=-0.3, y=1.02)

    # d per-cancer heatmap
    ax = fig.add_subplot(gs[2, :])
    canc = byc.groupby("cancer_type").n_samples.sum().sort_values(ascending=False)
    canc = canc[canc >= 15].index[:30]
    hm = byc.pivot_table(index="axis", columns="cancer_type", values="n_samples", aggfunc="sum").reindex(columns=canc)
    hm = hm.reindex(rep.axis)
    im = ax.imshow(np.log10(hm.fillna(0).values + 1), cmap=SEQ, aspect="auto", vmin=0, vmax=np.log10(60))
    ax.set_xticks(range(len(canc)), [textwrap.shorten(c, 26) for c in canc], rotation=60, ha="right", fontsize=5.0)
    ax.set_yticks(range(len(hm)), [f"{a.split('_')[0]} {r.ligand}→{r.receptor}" for a, r in zip(hm.index, rep.itertuples())],
                  fontsize=5.2)
    for t, cm in zip(ax.get_yticklabels(), rep.cm):
        t.set_color(CM_COL[cm])
    cb = fig.colorbar(im, ax=ax, fraction=0.015, pad=0.01)
    cb.set_ticks([0, 1, np.log10(60)], labels=["0", "9", "59"])
    cb.set_label("Supporting samples", fontsize=5.6); cb.ax.tick_params(labelsize=5.2)
    ax.set_title("Where each axis is observed: supporting samples per cancer label (labels with ≥15 supporting rows)",
                 loc="left")
    panel(ax, "d", x=-0.11, y=1.02)
    export(byc, "fig5_focus_axes_by_cancer.tsv")
    save(fig, "fig5_communication_axes")


# --------------------------------------------------------------------------
# Fig 6  spatial readouts
# --------------------------------------------------------------------------
def fig6() -> None:
    cs = pd.read_csv(SPAT / "communication_cross_slice_summary.tsv", sep="\t")
    mi = pd.read_csv(SPAT / "cm_morans_i.tsv", sep="\t")
    st = pd.read_csv(SPAT / "slice_analysis_status.tsv", sep="\t")
    t3 = pd.read_csv(ECH / "cosmx_col_ddr1_T3_cd8_exclusion.tsv", sep="\t")
    t3b = pd.read_csv(ECH / "cosmx_col_ddr1_T3b_cd8_positioning.tsv", sep="\t")
    t4 = pd.read_csv(ECH / "cosmx_col_ddr1_T4_receiver_emt.tsv", sep="\t")
    export(cs, "fig6_visium_cross_slice_summary.tsv")
    export(mi.drop(columns=[c for c in mi.columns if c in ("key",)]), "fig6_visium_cm_morans_i.tsv")
    export(t3, "fig6_cosmx_T3.tsv"); export(t3b, "fig6_cosmx_T3b.tsv"); export(t4, "fig6_cosmx_T4.tsv")
    cancer = st.drop_duplicates("dataset_id").set_index("dataset_id").cancer_type.to_dict()
    nsl = st[st.status == "PASS"].groupby("dataset_id").size().to_dict()

    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1], width_ratios=[1.5, 1, 1], hspace=0.75, wspace=0.75)

    ax = fig.add_subplot(gs[0, 0])
    piv = cs.pivot_table(index="axis_id", columns="dataset_id", values="median_permutation_z")
    piv = piv.loc[piv.median(axis=1).sort_values(ascending=False).index]
    im = ax.imshow(piv.values, cmap=DIV, norm=TwoSlopeNorm(0, -10, 10), aspect="auto")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=5.0, color="white" if abs(v) > 6 else INK)
    ax.set_xticks(range(piv.shape[1]), [f"{d} · {cancer.get(d, '')} · {nsl.get(d, 0)}" for d in piv.columns],
                  fontsize=5.0, rotation=40, ha="right")
    ylab = []
    for a in piv.index:
        p = a.split("_")
        ylab.append(f"{p[0]} {p[1]}→{p[2]}  {' '.join(p[3:]).replace(' to ', '→')}")
    ax.set_yticks(range(piv.shape[0]), ylab, fontsize=5.2)
    for t, a in zip(ax.get_yticklabels(), piv.index):
        t.set_color(CM_COL.get(a[:4], INK))
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Median permutation z", fontsize=5.4)
    cb.ax.tick_params(labelsize=5.2)
    ax.set_title("Visium: sender×receiver co-enrichment", loc="left", fontsize=6.6)
    panel(ax, "a", x=-0.7, y=1.02)

    ax = fig.add_subplot(gs[0, 1:])
    mi["cancer"] = mi.dataset_id.map(cancer)
    cms_show = sorted(mi.CM.unique())
    datasets = list(piv.columns)
    width = 0.8 / max(1, len(datasets))
    dcol = plt.get_cmap("tab10")
    for k, d in enumerate(datasets):
        for j, cm in enumerate(cms_show):
            vals = mi[(mi.dataset_id == d) & (mi.CM == cm)].morans_I.values
            if len(vals) == 0:
                continue
            x = j - 0.4 + width * (k + 0.5)
            ax.scatter(np.full(len(vals), x) + np.random.default_rng(k + j).normal(0, width / 6, len(vals)), vals, s=3,
                       color=dcol(k), alpha=0.7, lw=0)
            ax.plot([x - width / 3, x + width / 3], [np.median(vals)] * 2, color=INK, lw=0.8)
    ax.set_xticks(range(len(cms_show)), cms_show)
    for t in ax.get_xticklabels():
        if t.get_text() in CM_COL:
            t.set_color(CM_COL[t.get_text()]); t.set_fontweight("bold")
    ax.set_ylabel("Moran's I of projected module score")
    ax.set_ylim(0, 1)
    ax.legend(handles=[patches.Patch(color=dcol(k), label=f"{d} ({cancer.get(d, '')})") for k, d in enumerate(datasets)],
              frameon=False, fontsize=5.2, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.36))
    ax.set_title("Moran's I of projected module scores (each dot = section; all p ≤ 0.001)", loc="left", fontsize=6.6)
    panel(ax, "b", x=-0.12, y=1.02)

    # c CosMx distance
    ax = fig.add_subplot(gs[1, 0])
    t3b = t3b[t3b.tissue != "RCC4"]
    x = np.arange(len(t3b))
    ax.bar(x - 0.18, t3b.med_dist_to_interfaceHI_tumour, width=0.34, color="#c2532b", label="collagen×DDR1-high tumour")
    ax.bar(x + 0.18, t3b.med_dist_to_interfaceLO_tumour, width=0.34, color="#7c93ad", label="collagen×DDR1-low tumour")
    for xi, r in zip(x, t3b.itertuples()):
        ax.text(xi, max(r.med_dist_to_interfaceHI_tumour, r.med_dist_to_interfaceLO_tumour) + 18,
                f"n CD8 = {r.n_cd8:,}\np = {r.mwu_HI_gt_LO_p:.0e}", ha="center", fontsize=5.2)
    ax.set_xticks(x, t3b.tissue)
    ax.set_ylabel("Median CD8 distance (µm)")
    ax.set_ylim(0, 820)
    ax.legend(frameon=False, fontsize=5.3, loc="upper left")
    ax.set_title("ccRCC CosMx: CD8 distance to tumour cells", loc="left", fontsize=6.4)
    panel(ax, "c", x=-0.2, y=1.02)

    ax = fig.add_subplot(gs[1, 1])
    xs = np.arange(len(t3))
    ax.bar(xs - 0.18, t3.beta_col_exp, 0.34, color="#2a7f62", label="local collagen exposure")
    ax.bar(xs + 0.18, t3.beta_ddr1, 0.34, color="#c2532b", label="tumour DDR1")
    ax.axhline(0, color=INK, lw=0.5)
    ax.set_xticks(xs, t3.tissue)
    ax.set_ylabel("β (CD8 neighbours per tumour cell)")
    ax.set_ylim(-0.24, 0)
    ax.legend(frameon=False, fontsize=5.3, loc="lower left")
    ax.set_title("T3: CD8 neighbours", loc="left")
    panel(ax, "d", x=-0.32, y=1.02)

    ax = fig.add_subplot(gs[1, 2])
    bcols = [c for c in t4.columns if c.startswith("beta")]
    xs = np.arange(len(t4))
    for k, c in enumerate(bcols[:2]):
        ax.bar(xs + (k - 0.5) * 0.34, t4[c], 0.34, color=["#2a7f62", "#c2532b"][k],
               label={0: "local collagen", 1: "tumour DDR1"}[k])
    ax.axhline(0, color=INK, lw=0.5)
    ax.set_xticks(xs, t4.tissue)
    ax.set_ylabel("β (tumour EMT/MMP score)")
    ax.set_ylim(0, 0.068)
    ax.legend(frameon=False, fontsize=5.3, loc="upper right")
    ax.set_title("T4: receiver EMT/MMP", loc="left")
    panel(ax, "e", x=-0.32, y=1.02)
    save(fig, "fig6_spatial_readouts")


# --------------------------------------------------------------------------
# Fig 7  bulk and ICI cohorts
# --------------------------------------------------------------------------
def forest(ax, labels, hr, lo, hi, colors, xlab, xlim=None):
    y = np.arange(len(labels))[::-1]
    for yi, h, l, u, c in zip(y, hr, lo, hi, colors):
        ax.plot([l, u], [yi, yi], color=c, lw=0.9)
        ax.scatter(h, yi, s=13, color=c, zorder=3, lw=0)
    ax.axvline(1, color=GREY, lw=0.6, ls="--")
    ax.set_yticks(y, labels)
    ax.set_xscale("log")
    ax.set_xlabel(xlab)
    if xlim:
        ax.set_xlim(*xlim)
    lo_, hi_ = ax.get_xlim()
    cand = [0.25, 0.5, 1, 2, 4, 8] if hi_ / lo_ > 4 else [0.8, 0.9, 1, 1.1, 1.2, 1.3]
    ax.set_xticks([c for c in cand if lo_ <= c <= hi_])
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def fig7() -> None:
    strat = pd.read_csv(BULK / "tcga_cox_v1" / "pancancer_stratified_cox.tsv", sep="\t")
    per = pd.read_csv(BULK / "tcga_cox_v1" / "per_cancer_cox.tsv", sep="\t")
    resp = pd.read_csv(BULK / "ici_bulk_cox_v1" / "response_tests.tsv", sep="\t")
    surv = pd.read_csv(BULK / "ici_bulk_cox_v1" / "survival_tests.tsv", sep="\t")
    prox = pd.read_csv(ECH / "bulk_proxy_TCGA_models.tsv", sep="\t")
    proxi = pd.read_csv(ECH / "bulk_proxy_ICI_survival_models.tsv", sep="\t")
    for df, n in [(strat, "fig7_tcga_stratified_cox.tsv"), (per, "fig7_tcga_per_cancer_cox.tsv"),
                  (resp, "fig7_ici_response_tests.tsv"), (surv, "fig7_ici_survival_tests.tsv"),
                  (prox, "fig7_tcga_lr_proxy_models.tsv"), (proxi, "fig7_ici_lr_proxy_survival.tsv")]:
        export(df, n)

    fig = plt.figure(figsize=(7.2, 6.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], width_ratios=[1, 1.35, 1.1], hspace=0.55, wspace=1.15)

    ax = fig.add_subplot(gs[0, 0])
    s = strat[strat.endpoint == "OS"].copy()
    s = s.sort_values("HR")
    cols = [CM_COL.get(x, "#5b6573" if x.startswith("CM") else "#7b5ea7") for x in s.signature]
    forest(ax, [f"{x}" for x in s.signature], s.HR, s.CI_low, s.CI_high, cols[::-1][::-1],
           "Hazard ratio per SD (TCGA OS)", (0.85, 1.35))
    for t, x in zip(ax.get_yticklabels(), s.signature[::-1]):
        pass
    for t in ax.get_yticklabels():
        if t.get_text() in CM_COL:
            t.set_color(CM_COL[t.get_text()]); t.set_fontweight("bold")
    ax.set_title("TCGA stratified Cox\nn=9,034 · 2,664 deaths", loc="left", fontsize=6.4)
    panel(ax, "a", x=-0.55, y=1.02)

    ax = fig.add_subplot(gs[0, 1:])
    sigs = ["CM01", "CM07", "CM08", "desmoplasia_myCAF", "TLS_niche3"]
    p = per[(per.endpoint == "OS") & per.signature.isin(sigs)]
    hrm = p.pivot_table(index="signature", columns="cancer", values="logHR").reindex(sigs)
    qm = p.pivot_table(index="signature", columns="cancer", values="BH_q").reindex(sigs)
    order = hrm.loc[["CM01", "CM07", "CM08"]].mean().sort_values().index
    hrm, qm = hrm[order], qm[order]
    im = ax.imshow(hrm.values, cmap=DIV, norm=TwoSlopeNorm(0, -0.7, 0.7), aspect="auto")
    for i in range(hrm.shape[0]):
        for j in range(hrm.shape[1]):
            if qm.values[i, j] < 0.05:
                ax.text(j, i, "•", ha="center", va="center", fontsize=6.5, color="white" if abs(hrm.values[i, j]) > 0.35 else INK)
    ax.set_xticks(range(hrm.shape[1]), hrm.columns, rotation=90, fontsize=5.4)
    ax.set_yticks(range(len(sigs)), ["CM01 signature", "CM07 signature", "CM08 signature", "desmoplasia myCAF (21 genes)",
                                     "TLS niche (Meylan)"], fontsize=5.6)
    for t, sgn in zip(ax.get_yticklabels(), sigs):
        if sgn in CM_COL:
            t.set_color(CM_COL[sgn]); t.set_fontweight("bold")
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("log HR per SD (OS)", fontsize=5.6); cb.ax.tick_params(labelsize=5.2)
    ax.set_title("Cancer-specific direction of the same signatures (• BH q < 0.05 within signature)", loc="left", fontsize=6.4)
    panel(ax, "b", x=-0.22, y=1.02)

    ax = fig.add_subplot(gs[1, 0])
    r = resp[resp.signature.isin(["CM01", "CM02", "CM04", "CM07", "CM08", "TLS_niche3", "desmoplasia_myCAF"])]
    rp = r.pivot_table(index="signature", columns="cohort", values="cliffs_delta")
    rq = r.pivot_table(index="signature", columns="cohort", values="BH_q")
    nmap = r.groupby("cohort").apply(lambda g: f"{int(g.n_R.iloc[0])}R/{int(g.n_NR.iloc[0])}NR").to_dict()
    im = ax.imshow(rp.values, cmap=DIV, norm=TwoSlopeNorm(0, -0.6, 0.6), aspect="auto")
    for i in range(rp.shape[0]):
        for j in range(rp.shape[1]):
            v = rp.values[i, j]
            star = "*" if rq.values[i, j] < 0.1 else ""
            ax.text(j, i, f"{v:+.2f}{star}", ha="center", va="center", fontsize=4.9, color="white" if abs(v) > 0.4 else INK)
    ax.set_xticks(range(rp.shape[1]), [f"{c.split('_')[0]}{' pre' if 'pre' in c else ''}\n{nmap[c]}" for c in rp.columns],
                  fontsize=4.8, rotation=45, ha="right")
    ax.set_yticks(range(rp.shape[0]), rp.index, fontsize=5.6)
    for t in ax.get_yticklabels():
        if t.get_text() in CM_COL:
            t.set_color(CM_COL[t.get_text()]); t.set_fontweight("bold")
    ax.set_title("ICI response: Cliff's δ (R vs NR)\n* BH q < 0.1", loc="left", fontsize=6.4)
    panel(ax, "c", x=-0.6, y=1.02)

    ax = fig.add_subplot(gs[1, 1])
    sv = surv[surv.signature.isin(["CM01", "CM04", "CM07", "CM08", "TLS_niche3", "desmoplasia_myCAF"])].copy()
    sv = sv.sort_values(["cohort", "signature"])
    abbr = {"Hugo2016": "Hugo", "IMvigor210": "IMv210", "Jung2019": "Jung", "Rose2021": "Rose", "Riaz2017": "Riaz"}
    sname = {"desmoplasia_myCAF": "desmo.", "TLS_niche3": "TLS"}
    labels = [f"{abbr.get(c.split('_')[0], c)} {e} · {sname.get(s, s)}" for c, e, s in
              zip(sv.cohort, sv.endpoint, sv.signature)]
    cols = [CM_COL.get(x, "#5b6573") for x in sv.signature]
    forest(ax, labels, sv.HR, sv.CI_low, sv.CI_high, cols, "HR per SD", (0.25, 8))
    ax.tick_params(axis="y", labelsize=4.8)
    ax.set_title("ICI survival per cohort", loc="left", fontsize=6.4)
    panel(ax, "d", x=-0.62, y=1.02)

    ax = fig.add_subplot(gs[1, 2])
    pick = prox[(prox.endpoint == "OS") & prox.proxy.str.startswith(("min_", "mean_"))].copy()
    pick = pick[pick.proxy.isin(["min_SPP1_macrophage", "min_LGALS1_myCAF_CD8", "mean_LGALS1_myCAF_CD8",
                                 "min_collagen_DDR1", "min_COL1A1_only_DDR1"]) | pick.axis.str.contains("SPP1|LGALS1|DDR1")]
    pick = pick.drop_duplicates("proxy")
    pick = pick[pick.proxy.str.startswith("min_") | (pick.proxy == "mean_LGALS1_myCAF_CD8")]
    lab = [f"{p.replace('min_', '').replace('mean_', 'mean ')}\n| {par}" for p, par in zip(pick.proxy, pick.parent)]
    colmap = lambda a: CM_COL["CM01"] if "SPP1" in a else CM_COL["CM08"] if "LGALS1" in a else CM_COL["CM07"]
    y = np.arange(len(pick))[::-1]
    for yi, (hp, hpar, lrt, a) in enumerate(zip(pick.HR_proxy_in_C, pick.HR_parent_in_C, pick.lrt_C_vs_A_p, pick.axis)):
        yy = y[yi]
        ax.scatter(hpar, yy - 0.15, marker="s", s=12, color=GREY, lw=0, zorder=3)
        ax.scatter(hp, yy + 0.15, marker="o", s=14, color=colmap(a), lw=0, zorder=3)
        ax.text(1.42, yy, f"LRT p={lrt:.0e}", fontsize=4.9, va="center", color="#5b6573")
    ax.axvline(1, color=GREY, ls="--", lw=0.6)
    ax.set_yticks(y, lab, fontsize=4.9)
    ax.set_xlim(0.85, 1.4)
    ax.set_xlabel("HR in joint model C (TCGA OS)")
    ax.scatter([], [], marker="o", color=INK, s=12, label="LR proxy (min of L, R)")
    ax.scatter([], [], marker="s", color=GREY, s=12, label="parent signature")
    ax.legend(frameon=False, fontsize=5.0, loc="upper center", bbox_to_anchor=(0.45, -0.18), ncol=1)
    ax.set_title("LR proxy beyond its parent\nsignature (bulk proxy ≠ LIANA)", loc="left", fontsize=6.4)
    panel(ax, "e", x=-0.95, y=1.02)
    save(fig, "fig7_bulk_ici_readouts")


# --------------------------------------------------------------------------
# Sub-page figures
# --------------------------------------------------------------------------
def figS1_datasets() -> None:
    fig = plt.figure(figsize=(7.2, 4.6))
    ax = fig.add_axes([0.02, 0.50, 0.96, 0.46])
    ax.axis("off")
    panel(ax, "a", x=0.0, y=0.96)
    ax.text(0.03, 0.985, "Dataset generations and the question each one serves", fontsize=7.5, fontweight="bold",
            transform=ax.transAxes)
    gens = [
        ("2026-07", "Project_v3 discovery", "31 GSE · 773 samples\n1,963,745 cells", "first pan-cancer CoVarNet", "#7c93ad"),
        ("2026-08 → 09", "5M merge (0906 base)", "54 GSE · 43 cancers\n1,322 samples · 4.68 M cells",
         "K=9 CMs, LIANA, V7 identity", "#c2532b"),
        ("2026-10-10", "Pan-cancer CM v20261010", "172 subclusters\n1,290 libraries · 9 CMs",
         "re-labelled discovery", "#d9a441"),
        ("2026-10-14", "Immune atlas v1", "9 sources · 400 libraries\n2,670,057 cells",
         "ICI / immune-enriched CMs", "#3b6ea8"),
        ("2026-10-23", "Full-context archive", "65 source keys\n9,439,479 observations",
         "unified identity crosswalk", "#7b5ea7"),
        ("2026-10-30", "Qualified TME scope", "43 sources · 770 units\n2,968,368 records",
         "non-enriched TME modules", "#2a7f62"),
    ]
    xs = np.linspace(0.02, 0.84, len(gens))
    ax.plot([0.02, 0.98], [0.62, 0.62], color=LIGHT, lw=2, transform=ax.transAxes, zorder=0)
    for (d, t, n, role, col), x in zip(gens, xs):
        ax.scatter([x + 0.07], [0.62], s=40, color=col, transform=ax.transAxes, zorder=3)
        ax.text(x + 0.07, 0.74, d, ha="center", fontsize=6.2, fontweight="bold", color=col, transform=ax.transAxes)
        box(ax, x, 0.08, 0.14, 0.44, f"{t}\n\n{n}\n\n{role}", fc="#fbfbfc", ec=col, fs=5.6)
    ax.text(0.5, 0.86, "broader, mixed capture  →  qualified by library, capture material and patient",
            ha="center", fontsize=5.8, color="#5b6573", transform=ax.transAxes)

    md = V5F / "metadata_audit_20260829" / "visualization_r_20260829" / "tables"
    an = pd.read_csv(md / "anatomical_source_composition.tsv", sep="\t").sort_values("n_cells")
    ax = fig.add_axes([0.14, 0.06, 0.32, 0.36])
    ax.barh(range(len(an)), an.n_cells / 1e3, color="#7c93ad", height=0.7)
    ax.set_yticks(range(len(an)), an.anatomical_source_group, fontsize=5.0)
    ax.set_xlabel("Cells (thousands)")
    ax.set_title("5M merge: anatomical source", loc="left")
    panel(ax, "b", x=-0.42)

    imm = pd.DataFrame({
        "source": ["GSE176021", "GSE169246", "CAR product (CXG)", "GSE286827", "GSE305078", "GSE294109", "GSE316195",
                   "GSE292046", "GSE256326"],
        "cells": [868494, 489360, 396593, 375019, 252822, 130243, 82549, 52525, 22452],
        "libraries": [110, 78, 59, 57, 32, 25, 22, 8, 9],
        "endpoint_cells": [0, 484326, 224628, 0, 0, 0, 0, 0, 0],
    })
    export(imm, "figS1_immune_atlas_sources.tsv")
    ax = fig.add_axes([0.62, 0.06, 0.34, 0.36])
    imm = imm.iloc[::-1]
    ax.barh(range(len(imm)), imm.cells / 1e3, color="#9cc3d5", height=0.7, label="all cells")
    ax.barh(range(len(imm)), imm.endpoint_cells / 1e3, color="#3b6ea8", height=0.7, label="source-internal endpoint-eligible")
    for i, (c, l) in enumerate(zip(imm.cells, imm.libraries)):
        ax.text(c / 1e3 + 10, i, f"{l} lib", va="center", fontsize=5.0)
    ax.set_yticks(range(len(imm)), imm.source, fontsize=5.4)
    ax.set_xlabel("Cells (thousands)")
    ax.set_xlim(0, 1050)
    ax.legend(frameon=False, fontsize=5.0, loc="lower right")
    ax.set_title("Immune atlas v1: sources", loc="left")
    panel(ax, "c", x=-0.35)
    save(fig, "figS1_dataset_generations")


def figS2_liana() -> None:
    rep = pd.read_csv(CHAIN / "focus_axes_replication_summary.tsv", sep="\t")
    sweep = pd.read_csv(CHAIN / "collagen_ddr_sweep_summary.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 5.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[0.85, 1.15], wspace=1.0, hspace=0.45)

    ax = fig.add_subplot(gs[0, 0])
    stages = [("Per-sample LIANA rows\n(within shared CM top-10)", 11_210_000),
              ("rank05 events\n(magnitude & specificity ≤ 0.05)", 382_000),
              ("Consensus\n(CM, sender, receiver, LR)", 29_605),
              ("Research priorities\n≥10 eligible, ≥3 GSE, ≥3 cancers", 139),
              ("Focus axes\n(CM01 / 07 / 08)", 12)]
    y = np.arange(len(stages))[::-1]
    vals = np.log10([s[1] for s in stages])
    ax.barh(y, vals, color=["#c9d6e3", "#9cc3d5", "#5d97b8", "#2c6e8f", "#14324a"], height=0.62)
    for yi, (lab, v) in zip(y, stages):
        ax.text(0.1, yi, f"{v:,}", va="center", fontsize=5.8, color="white" if v < 1e5 else INK, fontweight="bold")
    ax.set_yticks(y, [s[0] for s in stages], fontsize=5.4)
    ax.set_xlabel("log10 count")
    ax.set_title("From per-sample LIANA to focus axes", loc="left")
    panel(ax, "a", x=-0.75)

    ax = fig.add_subplot(gs[1, :])
    xs_ = rep["support_frac_standard_0.05"].values
    ys_ = rep["robust_strict_vs_std"].values
    ax.scatter(xs_, ys_, s=np.sqrt(rep.n_eligible_samples) * 2.2, color=[CM_COL[c] for c in rep.cm], alpha=0.8, lw=0.4,
               ec="white", zorder=3)
    order = np.argsort(ys_)
    slots = np.linspace(0.38, 1.02, len(order))
    for slot, i in zip(slots, order):
        r = rep.iloc[i]
        ax.plot([xs_[i], 1.08], [ys_[i], slot], color=CM_COL[r.cm], lw=0.35, alpha=0.6)
        ax.text(1.09, slot, f"{r.cm} {r.ligand}–{r.receptor} ({short_state(r.source)}→{short_state(r.target)})",
                fontsize=4.6, color=CM_COL[r.cm], va="center")
    ax.axhline(0.6, color=LIGHT, lw=0.7); ax.axvline(0.8, color=LIGHT, lw=0.7)
    ax.set_xlabel("Support fraction (rank ≤ 0.05)")
    ax.set_ylabel("Strict / standard retention")
    ax.set_xlim(0.5, 1.45); ax.set_ylim(0.35, 1.05)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    ax.spines["bottom"].set_bounds(0.5, 1.0)
    ax.text(0.52, 0.37, "lower-left: threshold-sensitive", fontsize=4.9, color="#5b6573")
    ax.set_title("Breadth vs threshold robustness", loc="left")
    panel(ax, "b", x=-0.1)

    ax = fig.add_subplot(gs[0, 1])
    sw = sweep.sort_values("n_samples").tail(8)
    lab = [f"{a}–{b} → {short_state(t)}" for a, b, t in zip(sw.ligand_complex, sw.receptor_complex, sw.target)]
    ax.barh(range(len(sw)), sw.n_samples, color=["#c2532b" if r == "DDR1" else "#e8b49a" for r in sw.receptor_complex],
            height=0.65)
    for i, (n, g) in enumerate(zip(sw.n_samples, sw.n_gse)):
        ax.text(n + 3, i, f"{g} GSE", va="center", fontsize=5.0)
    ax.set_yticks(range(len(sw)), lab, fontsize=5.2)
    ax.set_xlabel("Supporting samples")
    ax.set_title("Collagen–DDR pairs present\n(only COL1A1 in resource)", loc="left", fontsize=6.4)
    panel(ax, "c", x=-0.9)
    save(fig, "figS2_liana_cm")


def figS3_annotation() -> None:
    tiers = pd.read_csv(V7 / "malignancy_tier_v7_by_cancer.tsv", sep="\t")
    w = pd.read_csv(V5F / "tables" / "cm_membership_weights.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 3.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.5, 1], wspace=0.75)

    ax = fig.add_subplot(gs[0, 0])
    cats = ["Strict major", "Non-strict major", "Strict subtype", "Non-strict /\ndescriptive subtype"]
    vals = [86.76, 13.24, 52.50, 47.50]
    ax.bar([0, 1, 2.6, 3.6], vals, color=["#2c6e8f", "#9cc3d5", "#2c6e8f", "#9cc3d5"], width=0.8)
    for x, v in zip([0, 1, 2.6, 3.6], vals):
        ax.text(x, v + 1.5, f"{v:.1f}%", ha="center", fontsize=5.8)
    ax.set_xticks([0, 1, 2.6, 3.6], cats, rotation=35, ha="right", fontsize=5.4)
    ax.set_ylabel("% of 4,676,787 cells (V7)")
    ax.set_ylim(0, 100)
    ax.set_title("V7 eligibility layers", loc="left")
    panel(ax, "a", x=-0.38)

    ax = fig.add_subplot(gs[0, 1])
    tiers = tiers.set_index(tiers.columns[0])
    num = tiers.select_dtypes("number")
    keep = [c for c in ["high_confidence_malignant", "cluster_supported_malignant", "sensitive_union_only", "uncertain",
                        "cnv_no_malignancy_detected", "not_evaluable"] if c in num.columns]
    cand = num[keep]
    frac = cand.div(cand.sum(axis=1), axis=0)
    frac = frac[cand.sum(axis=1) > 2000].sort_values("high_confidence_malignant")
    cols = ["#8c1d18", "#d9603b", "#f0a35e", "#e9d8a6", "#5fa58a", "#b0b7c1"]
    left = np.zeros(len(frac))
    for c, col in zip(keep, cols):
        ax.barh(range(len(frac)), frac[c] * 100, left=left, color=col, height=0.75, lw=0, label=c.replace("_", " "))
        left += frac[c].values * 100
    ax.set_yticks(range(len(frac)), [textwrap.shorten(c, 30) for c in frac.index], fontsize=4.6)
    ax.set_xlabel("% of CNV candidates")
    ax.set_xlim(0, 100)
    ax.legend(frameon=False, fontsize=4.9, ncol=3, loc="upper left", bbox_to_anchor=(-0.4, -0.13), handlelength=0.8)
    ax.set_title("V7 CNV evidence among candidates", loc="left")
    panel(ax, "b", x=-0.62)
    export(frac.reset_index(), "figS3_v7_malignancy_fraction_by_cancer.tsv")

    ax = fig.add_subplot(gs[0, 2])
    t10 = w[w.top10].copy()
    t10["class"] = np.where(t10.subCluster.str.contains("Malignan|uncertain|not_evaluable"), "malignancy-split node",
                            np.where(t10.subCluster.str.contains("unresolved") | t10.subCluster.isin(["Unknown", "Ambiguous"]),
                                     "unresolved / Unknown", "named fine state"))
    comp = t10.groupby(["cm", "class"]).weight.sum().unstack(fill_value=0)
    comp = comp.div(comp.sum(axis=1), axis=0)[["named fine state", "malignancy-split node", "unresolved / Unknown"]]
    cms = comp.index[::-1]
    left = np.zeros(len(cms))
    for c, col in zip(comp.columns, ["#2c6e8f", "#d16d8e", "#c9ced6"]):
        ax.barh(range(len(cms)), comp.loc[cms, c] * 100, left=left, color=col, height=0.7, label=c)
        left += comp.loc[cms, c].values * 100
    ax.set_yticks(range(len(cms)), cms)
    for t in ax.get_yticklabels():
        if t.get_text() in CM_COL:
            t.set_color(CM_COL[t.get_text()]); t.set_fontweight("bold")
    ax.set_xlabel("% of top-10 member weight")
    ax.legend(frameon=False, fontsize=5.0, loc="upper left", bbox_to_anchor=(-0.2, -0.12))
    ax.set_title("K=9 member node types", loc="left")
    panel(ax, "c", x=-0.3)
    export(comp.reset_index(), "figS3_cm_member_node_class.tsv")
    save(fig, "figS3_annotation_layers")


def nmf_kl(V, k, iters=400, seed=0, theta=0.0):
    """Minimal (ns)NMF with KL multiplicative updates; S smooths H (Pascual-Montano 2006)."""
    rng = np.random.default_rng(seed)
    n, m = V.shape
    W = rng.random((n, k)) + 0.1
    H = rng.random((k, m)) + 0.1
    S = (1 - theta) * np.eye(k) + theta / k * np.ones((k, k))
    eps = 1e-9
    for _ in range(iters):
        WS = W @ S
        H *= (WS.T @ (V / (WS @ H + eps))) / (WS.sum(0)[:, None] + eps)
        SH = S @ H
        W *= ((V / (W @ SH + eps)) @ SH.T) / (SH.sum(1)[None, :] + eps)
    return W, S @ H


def simulate_rare(seed: int):
    """12 subtypes x 200 samples: 3 true modules (S1-S9) + 3 rare noise-only subtypes (S10-S12)."""
    rng = np.random.default_rng(100 + seed)
    tW = np.zeros((12, 3))
    tW[0:3, 0] = [1, 0.8, 0.6]
    tW[3:6, 1] = [1, 0.7, 0.5]
    tW[6:9, 2] = [1, 0.9, 0.4]
    X = tW @ rng.gamma(1.5, 1.0, (3, 200)) + rng.gamma(0.4, 0.15, (12, 200))
    X[9:12] = rng.poisson(0.3, (3, 200)) * rng.uniform(0.5, 2, (3, 200)) * 0.15
    Xm = (X - X.min(1, keepdims=True)) / (X.max(1, keepdims=True) - X.min(1, keepdims=True))
    return X, Xm


def figS4_nsnmf() -> None:
    fits = {"raw": [], "minmax": []}
    stats = []
    for seed in range(10):
        X, Xm = simulate_rare(seed)
        for key, M in (("raw", X), ("minmax", Xm)):
            W, _ = nmf_kl(M + 1e-6, 3, seed=seed, theta=0.5)
            Wn = W / W.sum(0, keepdims=True)
            fits[key].append(Wn)
            top3 = np.argsort(-Wn, axis=0)[:3]
            stats.append({"seed": seed, "input": key, "rare_weight_share": float(Wn[9:].sum(0).mean()),
                          "rare_in_top3": int(np.isin(top3, [9, 10, 11]).sum())})
    st = pd.DataFrame(stats)
    export(st, "figS4_simulation_rare_nodes.tsv")

    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(2, 4, width_ratios=[1.25, 0.7, 0.7, 0.85], height_ratios=[1.35, 0.8], wspace=0.95,
                          hspace=0.45)
    ax = fig.add_subplot(gs[0, 0])
    ax.axis("off")
    panel(ax, "a", x=-0.05, y=1.0)
    ax.text(0.0, 0.97, "nsNMF in CoVarNet", fontsize=7, fontweight="bold", transform=ax.transAxes)
    txt = ("X: subtype × sample frequency\n   (within parent lineage; NA→0;\n   per-subtype min–max)\n\n"
           "X ≈ W · S · H  (KL loss)\n"
           "S = (1−θ)I + (θ/K)·11ᵀ\n\n"
           "W → top-N members per CM\n"
           "H → sample activity\n"
           "K ← cophenetic heuristic\n"
           "edges ← member correlation")
    ax.text(0.0, 0.86, txt, fontsize=5.8, va="top", transform=ax.transAxes, linespacing=1.35)
    for k, key in enumerate(["raw", "minmax"]):
        ax = fig.add_subplot(gs[0, 1 + k])
        im = ax.imshow(fits[key][0], cmap=SEQ, aspect="auto", vmin=0, vmax=0.4)
        ax.set_yticks(range(12), [f"S{i + 1}" + ("*" if i >= 9 else "") for i in range(12)], fontsize=5.4)
        ax.set_xticks(range(3), [f"F{i + 1}" for i in range(3)])
        ax.axhline(8.5, color="#c2532b", lw=0.8, ls="--")
        ax.set_title({"raw": "raw freq.", "minmax": "min–max"}[key], loc="left", fontsize=6.4)
        panel(ax, "bc"[k], x=-0.55, y=1.02)
    cb = fig.colorbar(im, ax=ax, fraction=0.1, pad=0.12)
    cb.set_label("W (col.-norm.)", fontsize=5.2); cb.ax.tick_params(labelsize=5)
    fig.text(0.33, 0.405, "* S10–S12: rare subtypes carrying only noise; K = 3 true modules", fontsize=5.4, color="#5b6573")

    ax = fig.add_subplot(gs[0, 3])
    for j, key in enumerate(["raw", "minmax"]):
        v = st[st.input == key].rare_weight_share * 100
        ax.scatter(np.full(len(v), j) + np.random.default_rng(j).normal(0, 0.06, len(v)), v, s=8,
                   color=["#7c93ad", "#c2532b"][j], lw=0)
        ax.plot([j - 0.22, j + 0.22], [v.mean()] * 2, color=INK, lw=0.9)
        ax.text(j, v.max() + 2, f"mean {v.mean():.1f}%", ha="center", fontsize=5.4)
    ax.set_xticks([0, 1], ["raw", "min–max"])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, 45)
    ax.set_ylabel("Rare-node share of W (%)")
    ax.set_title("10 simulations", loc="left", fontsize=6.4)
    panel(ax, "d", x=-0.6, y=1.02)

    ax = fig.add_subplot(gs[1, 1:])
    data = [("Annotation scheme\n(CT vs HPCA, K4)", 0.715, 0.090),
            ("Source-equal weights\n(CT, K4)", 0.920, 0.514),
            ("Equal cell budget\n(CT, K4)", np.nan, 0.472),
            ("Seed\n(T/NK–myeloid, K2)", np.nan, 0.999)]
    y = np.arange(len(data))[::-1]
    ax.barh(y + 0.18, [d[1] for d in data], height=0.34, color="#2c6e8f", label="basis cosine (W)")
    ax.barh(y - 0.18, [d[2] for d in data], height=0.34, color="#c2532b", label="activity correlation (H)")
    for yi, d in zip(y, data):
        for off, v in ((0.18, d[1]), (-0.18, d[2])):
            if np.isfinite(v):
                ax.text(v + 0.02, yi + off, f"{v:.3f}", va="center", fontsize=5.0)
    ax.set_yticks(y, [d[0] for d in data], fontsize=5.2)
    ax.set_xlim(0, 1.2)
    ax.set_xlabel("Median similarity to reference fit")
    ax.legend(frameon=False, fontsize=5.0, loc="center right")
    ax.set_title("Real fits (2026-10-30 TME): same theme, different sample ranking", loc="left", fontsize=6.4)
    panel(ax, "e", x=-0.25, y=1.02)
    save(fig, "figS4_nsnmf_limits")


def figS5_spatial_first() -> None:
    fig = plt.figure(figsize=(7.2, 2.9))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    panel(ax, "a", x=0.01, y=0.92)
    ax.text(0.04, 0.935, "Current order (0906)", fontsize=7, fontweight="bold", transform=ax.transAxes)
    cur = ["scRNA modules\n(sample frequency)", "LIANA within CM\n(all LR, all pairs)", "Pick axes", "Spatial: score\nprojection / LR z",
           "Bulk proxy"]
    xs = np.linspace(0.04, 0.80, len(cur))
    for t, x in zip(cur, xs):
        box(ax, x, 0.66, 0.15, 0.2, t, fc="#f4f5f7", ec="#b7c0cb", fs=5.9)
    for a, b in zip(xs[:-1], xs[1:]):
        arrow(ax, a + 0.152, 0.76, b - 0.002, 0.76, color=GREY)
    panel(ax, "b", x=0.01, y=0.48)
    ax.text(0.04, 0.495, "Proposed order: module → spatial niche → communication", fontsize=7, fontweight="bold",
            transform=ax.transAxes)
    new = [("Frozen module W\n(state vocabulary)", "#eaf4ef"),
           ("Map states on single-cell\nspatial data (Xenium/CosMx/\nVisium HD) via label transfer", "#eaf4ef"),
           ("Local module activity\nin 25/50/100 µm windows\nvs density-matched null", "#fbefe9"),
           ("Niche contact graph:\nwhich members actually\ntouch / co-localise", "#fbefe9"),
           ("LR restricted to\ncontact pairs + receiver\nresponse (TF / programme)", "#eef3f8")]
    xs = np.linspace(0.04, 0.80, len(new))
    for (t, fc), x in zip(new, xs):
        box(ax, x, 0.08, 0.15, 0.3, t, fc=fc, ec="#b7c0cb", fs=5.6)
    for a, b in zip(xs[:-1], xs[1:]):
        arrow(ax, a + 0.152, 0.23, b - 0.002, 0.23, color=GREY)
    ax.text(0.04, 0.01, "Unit of inference throughout: patient / section (not spot or cell); replication across "
            "independent sections and cancers.", fontsize=5.6, color="#5b6573", transform=ax.transAxes)
    save(fig, "figS5_spatial_first_design")


def figS6_improvement() -> None:
    scopes = pd.DataFrame({
        "scope": ["whole_TME", "immune_enriched", "CAR_product", "T_lineage_biopsy", "whole_TME_primary_tumor",
                  "T_whole_TME_primary_tumor", "whole_TME_metastatic_tumor", "immune_enriched_primary_tumor",
                  "T_immune_enriched_primary_tumor", "immune_enriched_blood", "T_immune_enriched_blood",
                  "immune_enriched_normal_adjacent", "T_immune_enriched_normal_adjacent",
                  "immune_enriched_unknown_material", "T_immune_enriched_unknown_material"],
        "rank": [2, 9, 12, 12, 12, 4, 2, 3, 3, 5, 5, 3, 3, 2, 2],
        "supported": [2, 9, 12, 12, 12, 4, 0, 3, 3, 0, 0, 3, 2, 0, 0],
        "edges": [70, 51, 36, 10, 95, 7, 0, 13, 2, 0, 0, 15, 3, 1, 0],
    })
    export(scopes, "figS6_immune_atlas_cm_scopes.tsv")
    fig = plt.figure(figsize=(7.2, 3.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.55)
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(scopes))[::-1]
    ax.barh(y + 0.18, scopes["rank"], height=0.34, color="#9cc3d5", label="NMF rank")
    ax.barh(y - 0.18, scopes.supported, height=0.34, color="#2c6e8f", label="CMs with network-supported members")
    for yi, e in zip(y, scopes.edges):
        ax.text(12.6, yi, f"{e}", va="center", fontsize=5.2, color="#5b6573")
    ax.text(12.6, len(scopes) - 0.3, "edges", fontsize=5.2, color="#5b6573")
    ax.set_yticks(y, scopes.scope.str.replace("_", " "), fontsize=5.2)
    ax.set_xlim(0, 14)
    ax.set_xlabel("Count")
    ax.legend(frameon=False, fontsize=5.2, loc="upper center", bbox_to_anchor=(0.4, -0.12))
    ax.set_title("Immune atlas v1: modules per material scope", loc="left")
    panel(ax, "a", x=-0.62)

    ax = fig.add_subplot(gs[0, 1]); ax.axis("off")
    panel(ax, "b", x=-0.05, y=1.0)
    ax.text(0.0, 0.98, "Three layers to compare", fontsize=7, fontweight="bold", transform=ax.transAxes)
    layers = [
        ("Composition CM\n(non-epithelial nodes only)", "immune–stromal–vascular frequency;\nmalignant programme as context", "#eaf4ef"),
        ("Composition CM\n(+ malignant / epithelial nodes)", "tumour-specific nodes per cancer;\nCNV-qualified, split by cancer", "#fbefe9"),
        ("State layer inside identity\n(immune datasets)", "patient × identity × programme:\nmean, upper quantile, active fraction", "#eef3f8"),
    ]
    for i, (t, s, fc) in enumerate(layers):
        yy = 0.66 - i * 0.3
        box(ax, 0.0, yy + 0.08, 0.9, 0.14, t.replace("\n", " "), fc=fc, ec="#b7c0cb", fs=5.8, bold=True)
        ax.text(0.02, yy + 0.03, s.replace("\n", " "), fontsize=5.5, va="center", transform=ax.transAxes,
                color="#475260")
    ax.text(0.0, 0.0, "Compare on common samples, same K and seeds;\nreport matched / split / merged modules.",
            fontsize=5.6, color="#5b6573", transform=ax.transAxes)
    save(fig, "figS6_cm_improvement")


def manifest() -> None:
    import hashlib

    rows = []
    for p in sorted(OUT.glob("*.png")):
        rows.append({"file": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest(), "bytes": p.stat().st_size})
    (OUT / "figure_manifest.json").write_text(json.dumps({"generated_by": "scripts/render_homepage_v2_figures.py",
                                                          "figures": rows}, indent=2))


if __name__ == "__main__":
    if not RAW.exists() or str(RAW) == "":
        raise SystemExit("Set RAW to the pulled server table directory")
    for f in (fig1, fig2, fig3, fig4, fig5, fig6, fig7, figS1_datasets, figS2_liana, figS3_annotation, figS4_nsnmf,
              figS5_spatial_first, figS6_improvement):
        f()
    manifest()
