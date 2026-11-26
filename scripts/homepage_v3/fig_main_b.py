"""Main figures 4-7: K=9 CM landscape, within-CM communication, spatial and bulk readouts."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle
from scipy.cluster.hierarchy import linkage, leaves_list

from style import (REPO, SRC2, INK, MUTED, GRID, CM_COL, OTHER_CM, LINEAGE_COL, DIV, SEQ, panel, export, save, cmc, box, arrow)

CMS = [f"CM0{i}" for i in range(1, 10)]


def _cm_ticks(ax, axis="x"):
    labs = ax.get_xticklabels() if axis == "x" else ax.get_yticklabels()
    for t in labs:
        t.set_color(cmc(t.get_text()))
        if t.get_text() in CM_COL:
            t.set_fontweight("bold")


def _short(s: str) -> str:
    rep = {"Monocyte_Classical_FCN1_FCN1": "Mono FCN1", "Macrophage_SPP1_TREM2_SPP1": "Mac SPP1",
           "Fibroblast_myCAF_ACTA2": "myCAF ACTA2", "Malignant::Epithelial": "Malignant epi",
           "Monocyte_Inflammatory_IL1B_IL1B": "Infl. Mono IL1B", "Macrophage_Inflammatory_IL1B_IL1B": "Infl. Mac IL1B",
           "NK_XCL1_CD56bright_XCL1": "NK XCL1", "Bcell_Naive_IGHD": "B naive IGHD", "Endothelial_Blood_PECAM1": "Endo PECAM1",
           "Endothelial_Venous_ACKR1_ACKR1": "Endo ACKR1", "Tcell_CD8_cytotoxic_CCL5": "CD8 CCL5",
           "Mast_Activated_IL1RL1_IL1RL1": "Mast IL1RL1", "Dendritic_pDC_LILRA4": "pDC LILRA4",
           "Bcell_Naive_TCL1A": "B naive TCL1A", "NK_FCGR3A_cytotoxic_FCGR3A": "NK FCGR3A", "Tcell_CD4_naive_IL7R": "CD4 naive IL7R",
           "Tcell_Cycling_T_MKI67": "T cycling MKI67", "Platelet_coarse": "Platelet (coarse)"}
    if s in rep:
        return rep[s]
    s = s.replace("_unresolved", " (unres.)").replace("Malignancy_uncertain::", "Mal. uncertain::").replace("_", " ")
    return s


# ---------------------------------------------------------------------------
def fig4() -> None:
    mem = pd.read_csv(SRC2 / "fig4_cm_top10_membership.tsv", sep="\t")
    lin = pd.read_csv(SRC2 / "fig4_cm_lineage_composition.tsv", sep="\t").set_index("cm")
    act = pd.read_csv(SRC2 / "fig4_cm_cancer_activity_ge10.tsv", sep="\t")
    inf = pd.read_csv(SRC2 / "fig4_cm_infiltration_associations.tsv", sep="\t")
    ncls = pd.read_csv(SRC2 / "figS3_cm_member_node_class.tsv", sep="\t").set_index("cm")

    fig = plt.figure(figsize=(7.2, 8.6))
    # a: member dot matrix
    ax = fig.add_axes([0.20, 0.52, 0.30, 0.45])
    top = mem[mem["rank"] <= 6].copy()
    # row order: by CM with maximal weight, then weight
    best = mem.loc[mem.groupby("subCluster").weight_fraction.idxmax(), ["subCluster", "cm", "weight_fraction"]]
    rows = [r for r in best.sort_values(["cm", "weight_fraction"], ascending=[True, False]).subCluster if r in set(top.subCluster)]
    full = mem[mem.subCluster.isin(rows)]
    yi = {r: i for i, r in enumerate(rows)}
    for _, r in full.iterrows():
        x = CMS.index(r.cm); y = yi[r.subCluster]
        c = LINEAGE_COL.get(r.majorCluster, "#999")
        s = r.weight_fraction * 900
        if r.review_flag:
            ax.scatter(x, y, s=s, facecolor="white", edgecolor=c, lw=0.9, zorder=3)
        else:
            ax.scatter(x, y, s=s, color=c, lw=0, zorder=3)
    for k, cm in enumerate(CMS):
        if cm in CM_COL:
            ax.axvspan(k - 0.5, k + 0.5, color=CM_COL[cm], alpha=0.07, lw=0, zorder=0)
    ax.set_xticks(range(9), CMS, rotation=90); _cm_ticks(ax)
    ax.set_yticks(range(len(rows)), [_short(r) for r in rows], fontsize=4.9)
    ax.set_xlim(-0.6, 8.6); ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False)
    for s_, lab in [(0.05, "5%"), (0.1, "10%"), (0.2, "20%")]:
        ax.scatter([], [], s=s_ * 900, color="#8a94a3", label=lab)
    ax.scatter([], [], s=90, facecolor="white", edgecolor="#8a94a3", lw=0.9, label="review-flagged\n(unresolved)")
    ax.legend(title="share of module W", title_fontsize=5.2, loc="upper left", bbox_to_anchor=(1.02, 0.35),
              labelspacing=1.1, borderpad=0.2)
    ax.set_title("Top-6 members of each CM (union; dot = share of W in that CM)", fontsize=6.4, pad=26)
    panel(ax, "a", x=-0.62, y=1.07)

    # b: lineage composition + node class
    ax = fig.add_axes([0.70, 0.73, 0.28, 0.22])
    order = ["Tcell", "NK", "ILC", "Bcell", "Plasma", "Macrophage", "Monocyte", "Dendritic", "Mast", "Granulocyte",
             "Platelet", "Erythroid", "HSPC", "Fibroblast", "Endothelial", "Epithelial", "Neural", "Unknown", "Ambiguous"]
    for i, cm in enumerate(CMS[::-1]):
        left = 0
        for l in order:
            v = lin.loc[cm].get(l, 0) * 100
            if v > 0:
                ax.barh(i, v, left=left, color=LINEAGE_COL[l], height=0.72, lw=0.3, ec="white")
                left += v
        fl = (1 - ncls.loc[cm, "named fine state"]) * 100
        ax.text(102, i, f"{fl:.0f}%", va="center", fontsize=5, color=MUTED)
    ax.text(101, 8.85, "not fine-\nnamed", fontsize=4.5, color=MUTED, va="bottom")
    ax.set_yticks(range(9), CMS[::-1]); _cm_ticks(ax, "y")
    ax.set_xlim(0, 100); ax.set_xlabel("% of top-10 member weight")
    ax.set_title("Lineage make-up", fontsize=6.4)
    hs = [Rectangle((0, 0), 1, 1, fc=LINEAGE_COL[l]) for l in order]
    ax.legend(hs, order, ncol=4, fontsize=4.6, bbox_to_anchor=(-0.2, -0.28), loc="upper left", handlelength=0.8,
              columnspacing=0.6, handletextpad=0.3)
    panel(ax, "b", x=-0.22)

    # c: activity heatmap (clustered cancers)
    ax = fig.add_axes([0.20, 0.06, 0.27, 0.38])
    piv = act.pivot_table(index="cancer_type", columns="cm", values="median_activity")[CMS]
    z = (piv - piv.mean()) / piv.std(ddof=0)
    order_r = leaves_list(linkage(z.fillna(0).values, "average", metric="correlation"))
    z = z.iloc[order_r]
    im = ax.imshow(z.values, cmap=DIV, norm=TwoSlopeNorm(0, -2.5, 2.5), aspect="auto")
    ax.set_xticks(range(9), CMS, rotation=90); _cm_ticks(ax)
    lab = [c if len(c) < 30 else c[:28] + "…" for c in z.index]
    ax.set_yticks(range(len(z)), lab, fontsize=4.6)
    ax.set_title("Median activity by cancer label (≥10 samples)", fontsize=6.4)
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03, shrink=0.5); cb.set_label("z across cancers", fontsize=5.4)
    cb.ax.tick_params(labelsize=5)
    for s in ax.spines.values():
        s.set_visible(False)
    panel(ax, "c", x=-0.8)

    # d: infiltration associations
    ax = fig.add_axes([0.66, 0.06, 0.25, 0.38])
    feats = ["malignant_strict_fraction", "malignancy_uncertain_fraction", "immune_hot_proxy", "lymphoid_fraction",
             "cytotoxic_fraction", "inflammatory_myeloid_fraction", "spp1_macrophage_fraction", "myeloid_fraction",
             "caf_fraction", "stromal_activation_proxy"]
    labmap = inf.drop_duplicates("feature").set_index("feature").feature_label.to_dict()
    feats = [f for f in feats if f in labmap]
    for i, f in enumerate(feats):
        for j, cm in enumerate(CMS):
            r = inf[(inf.feature == f) & (inf.cm == cm)]
            if r.empty:
                continue
            r = r.iloc[0]
            rho = r.within_cancer_spearman
            sig = (r.within_cancer_fdr < 0.05) and abs(rho) >= 0.2
            ax.scatter(j, i, s=abs(rho) * 260 + 4, c=[rho], cmap=DIV, norm=TwoSlopeNorm(0, -0.7, 0.7),
                       edgecolor=INK if sig else "none", lw=0.7)
    ax.set_xticks(range(9), CMS, rotation=90); _cm_ticks(ax)
    ax.set_yticks(range(len(feats)), [labmap[f].replace(" fraction", "") for f in feats], fontsize=5.4)
    ax.set_xlim(-0.6, 8.6); ax.set_ylim(len(feats) - 0.5, -0.5)
    ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
    sm = plt.cm.ScalarMappable(cmap=DIV, norm=TwoSlopeNorm(0, -0.7, 0.7))
    cb = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.03, shrink=0.5); cb.set_label("within-cancer Spearman ρ", fontsize=5.4)
    cb.ax.tick_params(labelsize=5)
    ax.set_title("CM activity vs infiltration (n = 1,303)\nblack ring: FDR < 0.05 and |ρ| ≥ 0.2", fontsize=6.4)
    panel(ax, "d", x=-0.6)
    export(inf[inf.feature.isin(feats)][["cm", "feature_label", "within_cancer_spearman", "within_cancer_fdr"]],
           "fig4d_infiltration_plotted.tsv")
    save(fig, "fig4_cm_landscape", ["inventories/publication/2026-11-02/homepage_v2/fig4_*.tsv", "figS3_cm_member_node_class.tsv"])


# ---------------------------------------------------------------------------
def _axis_label(a: str) -> tuple[str, str]:
    p = a.split("_")
    cm, lig, rec = p[0], p[1], p[2]
    rest = " ".join(p[3:]).replace(" to ", " → ").replace("ITGA1B1", "")
    if lig == "COL1A1" and rec == "ITGA1B1":
        rec = "ITGA1:ITGB1"; rest = " ".join(p[3:]).replace(" to ", " → ")
    return cm, f"{lig}→{rec}  ({rest})"


def fig5() -> None:
    rep = pd.read_csv(SRC2 / "fig5_focus_axes_replication.tsv", sep="\t")
    byc = pd.read_csv(SRC2 / "fig5_focus_axes_by_cancer.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 6.4))

    # a: schematic of three CMs and their axes
    ax = fig.add_axes([0.01, 0.62, 0.40, 0.36]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    panel(ax, "a", x=0.0, y=0.98)
    ax.text(0.06, 0.985, "Focus axes in CM01 / 07 / 08", fontsize=6.6, va="bottom")
    # node x positions (axes coords) per CM row
    rows = {
        "CM01": (0.80, [("Mono FCN1", 0.08), ("Mac SPP1", 0.32), ("Endo", 0.56)],
                 [(0, 1, "S100A9→CD68", 0), (1, 2, "C1QA→CD93", 0.03), (2, 1, "APP→CD74", -0.03)],
                 "FCN1 mono → SPP1 mac\n→ endothelium"),
        "CM07": (0.50, [("Infl. Mono", 0.08), ("myCAF", 0.32), ("Malignant epi", 0.56)],
                 [(1, 0, "COL1A1→CD44", 0), (1, 2, "COL1A1→DDR1", 0)],
                 "myCAF → malignant\nepithelium"),
        "CM08": (0.18, [("myCAF", 0.08), ("CD8 T", 0.32), ("Infl. Mono", 0.56)],
                 [(0, 1, "LGALS1→PTPRC\nMYL9→CD69\nCOL1A1→CD44", 0), (2, 1, "S100A8→CD69", -0.0001)],
                 "inflammatory myeloid /\nmyCAF → CD8 T"),
    }
    for cm, (y, nd, ed, desc) in rows.items():
        c = CM_COL[cm]
        for (i, j, lab, dy) in ed:
            x0, x1 = nd[i][1], nd[j][1]
            sgn = 1 if x1 > x0 else -1
            ax.annotate("", xy=(x1 - sgn * 0.045, y + dy), xytext=(x0 + sgn * 0.045, y + dy),
                        arrowprops=dict(arrowstyle="-|>", color=c, lw=0.9, mutation_scale=7))
            ty = y + dy + (0.035 if dy >= 0 else -0.035)
            ax.text((x0 + x1) / 2, ty, lab, ha="center", va="bottom" if dy >= 0 else "top", fontsize=4.6, color=c,
                    linespacing=1.15)
        for (n, x) in nd:
            ax.add_patch(Circle((x, y), 0.035, fc=c, ec="white", lw=0.6, alpha=0.9))
            ax.text(x, y - 0.065, n, ha="center", fontsize=4.9, color=INK, va="top")
        ax.text(0.70, y + 0.02, cm, fontsize=7, fontweight="bold", color=c, va="center")
        ax.text(0.70, y - 0.04, desc, fontsize=4.8, color=MUTED, va="top")
    ax.text(0.32, 0.41, "+ COL1A1→DDR2 and COL1A1→ITGA1:ITGB1 (myCAF autocrine)", fontsize=4.4, color=CM_COL["CM07"], ha="center")
    ax.text(0.45, 0.005, "CM08 also: S100A9→CD68 (Infl. Mono → Infl. Mac)", fontsize=4.5, color=CM_COL["CM08"], ha="center")

    # b: support at three thresholds (dumbbell)
    ax = fig.add_axes([0.58, 0.62, 0.40, 0.33])
    r = rep.copy()
    r["cm"] = r.axis.str[:4]
    r = r.sort_values(["cm", "support_frac_standard_0.05"], ascending=[True, False]).reset_index(drop=True)
    y = np.arange(len(r))[::-1]
    for yi_, (_, row) in zip(y, r.iterrows()):
        c = cmc(row.cm)
        v = [row["support_frac_strict_0.01"], row["support_frac_standard_0.05"], row["support_frac_loose_0.10"]]
        ax.plot([v[0], v[2]], [yi_, yi_], color=c, lw=1.0, alpha=0.5)
        ax.scatter(v[0], yi_, marker="o", facecolor="white", edgecolor=c, s=14, lw=0.8, zorder=3)
        ax.scatter(v[1], yi_, marker="o", color=c, s=18, zorder=4)
        ax.scatter(v[2], yi_, marker="|", color=c, s=30, zorder=3)
        ax.text(1.03, yi_, f"{int(row['n_samples_standard_0.05'])}/{int(row.n_eligible_samples)}", fontsize=4.7,
                va="center", color=MUTED)
    ax.set_yticks(y, [_axis_label(a)[1] for a in r.axis], fontsize=4.7)
    for t, cm in zip(ax.get_yticklabels(), r.cm):
        t.set_color(cmc(cm))
    ax.set_xlim(0.0, 1.0); ax.set_xlabel("Fraction of testable samples supporting the axis")
    ax.grid(axis="x", color=GRID, lw=0.5); ax.set_axisbelow(True)
    ax.legend([Line2D([], [], marker="o", ls="", mfc="white", mec=MUTED), Line2D([], [], marker="o", ls="", color=MUTED),
               Line2D([], [], marker="|", ls="", color=MUTED, ms=6)],
              ["rank ≤ 0.01", "rank ≤ 0.05", "rank ≤ 0.10"], ncol=3, loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=5)
    ax.text(1.03, len(r) - 0.3, "n/N", fontsize=4.8, color=MUTED)
    panel(ax, "b", x=-0.72, y=1.07)

    # c: breadth across GSE and cancer labels
    ax = fig.add_axes([0.07, 0.08, 0.30, 0.43])
    for _, row in r.iterrows():
        ax.scatter(row["n_gse_standard_0.05"], row["n_cancer_standard_0.05"], s=row["n_samples_standard_0.05"] * 0.35,
                   color=cmc(row.cm), alpha=0.55, ec="white", lw=0.5)
    lbl = {"CM07_COL1A1_DDR1_myCAF_to_Malignant": "COL1A1→DDR1", "CM01_S100A9_CD68_FCN1Mono_to_SPP1Mac": "S100A9→CD68 (CM01)",
           "CM08_LGALS1_PTPRC_myCAF_to_CD8T": "LGALS1→PTPRC", "CM01_C1QA_CD93_SPP1Mac_to_Endo": "C1QA→CD93"}
    offs = {"CM07_COL1A1_DDR1_myCAF_to_Malignant": (5, -12), "CM01_S100A9_CD68_FCN1Mono_to_SPP1Mac": (-6, 8),
            "CM08_LGALS1_PTPRC_myCAF_to_CD8T": (-60, -12), "CM01_C1QA_CD93_SPP1Mac_to_Endo": (6, -4)}
    for _, row in r.iterrows():
        if row.axis in lbl:
            ax.annotate(lbl[row.axis], (row["n_gse_standard_0.05"], row["n_cancer_standard_0.05"]), fontsize=4.9,
                        xytext=offs[row.axis], textcoords="offset points", color=cmc(row.cm))
    ax.set_xlabel("GSEs with ≥1 supporting sample"); ax.set_ylabel("Cancer labels with ≥1 supporting sample")
    for s_ in [50, 150, 300]:
        ax.scatter([], [], s=s_ * 0.35, color="#b9c2cc", label=f"{s_}")
    ax.legend(title="supporting samples", title_fontsize=5, loc="center right", labelspacing=1.0)
    ax.set_xlim(21, 43); ax.set_ylim(20.5, 35.5)
    ax.set_title("Breadth of support (rank ≤ 0.05)", fontsize=6.4)
    ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
    panel(ax, "c", x=-0.22)

    # d: by-cancer supporting samples heatmap
    ax = fig.add_axes([0.50, 0.08, 0.40, 0.43])
    piv = byc.pivot_table(index="axis", columns="cancer_type", values="n_samples", aggfunc="sum").fillna(0)
    topc = piv.sum().sort_values(ascending=False).head(18).index
    piv = piv.loc[r.axis, topc]
    im = ax.imshow(np.log10(piv.values + 1), cmap=SEQ, aspect="auto", vmin=0, vmax=np.log10(60))
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = int(piv.values[i, j])
            if v:
                ax.text(j, i, v, ha="center", va="center", fontsize=4.0, color="white" if v > 15 else INK)
    ax.set_yticks(range(len(piv)), [_axis_label(a)[1].split("  ")[0] for a in piv.index], fontsize=4.8)
    for t, cm in zip(ax.get_yticklabels(), r.cm):
        t.set_color(cmc(cm))
    ax.set_xticks(range(len(topc)), [c if len(c) < 24 else c[:22] + "…" for c in topc], rotation=60, ha="right", fontsize=4.6)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, shrink=0.6)
    cb.set_ticks(np.log10([1, 3, 11, 31])); cb.set_ticklabels(["0", "2", "10", "30"]); cb.set_label("supporting samples", fontsize=5.2)
    cb.ax.tick_params(labelsize=5)
    ax.set_title("Supporting samples per cancer label (18 labels with most support)", fontsize=6.4)
    panel(ax, "d", x=-0.2)
    export(r[["axis", "cm", "n_eligible_samples", "support_frac_strict_0.01", "support_frac_standard_0.05",
              "support_frac_loose_0.10", "n_gse_standard_0.05", "n_cancer_standard_0.05", "top_gse_share"]], "fig5_axes_plotted.tsv")
    save(fig, "fig5_within_cm_communication", ["inventories/publication/2026-11-02/homepage_v2/fig5_*.tsv"])


# ---------------------------------------------------------------------------
CANCER = {"GSE175540": "ccRCC", "GSE203612": "breast", "GSE238264": "HCC", "GSE273952": "ccRCC", "GSE301720": "HNSCC"}
NSL = {"GSE175540": 24, "GSE203612": 10, "GSE238264": 7, "GSE273952": 12, "GSE301720": 8}


def fig6() -> None:
    cs = pd.read_csv(SRC2 / "fig6_visium_cross_slice_summary.tsv", sep="\t")
    mi = pd.read_csv(SRC2 / "fig6_visium_cm_morans_i.tsv", sep="\t")
    t3 = pd.read_csv(SRC2 / "fig6_cosmx_T3.tsv", sep="\t")
    t3b = pd.read_csv(SRC2 / "fig6_cosmx_T3b.tsv", sep="\t")
    t4 = pd.read_csv(SRC2 / "fig6_cosmx_T4.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(2, 6, height_ratios=[1.15, 1], hspace=0.75, wspace=1.6)

    ax = fig.add_subplot(gs[0, :3])
    piv = cs.pivot_table(index="axis_id", columns="dataset_id", values="median_permutation_z")
    piv = piv.loc[piv.median(axis=1).sort_values(ascending=False).index]
    fp = cs.pivot_table(index="axis_id", columns="dataset_id", values="fraction_positive_z").loc[piv.index, piv.columns]
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            ax.scatter(j, i, s=min(abs(v), 11) * 22 + 6, c=[v], cmap=DIV, norm=TwoSlopeNorm(0, -10, 10), lw=0)
            ax.text(j + 0.32, i + 0.28, f"{fp.values[i, j]*100:.0f}%", fontsize=3.9, color=MUTED)
    ax.set_xticks(range(piv.shape[1]), [f"{d}\n{CANCER[d]} · {NSL[d]}" for d in piv.columns], fontsize=4.8)
    ylab = []
    for a in piv.index:
        p = a.split("_")
        ylab.append(f"{p[0]}  {p[1]}→{p[2]}")
    ax.set_yticks(range(len(piv)), ylab, fontsize=5.2)
    for t, a in zip(ax.get_yticklabels(), piv.index):
        t.set_color(cmc(a[:4]))
    ax.set_xlim(-0.5, piv.shape[1] - 0.2); ax.set_ylim(len(piv) - 0.5, -0.6)
    ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
    sm = plt.cm.ScalarMappable(cmap=DIV, norm=TwoSlopeNorm(0, -10, 10))
    cax = ax.inset_axes([0.62, -0.26, 0.36, 0.035])
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal"); cb.set_label("median permutation z", fontsize=5.0)
    cb.ax.tick_params(labelsize=4.8)
    ax.set_title("Visium: ligand×sender vs receptor×receiver co-enrichment\n(dot = median z over sections; % = sections with z > 0)", fontsize=6.2)
    panel(ax, "a", x=-0.3)

    ax = fig.add_subplot(gs[0, 3:])
    dsets = list(CANCER)
    pal = ["#2f63a3", "#d4a017", "#2a7f62", "#7aa6d6", "#9e7b9b"]
    for k, cm in enumerate(["CM01", "CM07", "CM08"]):
        for j, d in enumerate(dsets):
            v = mi[(mi.CM == cm) & (mi.dataset_id == d)].morans_I.values
            x = k * 6 + j
            parts = ax.violinplot([v], positions=[x], widths=0.85, showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(pal[j]); b.set_alpha(0.35); b.set_edgecolor("none")
            ax.scatter(np.full(len(v), x) + np.random.default_rng(j).uniform(-0.18, 0.18, len(v)), v, s=3, color=pal[j], lw=0)
            ax.plot([x - 0.3, x + 0.3], [np.median(v)] * 2, color=INK, lw=0.8)
    ax.set_xticks([2, 8, 14], ["CM01", "CM07", "CM08"]); _cm_ticks(ax)
    ax.set_ylabel("Moran's I of projected CM score"); ax.set_ylim(0, 1)
    ax.legend([Rectangle((0, 0), 1, 1, fc=p, alpha=0.6) for p in pal], [f"{d} ({CANCER[d]})" for d in dsets], ncol=2,
              fontsize=4.8, loc="lower left", bbox_to_anchor=(0, -0.42))
    ax.set_title("Spatial autocorrelation of module scores\n(183 section×CM tests; all permutation p ≤ 0.001)", fontsize=6.2)
    panel(ax, "b", x=-0.18)

    ax = fig.add_subplot(gs[1, :2])
    nl = pd.read_csv(REPO / "inventories/analysis/2026-11-06/CM07_DDR1_exploration_v1/A_T3b_density_null.tsv", sep="\t")
    export(nl[["tissue", "n_cd8", "frac_hi", "ratio_obs", "ratio_null", "ratio_null_fov", "p_vs_null_fov"]], "fig6c_T3b_density_null.tsv")
    x = np.arange(len(nl))
    ax.bar(x - 0.26, nl.ratio_obs, width=0.25, color="#c4502a", label="observed")
    ax.bar(x, nl.ratio_null, width=0.25, color="#d9dee5", label="density null (global)")
    ax.bar(x + 0.26, nl.ratio_null_fov, width=0.25, color="#8fa6bf", label="density null (within FOV)")
    for xi, r in zip(x, nl.itertuples()):
        ax.text(xi - 0.26, r.ratio_obs + 0.12, f"{r.ratio_obs:.2f}", ha="center", fontsize=4.6)
        ax.text(xi + 0.26, r.ratio_null_fov + 0.12, f"{r.ratio_null_fov:.2f}", ha="center", fontsize=4.6, color=MUTED)
    ax.axhline(1, color=MUTED, lw=0.5, ls="--")
    ax.set_xticks(x, nl.tissue); ax.set_ylabel("CD8 distance ratio\n(collagen×DDR1-high / -low tumour)"); ax.set_ylim(0, 6.6)
    ax.legend(fontsize=4.6, loc="upper left", handlelength=1)
    ax.set_title("CosMx ccRCC T3b, with density null\n(2026-11-06 re-check)", fontsize=6.0)
    panel(ax, "c", x=-0.42)

    for k, (tab, title, ylab, letter) in enumerate([(t3, "T3: CD8 neighbours of tumour cells", "β (per SD)", "d"),
                                                    (t4, "T4: tumour EMT/MMP programme", "β (per SD)", "e")]):
        ax = fig.add_subplot(gs[1, 2 + k * 2: 4 + k * 2])
        x = np.arange(len(tab))
        for off, col, cl, nm in [(-0.25, "beta_col_exp", "#2a7f62", "local collagen"), (0, "beta_ddr1", "#c4502a", "tumour DDR1"),
                                 (0.25, "beta_interaction", "#8a94a3", "interaction")]:
            ax.bar(x + off, tab[col], width=0.23, color=cl, label=nm)
        ax.axhline(0, color=INK, lw=0.5)
        ax.set_xticks(x, tab.tissue); ax.set_ylabel(ylab); ax.set_title(title, fontsize=6.0)
        if k == 0:
            ax.legend(fontsize=4.8, loc="lower left")
        panel(ax, letter, x=-0.4)
    save(fig, "fig6_spatial_readouts", ["inventories/publication/2026-11-02/homepage_v2/fig6_*.tsv"])


# ---------------------------------------------------------------------------
def fig7() -> None:
    st = pd.read_csv(SRC2 / "fig7_tcga_stratified_cox.tsv", sep="\t")
    per = pd.read_csv(SRC2 / "fig7_tcga_per_cancer_cox.tsv", sep="\t")
    resp = pd.read_csv(SRC2 / "fig7_ici_response_tests.tsv", sep="\t")
    surv = pd.read_csv(SRC2 / "fig7_ici_survival_tests.tsv", sep="\t")
    prox = pd.read_csv(SRC2 / "fig7_tcga_lr_proxy_models.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 7.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.8, 1], width_ratios=[1, 1.3, 1.05], hspace=0.6, wspace=1.25)

    # a: pan-cancer forest
    ax = fig.add_subplot(gs[0, 0])
    sigs = CMS + ["desmoplasia_myCAF", "TLS_cabrita12", "TLS_niche3"]
    d = st[st.endpoint == "OS"].set_index("signature").loc[sigs]
    y = np.arange(len(d))[::-1]
    for yi_, (s, r) in zip(y, d.iterrows()):
        c = cmc(s) if s.startswith("CM") else "#7b5ea7"
        ax.plot([r.CI_low, r.CI_high], [yi_, yi_], color=c, lw=1.0)
        ax.scatter(r.HR, yi_, color=c, s=16 if s in CM_COL else 9, zorder=3)
    ax.axvline(1, color=MUTED, lw=0.5, ls="--")
    ax.set_yticks(y, [s.replace("_", " ") for s in d.index], fontsize=5.4); _cm_ticks(ax, "y")
    ax.set_xlabel("HR per SD (OS)"); ax.set_title("TCGA pan-cancer, stratified by type\n9,034 tumours · 2,664 deaths", fontsize=6.2)
    panel(ax, "a", x=-0.55)

    # b: per-cancer heatmap for CM01/07/08 OS
    ax = fig.add_subplot(gs[0, 1:])
    pp = per[(per.endpoint == "OS") & (per.signature.isin(["CM01", "CM07", "CM08"]))]
    piv = pp.pivot_table(index="signature", columns="cancer", values="HR")
    q = pp.pivot_table(index="signature", columns="cancer", values="BH_q")
    cols = piv.columns[np.argsort(np.log(piv.loc["CM01"]).values)]
    piv = piv[cols]; q = q[cols]
    im = ax.imshow(np.log2(piv.values), cmap=DIV, norm=TwoSlopeNorm(0, -1, 1), aspect="auto")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            if q.values[i, j] < 0.05:
                ax.scatter(j, i, s=6, color=INK)
    ax.set_yticks(range(3), piv.index); _cm_ticks(ax, "y")
    ax.set_xticks(range(len(cols)), cols, rotation=90, fontsize=5)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label("log2 HR (OS)", fontsize=5.2); cb.ax.tick_params(labelsize=5)
    ax.set_title("Direction within each TCGA cancer type (dot: BH q < 0.05 within signature)", fontsize=6.2)
    panel(ax, "b", x=-0.06)

    # c: ICI response
    ax = fig.add_subplot(gs[1, 0])
    cohorts = ["Riaz2017_melanoma_all", "Riaz2017_melanoma_preonly", "Hugo2016_melanoma", "IMvigor210_urothelial"]
    sel = ["CM01", "CM07", "CM08", "desmoplasia_myCAF"]
    rr = resp[resp.signature.isin(sel)]
    pv = rr.pivot_table(index="signature", columns="cohort", values="cliffs_delta").loc[sel, cohorts]
    qv = rr.pivot_table(index="signature", columns="cohort", values="BH_q").loc[sel, cohorts]
    nr = rr.drop_duplicates("cohort").set_index("cohort")
    im = ax.imshow(pv.values, cmap=DIV, norm=TwoSlopeNorm(0, -0.5, 0.5), aspect="auto")
    for i in range(pv.shape[0]):
        for j in range(pv.shape[1]):
            ax.text(j, i, f"{pv.values[i, j]:+.2f}" + ("*" if qv.values[i, j] < 0.1 else ""), ha="center", va="center", fontsize=4.8)
    ax.set_yticks(range(len(sel)), [s.replace("_", " ") for s in sel]); _cm_ticks(ax, "y")
    ax.set_xticks(range(len(cohorts)), [f"{c.replace('_melanoma', '').replace('_urothelial', '').replace('_', ' ')}\n{int(nr.loc[c,'n_R'])} R / {int(nr.loc[c,'n_NR'])} NR" for c in cohorts],
                  rotation=45, ha="right", fontsize=4.8)
    ax.set_title("ICI response, Cliff's δ\n(R vs NR; * BH q < 0.1)", fontsize=6.2)
    panel(ax, "c", x=-0.55)

    # d: ICI survival forest
    ax = fig.add_subplot(gs[1, 1])
    ss = surv[surv.signature.isin(["CM01", "CM07", "CM08"])].copy()
    ss["lab"] = ss.cohort.str.replace("_", " ") + " · " + ss.signature
    ss = ss.sort_values(["cohort", "signature"]).reset_index(drop=True)
    y = np.arange(len(ss))[::-1]
    for yi_, r in zip(y, ss.itertuples()):
        c = cmc(r.signature)
        ax.plot([r.CI_low, r.CI_high], [yi_, yi_], color=c, lw=1.0)
        ax.scatter(r.HR, yi_, color=c, s=12, zorder=3)
        ax.text(1.02, yi_, f"n={r.n} q={r.BH_q:.2f}", fontsize=4.2, va="center", color=MUTED,
                transform=ax.get_yaxis_transform())
    ax.axvline(1, color=MUTED, lw=0.5, ls="--")
    ax.set_xscale("log"); ax.set_xlim(0.25, 5)
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xticks([0.25, 0.5, 1, 2, 4], ["0.25", "0.5", "1", "2", "4"])
    ax.set_yticks(y, [f"{r.cohort.split('_')[0]} ({r.endpoint}) · {r.signature}" for r in ss.itertuples()], fontsize=4.8)
    _cm_ticks(ax, "y")
    for t, r in zip(ax.get_yticklabels(), ss.itertuples()):
        t.set_color(cmc(r.signature))
    ax.set_xlabel("HR per SD"); ax.set_title("ICI survival (Cox)", fontsize=6.2)
    panel(ax, "d", x=-0.55)

    # e: proxy increment
    ax = fig.add_subplot(gs[1, 2])
    pr = prox[(prox.endpoint == "OS") & (prox.proxy.str.startswith(("min_", "mean_")))].copy()
    def _pl(p):
        kind, rest = p.split("_", 1)
        nm = {"collagen_DDR1": "collagen–DDR1", "COL1A1_only_DDR1": "COL1A1–DDR1", "LGALS1_myCAF_CD8": "LGALS1 CAF→CD8",
              "SPP1_macrophage": "SPP1 mac"}[rest]
        return f"{nm} ({kind})"
    pr["lab"] = pr.proxy.map(_pl)
    y = np.arange(len(pr))[::-1]
    for yi_, r in zip(y, pr.itertuples()):
        c = {"CM01": CM_COL["CM01"], "CM08": CM_COL["CM08"]}.get(r.parent, "#7b5ea7")
        ax.scatter(r.HR_proxy_in_C, yi_, color=c, s=16, zorder=3, label=None)
        ax.scatter(r.HR_parent_in_C, yi_, facecolor="white", edgecolor=c, s=14, zorder=3)
        ax.plot([r.HR_parent_in_C, r.HR_proxy_in_C], [yi_, yi_], color=c, lw=0.6, alpha=0.6)
        ax.text(1.36, yi_, f"LRT p={r.lrt_C_vs_A_p:.0e}" if r.lrt_C_vs_A_p < 1e-3 else f"LRT p={r.lrt_C_vs_A_p:.2f}",
                fontsize=4.3, va="center", color=MUTED)
    ax.axvline(1, color=MUTED, lw=0.5, ls="--")
    ax.set_yticks(y, list(pr.lab), fontsize=4.8)
    ax.set_xlim(0.8, 1.35); ax.set_xlabel("HR in model C (TCGA OS)")
    ax.legend([Line2D([], [], marker="o", ls="", color=MUTED), Line2D([], [], marker="o", ls="", mfc="white", mec=MUTED)],
              ["LR proxy", "parent signature"], fontsize=4.8, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2)
    ax.set_title("LR proxy beyond parent signature", fontsize=6.2, pad=14)
    panel(ax, "e", x=-0.75, y=1.1)
    save(fig, "fig7_bulk_ici_readouts", ["inventories/publication/2026-11-02/homepage_v2/fig7_*.tsv"])
