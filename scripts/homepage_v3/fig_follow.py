"""Figure 8 (follow-up discovery versions) and supplementary figures S1-S6 for the subpages."""
from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle

from style import (REPO, SRC2, INK, MUTED, GRID, CM_COL, LINEAGE_COL, DIV, SEQ, SEQ_WARM, panel, fig_letter, box, arrow,
                   export, save, cmc)

V4 = REPO / "inventories/analysis/2026-11-07/v4_CM_selection_v1"
R4 = REPO / "inventories/analysis/2026-11-19/WEHI_R4_CM_validation/results"
REL = REPO / "releases/communication_evidence_chain_20260906/tables"
BLUE = "#2f63a3"


# ---------------------------------------------------------------------------
# Figure 8  the same questions on the follow-up discovery versions
# ---------------------------------------------------------------------------
LEDGER_ROWS = ["Composition robustness\n(seed / source)", "Within-CM communication", "Spatial co-localisation",
               "TCGA outcome", "ICI response / survival"]
LEDGER_COLS = ["0906 K9\nCM01", "0906 K9\nCM07", "0906 K9\nCM08", "V4 K14\nCM03", "R4 K15\nCM12"]
# 3 = consistent in all tested units; 2 = mixed / context dependent; 1 = no support in pre-specified family; 0 = not tested
LEDGER = np.array([
    [0, 0, 0, 2, 2],
    [3, 2, 2, 1, 1],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
    [2, 1, 2, 1, 1],
])
LEDGER_NOTE = [
    ["no seed/LOSO refit on 0906", "no seed/LOSO refit on 0906", "no seed/LOSO refit on 0906",
     "seed W cos 0.990; LOSO edges 10/13", "13/15 fixed edges BH; no refit LOSO"],
    ["support 0.79–1.00 (abundance-like)", "COL1A1→DDR1 0.57", "LGALS1→PTPRC 0.92; strict 42%",
     "author-style response q ≥ 0.05", "LR events shared across CMs; patient LR 0 q<0.05"],
    ["C1QA→CD93 4/5; APP→CD74 3 neg", "VCAN→ITGB1 5/5; DDR1 only CosMx ccRCC", "S100A9→CD68 5/5; LGALS1 mixed",
     "ACKR1-EC neighbourhood B/NK +; triad open", "ESM1–RGS5 + in 4/4 Visium, 8/8 ST; FOLR2 hub mixed"],
    ["HR 0.97 pan; opposite by type", "HR 1.11 pan; collinear w/ desmoplasia", "HR 1.06 pan; LGG/UVM >1, SARC <1",
     "LGG OS q≈0.007 only", "KIRP model-specific; mural-driven"],
    ["Riaz all δ +0.30 (q 0.087) only", "no response / survival signal", "Rose HR 0.70 (q 0.069); IMvigor 1.00",
     "6 ICI families 0 BH", "Hugo/Riaz/Rose 480 tests 0 q<0.05"],
]


def _wrap(t, n=26):
    words, lines, cur = t.split(" "), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > n and cur:
            lines.append(cur); cur = w
        else:
            cur = (cur + " " + w).strip()
    lines.append(cur)
    return "\n".join(lines)


def fig8() -> None:
    rv = pd.read_csv(V4 / "ALL_CM_REVIEW.tsv", sep="\t")
    vis = pd.read_csv(R4 / "WEHI_R4_spatial_20261115/review/Visium_core_evidence.tsv", sep="\t")
    st = pd.read_csv(R4 / "CM_validation_20261118/integrated_evidence_review/CM12_three_edges_all_sections.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 7.6))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.15], hspace=0.85, wspace=0.75)

    # a: V4 K14 – network weight vs stability
    ax = fig.add_subplot(gs[0, 0])
    rec_col = {"首选验证候选": "#c4502a", "备选；优先查来源": "#d4a017"}
    for _, r in rv.iterrows():
        c = rec_col.get(r.recommendation, "#9aa5b1")
        ax.scatter(r.network_W_mass_fraction * 100, r.seed_W_cosine_median, s=12 + 60 * max(r.LOSO_W_cosine_median, 0) ** 4,
                   color=c, ec="white", lw=0.5, zorder=3)
        off = {"CM03": (-16, -9), "CM06": (5, 2), "CM09": (5, -2), "CM11": (-6, -10), "CM12": (5, -3), "CM05": (5, -2)}
        if r.cm in off:
            ax.annotate(r.cm, (r.network_W_mass_fraction * 100, r.seed_W_cosine_median), fontsize=4.9, xytext=off[r.cm],
                        textcoords="offset points", color=c if c != "#9aa5b1" else MUTED)
    ax.set_xlabel("% of factor W carried by its network"); ax.set_ylabel("seed W cosine (median of 5)")
    ax.set_ylim(0.1, 1.05); ax.set_xlim(-5, 105)
    ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.set_title("V4 K14 (560,871 non-epithelial cells)\n14 factors × 540 controlled refits", fontsize=6.0)
    ax.legend([Line2D([], [], marker="o", ls="", color="#c4502a"), Line2D([], [], marker="o", ls="", color="#d4a017"),
               Line2D([], [], marker="o", ls="", color="#9aa5b1")], ["first candidate", "back-up", "deferred / diagnostic"],
              fontsize=4.6, loc="lower right", markerscale=0.6, bbox_to_anchor=(1.0, 0.12))
    panel(ax, "a", x=-0.38)
    export(rv[["cm", "network_W_mass_fraction", "seed_W_cosine_median", "LOSO_W_cosine_median", "source_activity_R2",
               "recommendation"]], "fig8a_v4_k14_review.tsv")

    # b: CM03 edges
    ax = fig.add_subplot(gs[0, 1])
    e = pd.DataFrame({
        "edge": ["IGHD B – PECAM1 EC", "PECAM1 EC – XCL1 NK"],
        "raw r": [0.422, 0.410], "source-centred": [0.297, 0.227], "both parents ≥20 + centred": [0.414, 0.298],
        "LOSO min": [0.330, 0.348], "LOSO max": [0.504, 0.526], "within-source +": ["12/13", "10/13"]})
    export(e, "fig8b_v4_cm03_edges.tsv")
    metrics = ["raw r", "source-centred", "both parents ≥20 + centred"]
    cols = ["#9aa5b1", "#c4502a", "#7b5ea7"]
    for i, r in e.iterrows():
        for k, (m, c) in enumerate(zip(metrics, cols)):
            ax.barh(i * 1.2 + k * 0.28, r[m], height=0.26, color=c)
        ax.plot([r["LOSO min"], r["LOSO max"]], [i * 1.2 - 0.25] * 2, color=INK, lw=1.2)
        ax.text(r["LOSO max"] + 0.01, i * 1.2 - 0.25, "LOSO range", fontsize=4.5, va="center", color=MUTED)
        ax.text(0.58, i * 1.2 + 0.28, f"source +\n{r['within-source +']}", fontsize=4.6, va="center", color=MUTED)
    ax.set_yticks([0.28, 1.48], e.edge, fontsize=5.4); ax.set_xlim(0, 0.68); ax.invert_yaxis()
    ax.set_xlabel("Pearson r across samples")
    ax.legend([Rectangle((0, 0), 1, 1, fc=c) for c in cols], metrics, fontsize=4.6, loc="upper left", bbox_to_anchor=(-0.05, -0.32), ncol=2)
    ax.set_title("V4 CM03 = IGHD B – PECAM1 EC – XCL1 NK\n(no direct B–NK edge)", fontsize=6.0)
    panel(ax, "b", x=-0.6)

    # c: R4 Visium core edges
    ax = fig.add_subplot(gs[0, 2])
    vis["edge"] = vis.cm + " " + vis.node1.str.split("_").str[-1] + "–" + vis.node2.str.split("_").str[-1]
    order = ["CM12 ESM1–RGS5", "CM12 ESM1–FOLR2", "CM12 FOLR2–RGS5", "CM04 GJA5–MYH11", "CM01 ACKR1–PI16"]
    secs = list(vis.section.unique())
    mk = {"breast": "o", "ccRCC": "s"}
    for i, ed in enumerate(order):
        d = vis[vis.edge == ed]
        for j, r in enumerate(d.itertuples()):
            ax.scatter(r.partial_rank_total_and_epithelial_fraction, i + (j - 1.5) * 0.12, marker=mk[r.cancer], s=14,
                       color="#c4502a" if ed.startswith("CM12 ESM1–RGS5") else "#8a94a3", ec="white", lw=0.3)
    ax.axvline(0, color=MUTED, lw=0.5, ls="--")
    ax.set_yticks(range(len(order)), order, fontsize=5.2); ax.invert_yaxis()
    ax.set_xlabel("partial Spearman (total & epithelial adjusted)")
    ax.legend([Line2D([], [], marker="o", ls="", color=MUTED), Line2D([], [], marker="s", ls="", color=MUTED)],
              ["breast (2 sections)", "ccRCC (2 sections)"], fontsize=4.6, loc="upper left", bbox_to_anchor=(-0.05, -0.32), ncol=2, markerscale=0.7)
    ax.set_title("R4 K15 on Visium (4 sections, 10,160 spots)", fontsize=6.0)
    panel(ax, "c", x=-0.75)
    export(vis[vis.edge.isin(order)][["section", "cancer", "edge", "partial_rank_total_and_epithelial_fraction"]],
           "fig8c_r4_visium_edges.tsv")

    # d: legacy ST 8 sections, three CM12 edges
    ax = fig.add_subplot(gs[1, :2])
    st = st[st.posterior == "means"].copy()
    st["edge"] = st.node1.str.split("_").str[-1] + "–" + st.node2.str.split("_").str[-1]
    edges = ["ESM1–RGS5", "ESM1–FOLR2", "FOLR2–RGS5"]
    mets = [("raw_spearman", "raw"), ("partial_rank_total_melanoma_other_nuisance", "partial (nuisance adj.)"),
            ("parent_fraction_spearman", "within-parent fraction")]
    secs = sorted(st.section.unique())
    pal = plt.cm.tab10(np.linspace(0, 0.8, 8))
    for gi, ed in enumerate(edges):
        for mi, (col, lab) in enumerate(mets):
            x0 = gi * 4 + mi
            v = st[st.edge == ed].set_index("section")[col].reindex(secs)
            ax.scatter(np.full(8, x0) + np.linspace(-0.25, 0.25, 8), v, c=pal, s=10, zorder=3)
            ax.plot([x0 - 0.3, x0 + 0.3], [v.median()] * 2, color=INK, lw=0.9)
    ax.axhline(0, color=MUTED, lw=0.5, ls="--")
    ax.set_xticks([gi * 4 + mi for gi in range(3) for mi in range(3)], [m[1] for _ in edges for m in mets], rotation=30,
                  ha="right", fontsize=4.8)
    for gi, ed in enumerate(edges):
        ax.text(gi * 4 + 1, 1.02, ed, ha="center", fontsize=6, fontweight="bold",
                color="#c4502a" if gi == 0 else INK, transform=ax.get_xaxis_transform())
    ax.set_ylabel("Spearman across spots"); ax.set_ylim(-0.6, 1.0)
    ax.legend([Line2D([], [], marker="o", ls="", color=p, ms=4) for p in pal], [s.replace("_N30", "") for s in secs], ncol=4,
              fontsize=4.4, loc="lower left", bbox_to_anchor=(0, -0.62))
    ax.set_title("R4 CM12 core and hub edges across 8 melanoma legacy-ST sections (N = 30, posterior means)", fontsize=6.0, pad=12)
    panel(ax, "d", x=-0.1, y=1.08)
    export(st[["section", "edge", "raw_spearman", "partial_rank_total_melanoma_other_nuisance", "parent_fraction_spearman"]],
           "fig8d_r4_legacy_st.tsv")

    # e: version lineage
    ax = fig.add_subplot(gs[1, 2]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    vers = [("0906 · K9", "V5 identity · all lineages · 1,303 samples", "#c4502a"),
            ("10-10 · K9", "172 subclusters · 1,290 libraries", "#9aa5b1"),
            ("10-14 · immune atlas", "15 material scopes", "#7b5ea7"),
            ("10-30 · eligible TME", "43 sources · 4 panels", "#2a7f62"),
            ("11-06 · V4 K14", "non-epithelial · 164 samples · 45 nodes", BLUE),
            ("11-15 · R4 K15", "44 identities · 1,194 units", BLUE)]
    for i, (h, t, c) in enumerate(vers):
        y = 0.92 - i * 0.165
        ax.add_patch(Circle((0.06, y), 0.025, fc=c, ec="white"))
        if i < len(vers) - 1:
            ax.plot([0.06, 0.06], [y - 0.03, y - 0.135], color=GRID, lw=1.2)
        ax.text(0.13, y + 0.035, h, fontsize=5.6, fontweight="bold", color=c, va="center")
        ax.text(0.13, y - 0.035, t, fontsize=4.5, color=MUTED, va="center")
    ax.text(0.05, 1.08, "CM numbering is version-specific", fontsize=6.0, transform=ax.transAxes)
    panel(ax, "e", x=-0.1, y=1.08)

    # f: ledger
    ax = fig.add_subplot(gs[2, :])
    cmap = ListedColormap(["#eef0f3", "#d9dee5", "#f2d6b8", "#c4502a"])
    ax.imshow(LEDGER, cmap=cmap, vmin=-0.5, vmax=3.5, aspect="auto")
    for i in range(LEDGER.shape[0]):
        for j in range(LEDGER.shape[1]):
            ax.text(j, i, _wrap(LEDGER_NOTE[i][j]), ha="center", va="center", fontsize=4.5, linespacing=1.15,
                    color="white" if LEDGER[i, j] == 3 else INK, wrap=True)
    for k in range(1, 5):
        ax.axvline(k - 0.5, color="white", lw=2); ax.axhline(k - 0.5, color="white", lw=2)
    ax.set_xticks(range(5), LEDGER_COLS, fontsize=5.6)
    ax.set_yticks(range(5), LEDGER_ROWS, fontsize=5.6)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.legend([Rectangle((0, 0), 1, 1, fc=c) for c in ["#c4502a", "#f2d6b8", "#d9dee5", "#eef0f3"]],
              ["same direction in all tested units", "mixed / context-dependent", "no support in the pre-specified family",
               "not tested / not evaluable"], ncol=4, fontsize=5, loc="lower left", bbox_to_anchor=(0, -0.3))
    ax.set_title("Readout ledger: what each analysis layer recorded for the five highlighted modules (text = key number; see Section 9 and subpages)",
                 fontsize=6.0)
    panel(ax, "f", x=-0.12, y=1.04)
    export(pd.DataFrame([(LEDGER_ROWS[i].replace("\n", " "), LEDGER_COLS[j].replace("\n", " "), int(LEDGER[i, j]), LEDGER_NOTE[i][j])
                         for i in range(5) for j in range(5)], columns=["layer", "module", "category", "note"]),
           "fig8f_readout_ledger.tsv")
    save(fig, "fig8_followup_versions", ["inventories/analysis/2026-11-07/v4_CM_selection_v1/ALL_CM_REVIEW.tsv", "R4 Visium_core_evidence.tsv",
                                          "R4 CM12_three_edges_all_sections.tsv", "readout ledger = README sections 6–9"])


# ---------------------------------------------------------------------------
# S1 datasets
# ---------------------------------------------------------------------------
def figS1() -> None:
    src = pd.read_csv(SRC2 / "figS1_immune_atlas_sources.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 3.6))
    ax = fig.add_axes([0.02, 0.1, 0.5, 0.82]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig_letter(fig, 0.005, 0.93, "a")
    ax.text(0, 1.0, "Material decides what a proportion means", fontsize=6.8, va="top")
    cats = [("Non-enriched tumour tissue", "whole-TME composition and states", "#c4502a"),
            ("Enriched / sorted material\n(CD45⁺, CD3⁺, CAR product, blood)", "states within the enriched lineage only", "#7b5ea7"),
            ("Spatial data", "tissue validation of measurable identities", "#2a7f62"),
            ("Clinical bulk cohorts", "independent outcome / response evaluation", BLUE)]
    for i, (h, t, c) in enumerate(cats):
        y = 0.78 - i * 0.22
        box(ax, 0.0, y - 0.08, 0.5, 0.17, h, fc="white", ec=c, fs=5.4, color=c, weight="bold")
        arrow(ax, 0.51, y, 0.57, y, color=c)
        ax.text(0.59, y, t, fontsize=5.4, va="center", color=INK)
    ax2 = fig.add_axes([0.62, 0.16, 0.35, 0.72])
    fig_letter(fig, 0.56, 0.93, "b")
    s = src.sort_values("cells")
    y = np.arange(len(s))
    ax2.barh(y, s.cells / 1e3, color="#d0d6de", height=0.7, label="all cells")
    ax2.barh(y, s.endpoint_cells / 1e3, color="#7b5ea7", height=0.7, label="endpoint-eligible")
    for yi, r in zip(y, s.itertuples()):
        ax2.text(r.cells / 1e3 + 10, yi, f"{r.libraries} lib.", fontsize=4.8, va="center", color=MUTED)
    ax2.set_yticks(y, s.source, fontsize=5.4); ax2.set_xlabel("Cells (×10³)")
    ax2.legend(fontsize=5, loc="lower right")
    ax2.set_title("Immune atlas v1: 9 sources, 2,670,057 cells, 400 libraries", fontsize=6.2)
    export(src, "figS1_immune_atlas_sources.tsv")
    save(fig, "figS1_datasets", ["inventories/publication/2026-11-02/homepage_v2/figS1_immune_atlas_sources.tsv"])


# ---------------------------------------------------------------------------
# S2 LIANA
# ---------------------------------------------------------------------------
def figS2() -> None:
    rep = pd.read_csv(SRC2 / "fig5_focus_axes_replication.tsv", sep="\t")
    sw = pd.read_csv(REL / "collagen_ddr_sweep_summary.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 3.3))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.1, 0.8], wspace=1.25)
    ax = fig.add_subplot(gs[0])
    st = pd.DataFrame({"step": ["within-CM rows", "rank ≤ 0.05 events", "consensus LR", "priority candidates", "focus axes"],
                       "n": [11.21e6, 382e3, 29605, 139, 12]})
    export(st, "figS2_funnel.tsv")
    ax.barh(range(5)[::-1], st.n, color=["#c9d6e4", "#9cbad6", "#5f8fbf", "#2f63a3", "#c4502a"], height=0.6)
    ax.set_xscale("log"); ax.set_xlim(5, 5e7)
    for i, r in zip(range(5)[::-1], st.itertuples()):
        ax.text(r.n * 1.3, i, f"{r.n:,.0f}" if r.n < 1e5 else f"{r.n/1e6:.2f} M" if r.n > 1e6 else f"{r.n/1e3:.0f} k",
                va="center", fontsize=5.2)
    ax.set_yticks(range(5)[::-1], st.step, fontsize=5.4); ax.set_xlabel("count (log)")
    ax.set_title("From 1,301 per-sample runs to 12 axes", fontsize=6.2)
    panel(ax, "a", x=-0.7)

    ax = fig.add_subplot(gs[1])
    r = rep.copy(); r["cm"] = r.axis.str[:4]
    r["keep"] = r["n_samples_strict_0.01"] / r["n_samples_standard_0.05"]
    for _, row in r.iterrows():
        ax.scatter(row["support_frac_standard_0.05"], row.keep, color=cmc(row.cm), s=18, ec="white", lw=0.4)
        p = row.axis.split("_")
        if row["support_frac_standard_0.05"] < 0.96 or row.keep < 0.85:
            ax.annotate(f"{p[1]}→{p[2]}", (row["support_frac_standard_0.05"], row.keep), fontsize=4.4, xytext=(-3, 3),
                        textcoords="offset points", color=cmc(row.cm), ha="right")
    ax.text(0.985, 0.88, "CM01 axes,\nS100A9→CD68 (CM08)", fontsize=4.4, ha="right", color=MUTED)
    ax.set_xlabel("support fraction (rank ≤ 0.05)"); ax.set_ylabel("kept at rank ≤ 0.01 / kept at 0.05")
    ax.set_xlim(0.5, 1.05); ax.set_ylim(0.35, 1.05); ax.grid(color=GRID, lw=0.4); ax.set_axisbelow(True)
    ax.set_title("Breadth vs threshold robustness", fontsize=6.2)
    panel(ax, "b", x=-0.35)

    ax = fig.add_subplot(gs[2])
    sw = sw.sort_values("n_samples", ascending=True)
    lab = [f"{a}→{b} · {t.replace('_', ' ').replace('Malignant::Epithelial', 'Malignant epi')[:18]}" for a, b, t in zip(sw.ligand_complex, sw.receptor_complex, sw.target)]
    c = ["#c4502a" if (a == "COL1A1" and b == "DDR1" and "Malignant" in t) else "#8fa6bf"
         for a, b, t in zip(sw.ligand_complex, sw.receptor_complex, sw.target)]
    ax.barh(range(len(sw)), sw.n_samples, color=c, height=0.7)
    ax.set_yticks(range(len(sw)), lab, fontsize=4.4); ax.set_xlabel("supporting samples")
    ax.set_title("Collagen–DDR pairs present in the resource\n(DDR1 is paired only with COL1A1)", fontsize=6.0)
    panel(ax, "c", x=-1.3)
    save(fig, "figS2_liana_cm", ["fig5_focus_axes_replication.tsv", "releases/communication_evidence_chain_20260906/tables/collagen_ddr_sweep_summary.tsv"])


# ---------------------------------------------------------------------------
# S3 annotation
# ---------------------------------------------------------------------------
def figS3() -> None:
    tiers = pd.read_csv(SRC2 / "figS3_v7_malignancy_fraction_by_cancer.tsv", sep="\t").set_index("cancer_type")
    ncls = pd.read_csv(SRC2 / "figS3_cm_member_node_class.tsv", sep="\t").set_index("cm")
    fig = plt.figure(figsize=(7.2, 3.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.6, 0.9, 0.9], wspace=0.7)
    ax = fig.add_subplot(gs[0])
    cols = ["high_confidence_malignant", "cluster_supported_malignant", "sensitive_union_only", "uncertain",
            "cnv_no_malignancy_detected", "not_evaluable"]
    pal = ["#8c2f18", "#d0784e", "#f0b48f", "#e9d8a6", "#5fa58a", "#c9ced6"]
    t = tiers.sort_values("high_confidence_malignant").tail(24)
    left = np.zeros(len(t))
    for c, p in zip(cols, pal):
        ax.barh(range(len(t)), t[c] * 100, left=left, color=p, height=0.75, lw=0)
        left += t[c].values * 100
    ax.set_yticks(range(len(t)), [x if len(x) < 30 else x[:28] + "…" for x in t.index], fontsize=4.6)
    ax.set_xlim(0, 100); ax.set_xlabel("% of CNV candidate cells (V7)")
    ax.legend([Rectangle((0, 0), 1, 1, fc=p) for p in pal], [c.replace("_", " ") for c in cols], ncol=3, fontsize=4.5,
              loc="lower left", bbox_to_anchor=(-0.3, -0.25))
    ax.set_title("Per-sample CNV malignancy tiers by cancer (24 labels)", fontsize=6.2)
    panel(ax, "a", x=-0.5)

    ax = fig.add_subplot(gs[1])
    cms = [f"CM0{i}" for i in range(1, 10)]
    pal2 = ["#2f63a3", "#d16d8e", "#c9ced6"]
    left = np.zeros(9)
    for c, p in zip(ncls.columns, pal2):
        ax.barh(range(9)[::-1], ncls.loc[cms, c] * 100, left=left, color=p, height=0.7)
        left += ncls.loc[cms, c].values * 100
    ax.set_yticks(range(9)[::-1], cms)
    for tl in ax.get_yticklabels():
        tl.set_color(cmc(tl.get_text()))
    ax.set_xlim(0, 100); ax.set_xlabel("% of top-10 member W")
    ax.legend([Rectangle((0, 0), 1, 1, fc=p) for p in pal2], list(ncls.columns), fontsize=4.6, loc="lower left",
              bbox_to_anchor=(-0.2, -0.3))
    ax.set_title("What the K9 nodes are", fontsize=6.2)
    panel(ax, "b", x=-0.35)

    ax = fig.add_subplot(gs[2])
    vals = [2843051, 1218888]
    ax.pie(vals, colors=["#2f63a3", "#c9d6e4"], startangle=90, counterclock=False, wedgeprops=dict(width=0.38, ec="white"))
    ax.text(0, 0, "4,061,939\ncells", ha="center", va="center", fontsize=5.6)
    ax.legend([Rectangle((0, 0), 1, 1, fc=c) for c in ["#2f63a3", "#c9d6e4"]],
              ["original subtype kept\n2,843,051", "major class + marker\n1,218,888"], fontsize=4.8, loc="lower center",
              bbox_to_anchor=(0.5, -0.32))
    ax.set_title("V7 selective subtype restore\n(V4 labels, 2026-11-06)", fontsize=6.2)
    panel(ax, "c", x=-0.25, y=1.12)
    save(fig, "figS3_annotation_layers", ["figS3_*.tsv", "inventories/analysis/2026-11-06/v7_selective_subtype_restore_v4/README.md"])


# ---------------------------------------------------------------------------
# S4 nsNMF
# ---------------------------------------------------------------------------
def nmf_kl(V, k, iters=400, seed=0, theta=0.0):
    rng = np.random.default_rng(seed)
    n, m = V.shape
    W = rng.random((n, k)) + 0.1; H = rng.random((k, m)) + 0.1
    S = (1 - theta) * np.eye(k) + theta / k * np.ones((k, k)); eps = 1e-9
    for _ in range(iters):
        WS = W @ S
        H *= (WS.T @ (V / (WS @ H + eps))) / (WS.sum(0)[:, None] + eps)
        SH = S @ H
        W *= ((V / (W @ SH + eps)) @ SH.T) / (SH.sum(1)[None, :] + eps)
    return W, S @ H


def simulate_rare(seed):
    rng = np.random.default_rng(100 + seed)
    tW = np.zeros((12, 3)); tW[0:3, 0] = [1, 0.8, 0.6]; tW[3:6, 1] = [1, 0.7, 0.5]; tW[6:9, 2] = [1, 0.9, 0.4]
    X = tW @ rng.gamma(1.5, 1.0, (3, 200)) + rng.gamma(0.4, 0.15, (12, 200))
    X[9:12] = rng.poisson(0.3, (3, 200)) * rng.uniform(0.5, 2, (3, 200)) * 0.15
    Xm = (X - X.min(1, keepdims=True)) / (X.max(1, keepdims=True) - X.min(1, keepdims=True))
    return X, Xm


def figS4() -> None:
    sim = pd.read_csv(SRC2 / "figS4_simulation_rare_nodes.tsv", sep="\t")
    rv = pd.read_csv(V4 / "ALL_CM_REVIEW.tsv", sep="\t")
    X, Xm = simulate_rare(0)
    Wr, _ = nmf_kl(X, 3, theta=0.5); Wm, _ = nmf_kl(Xm, 3, theta=0.5)
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1], hspace=0.7, wspace=0.9)
    for k, (W, t) in enumerate([(Wr, "raw frequency input"), (Wm, "min–max per subcluster")]):
        ax = fig.add_subplot(gs[0, k])
        Wn = W / W.sum(0, keepdims=True)
        ax.imshow(Wn, cmap=SEQ, aspect="auto", vmin=0, vmax=0.35)
        ax.set_yticks(range(12), [f"S{i+1}" + (" (rare)" if i >= 9 else "") for i in range(12)], fontsize=5)
        ax.set_xticks(range(3), ["F1", "F2", "F3"]); ax.axhline(8.5, color="#c4502a", lw=0.8, ls="--")
        ax.set_title(f"W, {t}", fontsize=6.0)
        panel(ax, "ab"[k], x=-0.5)
    ax = fig.add_subplot(gs[0, 2])
    for j, inp in enumerate(["raw", "minmax"]):
        v = sim[sim.input == inp].rare_weight_share * 100
        ax.scatter(np.full(len(v), j) + np.linspace(-0.15, 0.15, len(v)), v, color=["#8fa6bf", "#c4502a"][j], s=10)
        ax.plot([j - 0.25, j + 0.25], [v.mean()] * 2, color=INK, lw=1)
        ax.text(j, v.max() + 2, f"mean {v.mean():.1f}%", ha="center", fontsize=5)
    ax.set_xticks([0, 1], ["raw", "min–max"]); ax.set_xlim(-0.6, 1.6); ax.set_ylabel("% of W on 3 rare subclusters")
    ax.set_title("10 simulations", fontsize=6.0)
    panel(ax, "c", x=-0.5)
    ax = fig.add_subplot(gs[0, 3])
    diag = pd.DataFrame({"contrast": ["two annotations\n(same cells, K4)", "source-equal\nweighting", "equal cell\nbudget", "seeds\n(same panel)"],
                         "basis_cos": [0.715, 0.920, np.nan, np.nan], "activity_r": [0.090, 0.514, 0.472, 0.999]})
    export(diag, "figS4_real_diagnostics.tsv")
    y = np.arange(4)[::-1]
    ax.barh(y + 0.18, diag.basis_cos, height=0.34, color="#7b5ea7", label="basis cosine")
    ax.barh(y - 0.18, diag.activity_r, height=0.34, color="#c4502a", label="activity r")
    ax.set_yticks(y, diag.contrast, fontsize=4.8); ax.set_xlim(0, 1.05)
    ax.legend(fontsize=4.6, loc="lower right", bbox_to_anchor=(1.0, -0.42))
    ax.set_title("2026-10-30 TME refits", fontsize=6.0)
    panel(ax, "d", x=-0.75)

    ax = fig.add_subplot(gs[1, :])
    r = rv.reset_index(drop=True)
    x = np.arange(len(r))
    for xi, row in zip(x, r.itertuples()):
        ax.plot([xi - 0.15, xi - 0.15], [row.seed_W_cosine_min, row.seed_W_cosine_median], color="#2f63a3", lw=1.2)
        ax.scatter(xi - 0.15, row.seed_W_cosine_median, color="#2f63a3", s=14, zorder=3)
        ax.plot([xi + 0.15, xi + 0.15], [row.LOSO_W_cosine_min, row.LOSO_W_cosine_median], color="#c4502a", lw=1.2)
        ax.scatter(xi + 0.15, row.LOSO_W_cosine_median, color="#c4502a", s=14, zorder=3)
        ax.bar(xi, row.network_W_mass_fraction, width=0.7, color="#e6e9ee", zorder=0)
    ax.set_xticks(x, ["V4-" + c for c in r.cm], fontsize=5.4); ax.set_ylim(0, 1.05)
    ax.set_ylabel("W cosine (dot median, line to min)\n/ grey bar: network share of W")
    ax.legend([Line2D([], [], marker="o", color="#2f63a3"), Line2D([], [], marker="o", color="#c4502a"), Rectangle((0, 0), 1, 1, fc="#e6e9ee")],
              ["5 extra seeds", "13 leave-one-source-out", "network share of W"], ncol=3, fontsize=5, loc="lower left",
              bbox_to_anchor=(0, 1.0))
    ax.set_title("V4 K14: an optimisation can be stable while the network carries little of the factor", fontsize=6.2, pad=16)
    panel(ax, "e", x=-0.07, y=1.1)
    save(fig, "figS4_nsnmf_limits", ["figS4_simulation_rare_nodes.tsv (+ re-run seed 0)", "inventories/analysis/2026-11-07/v4_CM_selection_v1/ALL_CM_REVIEW.tsv"])


# ---------------------------------------------------------------------------
# S5 spatial-first
# ---------------------------------------------------------------------------
def figS5() -> None:
    imgq = pd.read_csv(R4 / "CM_validation_20261118/imaging_core_review/primary_summary.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 4.4))
    ax = fig.add_axes([0.02, 0.52, 0.96, 0.44]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig_letter(fig, 0.005, 0.95, "a")
    ax.text(0, 1.0, "0906 order (top) and proposed spatial-first order (bottom)", fontsize=6.6, va="top")
    a = ["sc CM", "within-CM LIANA\n(all member pairs)", "pick axes", "Visium / CosMx\nreadout", "bulk proxy"]
    b = ["freeze CM recipe\n(state vocabulary)", "map states on\nsingle-cell spatial", "local CM activity vs\ndensity-matched null",
         "member contact /\nco-location graph", "LR only on contact\npairs + receiver response"]
    for row, items, c in [(0.62, a, "#9aa5b1"), (0.12, b, "#2a7f62")]:
        for i, t in enumerate(items):
            box(ax, i * 0.2, row, 0.17, 0.26, t, fc="white", ec=c, fs=5.2)
            if i < 4:
                arrow(ax, i * 0.2 + 0.172, row + 0.13, i * 0.2 + 0.198, row + 0.13, color=c)
    ax.text(0.99, 0.02, "inference unit = patient or whole section; buffers ≥ neighbourhood radius", fontsize=4.8,
            color=MUTED, ha="right")
    ax2 = fig.add_axes([0.08, 0.08, 0.38, 0.34])
    fig_letter(fig, 0.005, 0.43, "b")
    led = pd.DataFrame([("V4: Visium (3 sources)", 53, "#2f63a3"), ("V4: Xenium / CosMx imaging units", 13, "#2f63a3"),
                        ("R4: Visium sections", 4, "#c4502a"), ("R4: melanoma legacy ST sections", 8, "#c4502a"),
                        ("0906: Visium sections", 61, "#9aa5b1"), ("0906: CosMx ccRCC tissues", 3, "#9aa5b1")],
                       columns=["dataset", "n", "c"])
    export(led[["dataset", "n"]], "figS5_spatial_units.tsv")
    ax2.barh(range(len(led))[::-1], led.n, color=led.c, height=0.65)
    for i, r in zip(range(len(led))[::-1], led.itertuples()):
        ax2.text(r.n + 1, i, str(r.n), va="center", fontsize=5)
    ax2.set_yticks(range(len(led))[::-1], led.dataset, fontsize=5.2); ax2.set_xlabel("sections / imaging units")
    ax2.set_title("Spatial units used per discovery version", fontsize=6.2)
    ax3 = fig.add_axes([0.6, 0.08, 0.36, 0.34])
    fig_letter(fig, 0.53, 0.43, "c")
    d = imgq[(imgq.status == "TESTED") & (imgq.threshold == 1) & (imgq.radius_um == 30) & (imgq.tile_um == 200)].copy()
    d["lab"] = d.section + " " + d.direction.str.replace("_to_", "→").str.replace("_EC", "").str.replace("_mural", "")
    y = np.arange(len(d))[::-1]
    ax3.scatter(d.effect_fraction * 100, y, c=["#c4502a" if q < 0.05 else "#9aa5b1" for q in d.q_family], s=14)
    ax3.axvline(0, color=MUTED, lw=0.5, ls="--")
    ax3.set_yticks(y, d.lab, fontsize=4.6); ax3.set_xlabel("observed − expected neighbour fraction (%)")
    ax3.set_title("R4 breast imaging: ESM1 EC ↔ RGS5 mural (30 µm)\nred: family q < 0.05", fontsize=6.0)
    export(d[["section", "direction", "threshold", "radius_um", "effect_fraction", "p_exact", "q_family"]], "figS5_imaging_core.tsv")
    save(fig, "figS5_spatial_first", ["R4 imaging_core_review/primary_summary.tsv", "unit counts from analysis READMEs"])


# ---------------------------------------------------------------------------
# S6 CM discovery improvements
# ---------------------------------------------------------------------------
def figS6() -> None:
    sc = pd.read_csv(SRC2 / "figS6_immune_atlas_cm_scopes.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 3.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1.05], wspace=0.9)
    ax = fig.add_subplot(gs[0])
    s = sc.iloc[::-1]
    y = np.arange(len(s))
    ax.barh(y, s["rank"], color="#d9dee5", height=0.7, label="NMF rank")
    ax.barh(y, s.supported, color="#7b5ea7", height=0.7, label="CMs with network support")
    for yi, e in zip(y, s.edges):
        ax.text(13.3, yi, str(e), fontsize=4.8, va="center", color=MUTED)
    ax.text(13.3, len(s) - 0.2, "edges", fontsize=4.8, color=MUTED)
    ax.set_yticks(y, s.scope.str.replace("_", " "), fontsize=4.6); ax.set_xlim(0, 15)
    ax.legend(fontsize=4.8, loc="lower right", bbox_to_anchor=(1.0, -0.22), ncol=2)
    ax.set_title("Immune atlas v1: CMs per material scope", fontsize=6.2)
    panel(ax, "a", x=-0.9)

    ax = fig.add_subplot(gs[1])
    comp = pd.DataFrame({"metric": ["cells (M)", "samples (÷100)", "sources (÷10)", "K"],
                         "0906 K9": [4.675, 13.03, 5.4, 9], "V4 K14": [0.561, 1.64, 1.3, 14]})
    export(comp, "figS6_k9_vs_k14.tsv")
    y = np.arange(len(comp))[::-1]
    ax.barh(y + 0.18, comp["0906 K9"], height=0.34, color="#c4502a", label="0906 K9 (all lineages)")
    ax.barh(y - 0.18, comp["V4 K14"], height=0.34, color=BLUE, label="V4 K14 (non-epithelial)")
    for yi, a, b in zip(y, comp["0906 K9"], comp["V4 K14"]):
        if np.isfinite(a):
            ax.text(a + 0.2, yi + 0.18, f"{a:g}", fontsize=4.6, va="center")
        ax.text(b + 0.2, yi - 0.18, f"{b:g}", fontsize=4.6, va="center")
    ax.set_yticks(y, comp.metric, fontsize=5.2)
    ax.legend(fontsize=4.6, loc="lower right", bbox_to_anchor=(1.0, -0.3))
    ax.set_title("Removing epithelium: input change", fontsize=6.2)
    panel(ax, "b", x=-0.55)

    ax = fig.add_subplot(gs[2]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0, 1.0, "Three layers to compare", fontsize=6.2, va="top")
    layers = [("Composition CM,\nnon-epithelial", "cross-cancer comparable\nTME themes", BLUE),
              ("Composition CM +\nmalignant / epithelial", "tumour-state × TME\nco-variation (purity risk)", "#c4502a"),
              ("State layer within\ntrusted identities", "cNMF programmes: mean,\nhigh quantile, active fraction", "#7b5ea7")]
    for i, (h, t, c) in enumerate(layers):
        y0 = 0.68 - i * 0.3
        box(ax, 0, y0 - 0.02, 1, 0.26, "", fc="white", ec=c)
        ax.text(0.05, y0 + 0.17, h, fontsize=5.2, fontweight="bold", color=c, va="center", linespacing=1.1)
        ax.text(0.05, y0 + 0.04, t, fontsize=4.6, color=INK, va="center", linespacing=1.1)
    panel(ax, "c", x=-0.05, y=1.0)
    save(fig, "figS6_cm_improvement", ["figS6_immune_atlas_cm_scopes.tsv", "V4 K14 README numbers"])


if __name__ == "__main__":
    for f in [fig8, figS1, figS2, figS3, figS4, figS5, figS6]:
        f()
