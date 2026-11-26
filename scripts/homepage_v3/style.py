"""Shared style for the 2026-11-26 homepage (v3) figures.

Only aggregate tables already present in this repository are read. UMAP panels are
re-used from the 2026-11-02 rendering (cell coordinates stay on the research server).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Rectangle

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "assets" / "publication" / "20261126"
SRC2 = REPO / "inventories" / "publication" / "2026-11-02" / "homepage_v2"
SRC3 = REPO / "inventories" / "publication" / "2026-11-26" / "homepage_v3"
OUT.mkdir(parents=True, exist_ok=True)
SRC3.mkdir(parents=True, exist_ok=True)

for f in font_manager.findSystemFonts():
    if "Inter" in f and ("Inter-" in Path(f).name):
        try:
            font_manager.fontManager.addfont(f)
        except Exception:
            pass

FONT = "Inter" if any(f.name == "Inter" for f in font_manager.fontManager.ttflist) else "DejaVu Sans"

plt.rcParams.update({
    "font.family": FONT,
    "font.size": 6.5,
    "axes.titlesize": 7,
    "axes.labelsize": 6.5,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 5.8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.titlelocation": "left",
    "axes.titleweight": "medium",
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
    "figure.dpi": 150,
    "savefig.dpi": 320,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
})

INK = "#1f2933"
MUTED = "#6b7785"
GRID = "#e6e9ee"
CM_COL = {"CM01": "#1b7f5f", "CM07": "#c4502a", "CM08": "#2f63a3"}
OTHER_CM = "#8a94a3"
LINEAGE_COL = {
    "Tcell": "#3b6ea8", "NK": "#6a9fd4", "ILC": "#a3c4e8", "Bcell": "#7b5ea7", "Plasma": "#b39ddb",
    "Macrophage": "#c2532b", "Monocyte": "#e08a5c", "Dendritic": "#f0b48f", "Mast": "#d4a017",
    "Granulocyte": "#b8860b", "Platelet": "#8d6e63", "Erythroid": "#a1887f", "HSPC": "#795548",
    "Fibroblast": "#2a7f62", "Endothelial": "#5fa58a", "Epithelial": "#d16d8e", "Neural": "#9e7b9b",
    "Melanocytic": "#5d4037", "Unknown": "#c9ced6", "Ambiguous": "#a7afba",
}
DIV = LinearSegmentedColormap.from_list("div", ["#2b5c8a", "#9cbad6", "#f4f4f2", "#e3a38c", "#b23a24"])
SEQ = LinearSegmentedColormap.from_list("seq", ["#f6f8fa", "#bcd3e6", "#5f8fbf", "#1d3f66"])
SEQ_WARM = LinearSegmentedColormap.from_list("seqw", ["#fbf6f2", "#f0c4a8", "#d0784e", "#8c2f18"])

MANIFEST: list[dict] = []


def cmc(cm: str) -> str:
    return CM_COL.get(cm, OTHER_CM)


def panel(ax, letter: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left", color=INK)


def fig_letter(fig, x, y, letter):
    fig.text(x, y, letter, fontsize=9, fontweight="bold", va="bottom", ha="left", color=INK)


def box(ax, x, y, w, h, text, fc="#f3f5f8", ec="#c5ccd6", fs=6.0, color=INK, weight="normal", r=0.012, lw=0.6, ha="center"):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.0,rounding_size={r}", fc=fc, ec=ec, lw=lw,
                       transform=ax.transAxes, clip_on=False)
    ax.add_patch(p)
    tx = x + w / 2 if ha == "center" else x + 0.008
    ax.text(tx, y + h / 2, text, ha=ha, va="center", fontsize=fs, color=color, fontweight=weight,
            transform=ax.transAxes, linespacing=1.25)
    return p


def arrow(ax, x0, y0, x1, y1, color=MUTED, lw=0.8, style="-|>", ms=6, ls="-"):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle=style, mutation_scale=ms, color=color, lw=lw,
                        transform=ax.transAxes, clip_on=False, linestyle=ls)
    ax.add_patch(a)


def export(df, name: str) -> None:
    df = df.copy()
    for c in df.columns:
        if df[c].dtype == object or str(df[c].dtype).startswith("str"):
            df[c] = df[c].astype(str).str.replace("\n", " ", regex=False)
    df.to_csv(SRC3 / name, sep="\t", index=False)


def save(fig, name: str, sources: list[str]) -> None:
    png = OUT / f"{name}.png"
    pdf = OUT / f"{name}.pdf"
    fig.savefig(png)
    fig.savefig(pdf)
    plt.close(fig)
    MANIFEST.append({
        "figure": name,
        "png_sha256": hashlib.sha256(png.read_bytes()).hexdigest(),
        "pdf_sha256": hashlib.sha256(pdf.read_bytes()).hexdigest(),
        "sources": sources,
    })


def write_manifest() -> None:
    (OUT / "figure_manifest.json").write_text(json.dumps({
        "renderer": "scripts/homepage_v3/render_figures.py",
        "font": FONT,
        "note": "Aggregate tables only; UMAP raster panels re-used from assets/publication/20261102 (cell coordinates private).",
        "figures": MANIFEST,
    }, indent=2, ensure_ascii=False))
