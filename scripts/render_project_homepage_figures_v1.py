"""Render publication panels from reviewed aggregate statistics, without refitting.

Usage: python render_project_homepage_figures_v1.py aggregates.json new_output_dir
The --stdin-bundle mode renders on the analysis host and returns image bytes only.
"""
import argparse
import base64
import hashlib
import io
import json
from pathlib import Path


def render(data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch
    from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none", "figure.facecolor": "white"})
    navy, teal, muted = "#16324F", "#147D92", "#708090"
    images = {}

    def save(fig, name):
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=180, bbox_inches="tight", facecolor="white")
        images[name] = buf.getvalue()
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1.15, 1]})
    fig.subplots_adjust(top=.78, bottom=.15, wspace=.42)
    fig.suptitle("A source-balanced objective retains the CM02 core", x=.04, ha="left",
                 fontsize=18, weight="bold", color=navy, y=.97)
    fig.text(.04, .85, "2026-10-14 immune atlas v1  |  T_whole_TME_primary_tumor", color=muted)
    ax = axes[0]
    ax.set_axis_off()
    for j, node in enumerate(data["CM02_core_members"]):
        row, col = divmod(j, 2)
        x, y = col * .51, .79 - row * .23
        ax.add_patch(FancyBboxPatch((x, y), .46, .16, boxstyle="round,pad=.015",
                                    transform=ax.transAxes, facecolor="#E8F3F4", edgecolor="none"))
        ax.text(x+.23, y+.08, node, transform=ax.transAxes, ha="center", va="center",
                fontsize=11, color=navy)
    ax.text(0, -.07, "7 / 7 original network members remain in the top 10\nMembership display; no edges or loading magnitudes encoded.",
            transform=ax.transAxes, fontsize=9, color=muted, va="top")
    ax = axes[1]
    rows = [r for r in data["source_balance"] if r["mode"] == "equal_source_weight"]
    for j, r in enumerate(rows):
        lo, hi = r["A_cosine_range"]
        color = teal if r["module"] == "CM02" else "#A5B1BC"
        ax.hlines(j, lo, hi, color=color, lw=5)
        ax.plot((lo+hi)/2, j, "o", color=color, ms=9)
        ax.text(1.008, j, f"{lo:.3f}", va="center", fontsize=10, color=color)
    ax.set_yticks(range(len(rows)), [r["module"] for r in rows])
    ax.set_ylim(len(rows)-.5, -.5)
    ax.set_xlim(.83, 1.065)
    ax.set_xticks([.85, .90, .95, 1.00])
    ax.set_xlabel("Effective-basis cosine similarity to the original fit")
    ax.set_title("All four T-primary-tumor components", loc="left", fontsize=11)
    ax.grid(axis="x", alpha=.15)
    fig.text(.04, .015, "Ranges span four seeds on the same input. Cosine similarity measures component geometry.",
             fontsize=9, color=muted)
    save(fig, "cm02_source_balance.png")

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 3.8))
    fig.subplots_adjust(top=.72, bottom=.28, left=.06, right=.98, wspace=.45)
    fig.suptitle("A fixed regulatory–cytotoxic RNA proxy has clinical associations",
                 fontsize=16, weight="bold", color=navy, x=.04, ha="left", y=.97)
    for ax, row in zip(axes, data["clinical_CM02_main"]):
        v, lo, hi = row["effect"], row["CI_low"], row["CI_high"]
        ax.hlines(0, lo, hi, color=teal, lw=3)
        ax.plot(v, 0, "o", color=teal, ms=10)
        ax.axvline(1, color=muted, ls="--", lw=1)
        ax.set_xscale("log")
        is_response = row["endpoint"] == "response"
        ax.set_xlim((.75, 9) if is_response else (.22, 1.3))
        ax.xaxis.set_major_locator(FixedLocator([1, 2, 4, 8] if is_response else [.25, .5, 1]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_ylim(-.8, .8)
        ax.set_yticks([])
        ax.set_title(f"{row['cohort']}  |  {'ICI response' if is_response else 'Overall survival'}", loc="left", fontsize=12)
        metric = "OR" if is_response else "HR"
        ax.set_xlabel(f"{metric} per SD of the fixed joint score (95% CI)")
        ax.text(.02, .05, f"{metric} {v:.2f} [{lo:.2f}, {hi:.2f}]\nn={row['n']}; events={row['events']}; BH q={row['BH_q']:.4g}",
                transform=ax.transAxes, fontsize=10, color=navy)
    fig.text(.04, .08, "Joint score = min(regulatory-associated panel, cytotoxic panel); primary models include pan-T background.\nTwo distinct cohorts and endpoints; each retains its original test family. This panel presents association estimates.",
             fontsize=9, color=muted)
    save(fig, "cm02_clinical_associations.png")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), sharey=True)
    fig.subplots_adjust(top=.72, bottom=.23, left=.17, right=.98, wspace=.25)
    fig.suptitle("NSCLC regional lymph nodes organize relative T-cell protein states",
                 x=.025, ha="left", fontsize=16, weight="bold", color=navy, y=.97)
    fig.text(.025, .83, "Anchor: CD4-like FOXP3/CTLA4 relative-high candidate  |  11 registered patient units, 97 images", color=muted, fontsize=10)
    labels = ["PD1 / TIM3", "MKI67", "CCR7 / CD127"]
    rows = data["NSCLC_state_pairs"]
    for ax, field, title, xlabel in zip(axes,
            ["within_MKI67_MH_OR", "local_centered_enrichment_factor"],
            ["Association within MKI67 strata", "Enrichment beyond the local null"],
            ["Spatial odds ratio (95% bootstrap CI)", "Centered enrichment factor (95% bootstrap CI)"]):
        for j, row in enumerate(rows):
            z = row[field]
            ax.hlines(j, z["lower"], z["upper"], lw=2.5, color=teal)
            ax.plot(z["point"], j, "o", ms=8, color=teal)
            ax.annotate(f"{z['point']:.2f}", (z["point"], j), xytext=(0, 11), textcoords="offset points", ha="center", color=navy)
        ax.axvline(1, color=muted, ls="--", lw=1)
        ax.set_xscale("log")
        ax.set_xlim((.85, 4.1) if field == "within_MKI67_MH_OR" else (.96, 1.68))
        ax.xaxis.set_major_locator(FixedLocator([1, 2, 3, 4] if field == "within_MKI67_MH_OR" else [1, 1.2, 1.4, 1.6]))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda value, pos: f"{value:g}"))
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_title(title, loc="left", fontsize=12)
        ax.set_xlabel(xlabel)
        ax.set_yticks(range(3), labels)
        ax.set_ylim(2.65, -.75)
        ax.grid(axis="x", alpha=.15)
    fig.text(.025, .04, "All three prespecified endpoints shown. Main BH q=0.00209; local BH q=0.003; BY q<0.01.\nRelative-protein state comparisons; the two estimands retain their own definitions.",
             fontsize=9, color=muted)
    save(fig, "nsclc_state_organization.png")
    return images


def main():
    import sys
    if sys.argv[1:] == ["--stdin-bundle"]:
        data = json.load(sys.stdin)
        images = render(data)
        print(json.dumps({"kind": "RENDERED_ACCEPTED_AGGREGATES_NO_REFIT", "files": {
            n: {"SHA256": hashlib.sha256(b).hexdigest(), "base64": base64.b64encode(b).decode("ascii")}
            for n, b in images.items()}}))
        return
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    images = render(json.loads(args.input.read_text()))
    args.output.mkdir(parents=True, exist_ok=False)
    for name, content in images.items():
        with (args.output / name).open("xb") as handle:
            handle.write(content)
    print(json.dumps({"rendered": list(images), "new_inference": False}))


if __name__ == "__main__":
    main()
