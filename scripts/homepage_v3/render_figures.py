"""Render all 2026-11-26 homepage figures (main 1-8, supplementary S1-S6).

Usage (from repository root):  python scripts/homepage_v3/render_figures.py
Reads only aggregate tables in this repository; writes PNG + PDF to assets/publication/20261126/
and the plotted aggregate values to inventories/publication/2026-11-26/homepage_v3/.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import style  # noqa: E402
import fig_main_a, fig_main_b, fig_follow  # noqa: E402

for fn in [fig_main_a.fig1, fig_main_a.fig2, fig_main_a.fig3, fig_main_b.fig4, fig_main_b.fig5, fig_main_b.fig6,
           fig_main_b.fig7, fig_follow.fig8, fig_follow.figS1, fig_follow.figS2, fig_follow.figS3, fig_follow.figS4,
           fig_follow.figS5, fig_follow.figS6]:
    fn()
style.write_manifest()
print(f"{len(style.MANIFEST)} figures written")
