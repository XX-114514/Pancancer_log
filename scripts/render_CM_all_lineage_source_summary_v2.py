"""Source-only fixed-panel RNA attribution and both-threshold capture diagnostics."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    root = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    terminal = json.loads((root / 'private_terminal.json').read_text())
    for name in ['private_all_source_scope_major_panel_descriptions.json', 'private_all_source_scope_major_two_arm_overlaps.json']:
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == terminal['files'][name]['SHA256']
    d = pd.DataFrame(json.loads((root / 'private_all_source_scope_major_panel_descriptions.json').read_text()))
    j = pd.DataFrame(json.loads((root / 'private_all_source_scope_major_two_arm_overlaps.json').read_text()))
    scope = 'original_CM_library_material_membership'
    d = d[d.scope.eq(scope)]
    sources = sorted(d.loc[d.source_scope_cells.gt(0), 'source'].unique())
    panels = list(json.loads((root / 'private_frozen_plan.json').read_text()).get('display_panel_order', ['regulatory_associated', 'cytotoxic', 'memory_associated', 'GZMK_effector', 'pan_T', 'interferon', 'proliferation']))
    out = root / 'private_source_figures_v2'
    assert not out.exists(); out.mkdir(mode=0o700)
    fig = plt.figure(figsize=(18, 11))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.55, 1], hspace=.45, wspace=.35)
    ax = fig.add_subplot(grid[:, 0]); top = fig.add_subplot(grid[0, 1]); bottom = fig.add_subplot(grid[1, 1])
    colors = ['#336c99', '#df9e38', '#b8b8b8']
    labels = []; shares = []
    for src in sources:
        for panel in panels:
            g = d[d.source.eq(src) & d.panel.eq(panel)]
            labels.append(src + ' / ' + panel)
            if g.panel_status.eq('DESCRIPTIVE_RNA_CAPTURE_ONLY').all():
                v = g.set_index('major').CP10K_cell_weighted_share
                vals = [v.get('T', 0), v.get('NK', 0), v.loc[~v.index.isin(['T', 'NK'])].sum()]
                assert abs(sum(vals) - 1) < 1e-10
                shares.append(vals)
            else:
                shares.append([np.nan, np.nan, np.nan])
    values = np.asarray(shares); y = np.arange(len(labels)); left = np.zeros(len(labels))
    for k, label in enumerate(['T', 'NK', 'All other frozen majors']):
        ax.barh(y, values[:, k], left=left, color=colors[k], label=label, height=.8)
        left += np.nan_to_num(values[:, k])
    for i, vals in enumerate(values):
        if not np.isfinite(vals).all(): ax.text(.02, i, 'Unavailable complete panel', va='center', fontsize=8)
    for i in range(1, len(sources)): ax.axhline(i * len(panels) - .5, color='0.45', lw=.6)
    ax.set_yticks(y, labels, fontsize=8); ax.invert_yaxis(); ax.set_xlim(0, 1)
    ax.set_xlabel('Equal-gene average of per-gene CP10K source shares')
    ax.set_title('All seven fixed panels: captured-cell-weighted attribution', fontsize=11)
    ax.legend(loc='lower center', bbox_to_anchor=(.5, 1.01), ncol=3, fontsize=9, frameon=False)
    for axis, major in [(top, 'T'), (bottom, 'NK')]:
        for ti, threshold in enumerate([1, 2]):
            vals = []
            for src in sources:
                a = d[d.source.eq(src) & d.major.eq(major) & d.panel.eq('regulatory_associated')].iloc[0]
                b = j[j.source.eq(src) & j.scope.eq(scope) & j.major.eq(major) & j.panelA.eq('regulatory_associated') & j.panelB.eq('cytotoxic') & j.threshold.eq(threshold)].iloc[0]
                n = a['captured_cells_threshold' + str(threshold)]
                vals.append(b.both_captured_cells / n if n and b.both_panels_measurable else np.nan)
            axis.plot(vals, np.arange(len(sources)) + (ti - .5) * .18, ['o', 's'][ti], color=colors[ti], label='Gene count >= ' + str(threshold))
        axis.set_yticks(np.arange(len(sources)), sources, fontsize=10); axis.invert_yaxis(); axis.set_xlim(0, 1.04)
        axis.grid(axis='x', color='0.9'); axis.set_title(major + ': cytotoxic co-capture within regulatory-panel capture', fontsize=10)
        axis.set_xlabel('Fraction; each panel requires >= 2 distinct genes')
        axis.legend(loc='upper right' if major == 'T' else 'lower left', fontsize=9, frameon=False)
    for axis in [ax, top, bottom]: axis.spines[['top', 'right']].set_visible(False)
    fig.subplots_adjust(left=.23, right=.98, top=.89, bottom=.13)
    fig.suptitle('2026-10-14 immune atlas CM (v1): fixed-panel RNA sources\nOriginal CM material libraries; every source with nonempty membership retained', fontsize=14)
    fig.text(.03, .035, 'RNA attribution describes captured cells and sorting context, not tissue cell fractions or decomposition of bulk percentile/min scores.\nRight panels describe same-cell RNA capture; no formal Treg/CD8, receptor engagement, clinical test or direct cell-contact claim.', fontsize=10)
    for suffix in ['png', 'svg']: fig.savefig(out / ('private_original_CM_material_all_panels_and_overlap.' + suffix), dpi=180)
    plt.close(fig)
    files = {p.name:dict(bytes=p.stat().st_size, SHA256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in out.iterdir()}
    assert sum(x['bytes'] for x in files.values()) < 4 * 1024**2
    result = dict(status='SOURCE_SUMMARY_RENDERED_ACTUAL_VISUAL_REVIEW_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), parent_terminal_SHA256=hashlib.sha256((root / 'private_terminal.json').read_bytes()).hexdigest(), files=files)
    (out / 'private_figure_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], files=files)))


if __name__ == '__main__':
    main()
