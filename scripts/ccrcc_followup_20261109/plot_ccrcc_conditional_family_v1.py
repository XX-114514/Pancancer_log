"""Show all 96 comparisons, separating conditional excess from raw differences."""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

r = Path(sys.argv[1]).resolve(); o = r / 'followup_v1'
assert json.loads((o / 'independent_audit.json').read_text())['audit_pass']
table_path = o / 'spatial_complete96.tsv'
t = pd.read_csv(table_path, sep='\t')
assert len(t) == 96
destination = o / 'conditional_family96.png'
assert not destination.exists(), 'Preserve previous figures; use a new version for changes'
rows = [(typ, den) for typ in ['B', 'T', 'NK'] for den in ['allQC_nonEC', 'B_T_NK_only']]
cols = [(geo, rad) for geo in ['centroid', 'EC_cell_polygon'] for rad in [25, 50, 100, 200]]
vmax = max(.1, float(t.excess_over_null_log_odds.abs().max()))
fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)
for ax, tile in zip(axes, [250, 500]):
    values = np.zeros((6, 8)); labels = []
    for i, (typ, denominator) in enumerate(rows):
        for j, (geometry, radius) in enumerate(cols):
            z = t[(t.receiver == typ) & (t.denominator == denominator) &
                  (t.geometry == geometry) & (t.radius_um == radius) & (t.tile_um == tile)]
            assert len(z) == 1; z = z.iloc[0]
            values[i, j] = z.excess_over_null_log_odds
            marker = '**' if z.BY_q_complete96 < .05 else ('*' if z.BH_q_complete96 < .05 else '')
            ax.text(j, i, f'{z.BH_q_complete96:.3f}' + marker,
                    ha='center', va='center', fontsize=7, color='black')
    im = ax.imshow(values, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto', zorder=-1)
    ax.set_xticks(range(8), [f'{"Centroid" if g == "centroid" else "EC region"}\n{rad} µm' for g, rad in cols], fontsize=8)
    ax.set_yticks(range(6), [typ + ' | ' + ('all non-EC' if den == 'allQC_nonEC' else 'B/T/NK only') for typ, den in rows], fontsize=9)
    ax.axvline(3.5, color='white', lw=2)
    ax.set_title(f'Local null: {tile} µm tiles × EC RNA-depth strata')
    ax.set_xlabel('Text: BH q across all 96; * BH<0.05; ** also BY<0.05', fontsize=8)
fig.suptitle('ccRCC, one donor: conditional neighborhood excess\nACKR1-independent EC identity; exploratory, not patient-level inference', fontsize=12)
fig.colorbar(im, ax=axes.ravel().tolist(), label='Observed log-odds − local-null mean', fraction=.02, pad=.025)
fig.subplots_adjust(left=.12, right=.9, bottom=.18, top=.84, wspace=.1)
fig.savefig(destination, dpi=180); plt.close(fig)
receipt = {'plot_input_table_sha256': hashlib.sha256(table_path.read_bytes()).hexdigest(),
           'PNG_sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
           'plotted_comparisons': 96, 'BH_lt05': int((t.BH_q_complete96 < .05).sum()),
           'BY_lt05': int((t.BY_q_complete96 < .05).sum()),
           'new_hypothesis_tests': 0, 'spatial_complete96_unchanged': True}
(o / 'family96_figure_audit.json').write_text(json.dumps(receipt, indent=2))
print(json.dumps(receipt))
