"""Render only audited, non-row-level aggregate tables; no reannotation."""
from pathlib import Path
import sys
import os
import json
import hashlib
import numpy as np
import pandas as pd

r = Path(sys.argv[1]).resolve()
os.environ.setdefault('MPLCONFIGDIR', str(r / 'mpl_cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

assert json.loads((r / 'independent_forward_audit.json').read_text())['audit_pass']
out = r / 'figures_v2'; out.mkdir(exist_ok=True)
comp = pd.read_csv(r / 'CM03_existing_identity_spatial12.tsv', sep='\t')
lr = pd.read_csv(r / 'spatial_expression_communication_complete108.tsv', sep='\t')
assert len(comp) == 12 and len(lr) == 108
relations = ['EC_B_edges', 'EC_NK_edges', 'same_EC_B_and_NK']
row_keys = [(relation, radius) for relation in relations for radius in [50, 100]]
values = np.zeros((6, 2)); text = {}
for i, (relation, radius) in enumerate(row_keys):
    for j, tile in enumerate([250, 500]):
        z = comp[(comp.relation == relation) & (comp.radius_um == radius) & (comp.tile_um == tile)].iloc[0]
        values[i, j] = 100 * (z.observed_count / z.null_mean - 1)
        marker = '**' if z.BY_q_complete12 < .05 else ('*' if z.BH_q_complete12 < .05 else '')
        text[i, j] = f'{int(z.observed_count)} / {z.null_mean:.1f}\nBH q={z.BH_q_complete12:.3f}{marker}'
fig, ax = plt.subplots(figsize=(8, 7)); lim = max(5, np.max(abs(values)))
im = ax.imshow(values, cmap='RdBu_r', vmin=-lim, vmax=lim, aspect='auto')
for (i, j), label in text.items():
    ax.text(j, i, label, ha='center', va='center', fontsize=9,
            color='white' if abs(values[i, j]) > .7 * lim else 'black')
ax.set_xticks([0, 1], ['250 µm local tiles', '500 µm local tiles'])
names = {'EC_B_edges': 'EC–B edges', 'EC_NK_edges': 'EC–NK edges', 'same_EC_B_and_NK': 'EC has B AND NK'}
ax.set_yticks(range(6), [names[rel] + f' | {rad} µm' for rel, rad in row_keys])
ax.set_title('ccRCC CM03 spatial tests with existing identities\nNumbers: observed / local-null mean; complete12 BH/BY', fontsize=11)
fig.colorbar(im, ax=ax, label='Count excess over local null (%)', fraction=.035)
fig.text(.5, .015, '* BH<0.05; ** also BY<0.05. One donor, conditional exploratory tests.', ha='center', fontsize=8)
fig.tight_layout(rect=(0, .035, 1, 1)); fig.savefig(out / 'CM03_current_identity_spatial12.png', dpi=180); plt.close(fig)
senders = ['Blood_EC', 'Mural', 'Lymphatic']; receivers = ['B', 'T', 'NK']
pairs = [(s, t) for s in senders for t in receivers]
columns = [(radius, tile) for tile in [250, 500] for radius in [50, 100]]
fig, axes = plt.subplots(1, 3, figsize=(16, 8), sharey=True)
lim = max(1, float(lr.difference_pp.abs().max()))
for ax, (ligand, receptor) in zip(axes, [('CCL19', 'CCR7'), ('CXCL12', 'CXCR4'), ('CCL21', 'CCR7')]):
    values = np.full((9, 4), np.nan); texts = {}
    for i, (sender, receiver) in enumerate(pairs):
        for j, (radius, tile) in enumerate(columns):
            z = lr[(lr.sender == sender) & (lr.receiver == receiver) & (lr.ligand == ligand) &
                   (lr.radius_um == radius) & (lr.tile_um == tile)].iloc[0]
            if z.available:
                values[i, j] = z.difference_pp
                marker = '**' if z.BY_q_complete108 < .05 else ('*' if z.BH_q_complete108 < .05 else '')
                texts[i, j] = f'{z.BH_q_complete108:.3f}{marker}'
            else:
                texts[i, j] = 'NA'
    cmap = plt.get_cmap('RdBu_r').copy(); cmap.set_bad('#dddddd')
    im = ax.imshow(values, cmap=cmap, vmin=-lim, vmax=lim, aspect='auto')
    for (i, j), label in texts.items():
        ax.text(j, i, label, ha='center', va='center', fontsize=7,
                color='white' if np.isfinite(values[i, j]) and abs(values[i, j]) > .7 * lim else 'black')
    ax.set_xticks(range(4), [f'r={rad}\ntile={tile}' for rad, tile in columns], fontsize=8)
    ax.set_yticks(range(9), [s + ' → ' + t for s, t in pairs], fontsize=8)
    ax.set_title(ligand + '–' + receptor + (' (not measured)' if ligand == 'CCL21' else ''))
    ax.axvline(1.5, color='white', lw=2)
    ax.set_xlabel('Radius / local tile (µm); text = complete108 BH q', fontsize=8)
fig.suptitle('ccRCC: local ligand/receptor expression coupling with existing identities\nFixed 6NN graph; 999 expression-profile shuffles; * BH<0.05; ** also BY<0.05', fontsize=12)
fig.subplots_adjust(left=.13, right=.88, top=.88, bottom=.12, wspace=.15)
cax = fig.add_axes([.915, .18, .012, .62])
fig.colorbar(im, cax=cax, label='Observed − local-null dual-positive edge fraction (pp)')
fig.savefig(out / 'current_identity_spatial_LR108.png', dpi=180); plt.close(fig)
files = ['CM03_current_identity_spatial12.png', 'current_identity_spatial_LR108.png']
receipt = {'layout_revision': 2, 'renderer_Matplotlib': matplotlib.__version__, 'input_audit_pass': True, 'new_hypothesis_tests': 0,
           'identity_changed': False,
           'input_tables': {n: hashlib.sha256((r / n).read_bytes()).hexdigest() for n in
                            ['CM03_existing_identity_spatial12.tsv', 'spatial_expression_communication_complete108.tsv']},
           'figures': {n: hashlib.sha256((out / n).read_bytes()).hexdigest() for n in files}}
(out / 'figure_manifest.json').write_text(json.dumps(receipt, indent=2)); print(json.dumps(receipt))
