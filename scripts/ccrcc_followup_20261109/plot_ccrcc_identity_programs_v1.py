"""Compact multi-marker and protein review; diagnostic groups can overlap."""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ccrcc_followup_v1 as base

r = Path(sys.argv[1]).resolve(); o = r / 'identity_refinement_v2'
assert json.loads((o / 'independent_audit.json').read_text())['audit_pass']
assert not (o / 'identity_programs_marker.png').exists()
ids, c, xy, xr, rn, xp, pn, scale, total, ng, a = base.load(r)
panels = {
    'B': ['CD19', 'MS4A1', 'CD79A', 'BANK1'],
    'ASC': ['MZB1', 'TENT5C', 'SDC1', 'SLAMF7', 'PRDM1', 'CD38'],
    'pDC_review': ['IL3RA', 'LILRA4', 'TCF4', 'IRF8', 'PLD4', 'GZMB'],
    'T': ['CD3D', 'CD3E', 'TRAC', 'CD2'],
    'NK_review': ['KLRD1', 'GNLY', 'NKG7', 'KLRC1', 'KLRF1', 'NCAM1', 'PRF1'],
    'EC': ['PECAM1', 'VWF', 'EMCN', 'ERG', 'CD34', 'ADGRL4', 'CD93', 'CLEC14A'],
    'Lymphatic': ['PROX1', 'LYVE1', 'PDPN', 'MMRN1'],
    'Mural': ['PDGFRB', 'ACTA2', 'MYH11'],
    'Fibroblast': ['PDGFRA', 'DCN', 'LUM', 'COL1A1'],
    'Renal_epithelial': ['EPCAM', 'CDH16', 'LRP2', 'CA9'],
    'Podocyte_review': ['NPHS1', 'WT1', 'SYNPO'],
    'Exposure_only': ['ACKR1']}
base.PANELS = panels
for g in sorted(set(sum(panels.values(), []))):
    if g in rn and g not in a:
        a[g] = xr[:, list(rn).index(g)].toarray().ravel()
spec = [('identity_refinement_v2/B_ASC/RNA_only_refinement_PRIVATE.h5ad', [('0', 'B program'), ('1', 'ASC program / discordance'), ('2', 'pDC-like program')]),
        ('identity_refinement_v2/T_NK/RNA_only_refinement_PRIVATE.h5ad', [('0', 'T IL7R program'), ('1', 'T CD8A program'), ('2', 'NK-like program')]),
        ('followup_v1/clusters_vascular_stromal/RNA_only_cluster_object_PRIVATE.h5ad',
         [('0', 'EC program'), ('9', 'Lymphatic program'), ('3', 'Mural ACTA2/MYH11'), ('6', 'Mural PDGFRB'),
          ('4', 'Epithelial CDH16'), ('7', 'Renal epithelial LRP2'), ('10', 'Podocyte-like program')])]
markers, proteins, index_rows = [], [], []
for path, groups in spec:
    obj = ad.read_h5ad(r / path, backed='r')
    index = obj.obs.source_index_PRIVATE.to_numpy(int); labels = obj.obs.seed_0.astype(str).to_numpy()
    for cluster, description in groups:
        ix = index[labels == cluster]; name = description + ' | C' + cluster
        group = np.repeat(name, len(ix))
        markers.append(base.marker_table({g: v[ix] for g, v in a.items()}, group, [name]))
        proteins.append(base.protein_table(xp[ix], pn, group, [name]))
        index_rows.append({'review_group': name, 'cells': len(ix), 'parent_object': path,
                           'seed0_cluster': cluster, 'purpose': 'diagnostic_review_not_frozen_label'})
    obj.file.close()
mt = pd.concat(markers, ignore_index=True); pt = pd.concat(proteins, ignore_index=True)
mt.to_csv(o / 'identity_programs_marker.tsv', sep='\t', index=False, na_rep='NA')
pt.to_csv(o / 'identity_programs_protein.tsv', sep='\t', index=False, na_rep='NA')
pd.DataFrame(index_rows).to_csv(o / 'identity_programs_index.tsv', sep='\t', index=False)
base.marker_plot(mt, o / 'identity_programs_marker.png', 'ccRCC identity programs: markers + expected identity (pools can overlap; not frozen labels)')
keep = ['CD20', 'CD3E', 'CD31', 'CD68', 'GranzymeB', 'alphaSMA', 'PanCK', 'CD45']
group_order = [row['review_group'] for row in index_rows]
mat = pt[pt.antibody.isin(keep)].pivot(index='group', columns='antibody', values='median_mean_intensity').loc[group_order, keep]
values = mat.to_numpy(); log = np.log1p(values)
colour = (log - log.mean(0, keepdims=True)) / np.maximum(log.std(0, keepdims=True), 1e-8)
fig, ax = plt.subplots(figsize=(12, 9))
im = ax.imshow(colour, aspect='auto', cmap='coolwarm', vmin=-2, vmax=2)
for i in range(values.shape[0]):
    for j in range(values.shape[1]):
        ax.text(j, i, f'{values[i, j]:.1f}', ha='center', va='center', fontsize=8)
ax.set_yticks(range(len(mat)), mat.index, fontsize=8)
ax.set_xticks(range(len(keep)), [g + ' | ' + {'GranzymeB': 'cytotoxic/pDC review'}.get(g, base.PROTEIN_EXPECTED.get(g, 'review')) for g in keep], rotation=55, ha='right', fontsize=8)
ax.set_title('Protein cross-review of RNA programs\nNumbers: group median mean intensity / stored scaling factor; compare within each antibody only', fontsize=10)
fig.colorbar(im, ax=ax, label='Colour standardized within antibody', fraction=.025)
fig.tight_layout(); fig.savefig(o / 'identity_programs_protein.png', dpi=180); plt.close(fig)
receipt = {'diagnostic_groups': len(index_rows), 'pools_can_overlap': True,
           'source_protein_scaling_factor': scale, 'new_statistical_tests': 0,
           'label_freeze_performed': False,
           'files': {n: hashlib.sha256((o / n).read_bytes()).hexdigest() for n in
                     ['identity_programs_marker.tsv', 'identity_programs_protein.tsv', 'identity_programs_index.tsv',
                      'identity_programs_marker.png', 'identity_programs_protein.png']}}
(o / 'identity_programs_figure_receipt.json').write_text(json.dumps(receipt, indent=2)); print(json.dumps(receipt))
