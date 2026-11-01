# Marker-gene overlays on the integration v2 UMAP (log1p-normalised decontX expression), to read the
# compartment structure without any annotation. Cell order = integration order (shards sorted by
# source, as in integrate_v2.py); verified against the stored UMAP length.
# Usage: plot_umap_markers_v2.py RUN_DIR
import sys, os, glob
import numpy as np, anndata as ad, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

run = sys.argv[1]
MARKERS = ['PTPRC', 'CD3E', 'NKG7', 'MS4A1', 'JCHAIN', 'CD68', 'LYZ', 'TPSAB1', 'EPCAM', 'KRT19', 'COL1A1', 'ACTA2', 'PECAM1', 'MKI67', 'PMEL', 'GFAP']
U = np.load(f'{run}/report/X_umap.npy')
vals = {m: [] for m in MARKERS}
for f in sorted(glob.glob(f'{run}/shards/*.h5ad')):
    a = ad.read_h5ad(f, backed='r')
    idx = [a.var_names.get_loc(m) if m in a.var_names else -1 for m in MARKERS]
    X = a.X[:, [i for i in idx if i >= 0]]
    X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
    k = 0
    for m, i in zip(MARKERS, idx):
        if i >= 0:
            vals[m].append(X[:, k]); k += 1
        else:
            vals[m].append(np.full(a.n_obs, np.nan, np.float32))
    a.file.close()
vals = {m: np.concatenate(v) for m, v in vals.items()}
assert len(vals['PTPRC']) == len(U), (len(vals['PTPRC']), len(U))
order = np.random.default_rng(0).permutation(len(U))
fig, axes = plt.subplots(4, 4, figsize=(20, 20))
for ax, m in zip(axes.ravel(), MARKERS):
    v = vals[m][order]
    ax.scatter(U[order, 0], U[order, 1], c=np.nan_to_num(v), s=0.02, cmap='viridis', vmin=0, vmax=np.nanpercentile(v, 99.5) or 1, rasterized=True, linewidths=0)
    ax.set_title(m, fontsize=14); ax.set_xticks([]); ax.set_yticks([])
fig.suptitle(f'Integration v2 UMAP ({len(U):,} cells): marker expression (log1p, decontX)', fontsize=16)
fig.tight_layout(); fig.savefig(f'{run}/figures/umap_markers.png', dpi=110)
print('saved', f'{run}/figures/umap_markers.png')
