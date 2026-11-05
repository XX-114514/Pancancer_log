#!/usr/bin/env python3
"""Bounded display-only UMAP from existing batch-benchmark embeddings.

Run in the existing scanpy environment. Never reads expression matrices, changes
integration results, exports cell rows, or fits a batch-correction/CM model.
UMAP is newly fitted on the same fixed 20,000 cells for every displayed method.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
from datetime import datetime, timezone

import h5py
import numpy as np
from anndata.io import read_elem
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from umap import UMAP

METHODS = [
    ('unintegrated', 'Unintegrated PCA'),
    ('harmony_source_unit', 'Harmony: source + unit'),
    ('harmony_source_t2', 'Harmony: source, theta = 2'),
    ('scvi_source', 'scVI: source'),
    ('scvi_unit', 'scVI: unit'),
    ('scanorama', 'Scanorama'),
]


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--benchmark-root', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--n', type=int, default=20000)
    args = ap.parse_args()
    assert 1000 <= args.n <= 20000, 'Display-only size limit'
    args.output.mkdir(parents=True, exist_ok=True)
    assert not (args.output / 'umap_receipt.json').exists(), 'Do not overwrite completed rendering'
    seed = 20261105
    with h5py.File(args.benchmark_root / 'sub.h5ad', 'r') as f:
        obs = read_elem(f['obs'])
    idx = np.sort(np.random.default_rng(seed).choice(len(obs), args.n, replace=False))
    display = obs.iloc[idx]
    sources = sorted(obs.source.astype(str).unique())
    labels = sorted(obs.label_coarse.astype(str).unique())
    source_colors = {s: plt.cm.hsv(i / len(sources)) for i, s in enumerate(sources)}
    palette = ['#0072B2', '#E69F00', '#009E73', '#999999', '#CC79A7', '#56B4E9', '#D55E00', '#332288', '#88CCEE', '#44AA99']
    label_colors = dict(zip(labels, palette))
    label_colors['ambiguous'] = '#b8bec5'
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                         'axes.titlesize': 10, 'pdf.fonttype': 42,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(4, 3, figsize=(12, 12))
    fig.subplots_adjust(left=.035, right=.995, top=.90, bottom=.11, wspace=.08, hspace=.20)
    fig.suptitle('TME v1 batch benchmark | same 20,000 cells in each panel', fontsize=14, y=.985)
    fig.text(.5, .953, 'New display UMAP only; original benchmark metrics use larger cohorts', ha='center', fontsize=9)
    receipt = {'generated_at_utc': datetime.now(timezone.utc).isoformat(),
               'scope': 'TME v1 batch_benchmark_v1; not final TME run_v2',
               'benchmark_cells': len(obs), 'display_cells': args.n, 'seed': seed,
               'selection': 'uniform without replacement, same sorted row indices for all methods',
               'row_index_sha256': hashlib.sha256(idx.astype('<i8').tobytes()).hexdigest(),
               'obs_index_sha256': hashlib.sha256('\n'.join(obs.index.astype(str)).encode()).hexdigest(),
               'source_count': len(sources), 'source_colors': {s: matplotlib.colors.to_hex(c) for s, c in source_colors.items()},
               'display_label_counts': display.label_coarse.astype(str).value_counts().to_dict(),
               'umap': {'n_neighbors': 15, 'min_dist': .5, 'metric': 'euclidean', 'n_epochs': 300, 'random_state': seed, 'n_jobs': 1},
               'interpretation': 'Display-only induced 20k neighborhoods; not full-cohort UMAP, not paired trajectories, marker labels are not independent truth.',
               'packages': {x: importlib.metadata.version(x) for x in ['numpy', 'umap-learn', 'matplotlib', 'anndata', 'h5py']},
               'input_embeddings': {}, 'method_order': [m for m, _ in METHODS]}
    for k, (method, title) in enumerate(METHODS):
        path = args.benchmark_root / ('emb_' + method + '.npy')
        matrix = np.load(path, mmap_mode='r')
        assert matrix.shape[0] == len(obs)
        values = np.asarray(matrix[idx], dtype=np.float32)
        assert np.isfinite(values).all()
        coord = UMAP(**receipt['umap']).fit_transform(values)
        assert coord.shape == (args.n, 2) and np.isfinite(coord).all()
        row, col = (k // 3) * 2, k % 3
        order = np.random.default_rng(seed + 1).permutation(args.n)
        for off, field, cmap in [(0, 'source', source_colors), (1, 'label_coarse', label_colors)]:
            ax = axes[row + off, col]
            colors = [cmap[s] for s in display[field].astype(str).iloc[order]]
            ax.scatter(coord[order, 0], coord[order, 1], c=colors, s=.65, alpha=.75, linewidths=0, rasterized=True)
            ax.set_xticks([]); ax.set_yticks([])
            ax.spines[['top', 'right', 'left', 'bottom']].set_visible(False)
            ax.set_title((f'{chr(97 + k)}  {title}\nSource' if off == 0 else 'Marker-derived broad identity'), pad=4)
        receipt['input_embeddings'][method] = {'logical_path': 'tme_v1/batch_benchmark_v1/' + path.name,
                                              'sha256': digest(path), 'shape': list(matrix.shape)}
        print('UMAP_COMPLETE', method, flush=True)
    handles = [Line2D([0], [0], marker='o', linestyle='', color=label_colors[l], markersize=4, label=l) for l in labels]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .038), ncol=5, frameon=False, fontsize=8)
    fig.text(.5, .018, 'Source colors are fixed across methods; complete source legend supplied separately. Axes are arbitrary UMAP coordinates.', ha='center', fontsize=8)
    for ext in ['png', 'pdf']:
        fig.savefig(args.output / ('fig2_batch_umap.' + ext), dpi=240)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(12, 4.6)); ax.axis('off')
    handles = [Line2D([0], [0], marker='o', linestyle='', color=source_colors[s], markersize=5, label=s) for s in sources]
    ax.legend(handles=handles, loc='center', ncol=6, frameon=False, fontsize=7)
    ax.set_title('Fixed source colors for the TME v1 benchmark display', fontsize=12)
    for ext in ['png', 'pdf']:
        fig.savefig(args.output / ('figS1_source_legend.' + ext), dpi=240, bbox_inches='tight')
    plt.close(fig)
    receipt['outputs_sha256'] = {p.name: digest(p) for p in args.output.glob('*.png')}
    receipt['outputs_sha256'].update({p.name: digest(p) for p in args.output.glob('*.pdf')})
    (args.output / 'umap_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('DISPLAY_RENDER_COMPLETE', flush=True)


if __name__ == '__main__':
    main()
