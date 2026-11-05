"""Independent CellTypist candidate labels; never modify the integration shards."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('NUMBA_NUM_THREADS', '2')
os.environ.setdefault('NUMBA_CACHE_DIR', f'/tmp/celltypist72-numba-{os.getuid()}')
os.environ.setdefault('MPLCONFIGDIR', f'/tmp/celltypist72-matplotlib-{os.getuid()}')
os.makedirs(os.environ['NUMBA_CACHE_DIR'], exist_ok=True)
os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)

import anndata as ad
import celltypist
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--shard', type=Path, required=True)
    p.add_argument('--models', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--repo', type=Path, required=True)
    p.add_argument('--max-cells', type=int, default=0,
                   help='0 means all cells; positive values create a reproducible pilot')
    args = p.parse_args()
    if args.max_cells < 0:
        p.error('--max-cells must be nonnegative')
    args.out.mkdir(parents=True, exist_ok=True)
    # Version-pinned repository wrapper preserves official inference, changing
    # only gene-name array representation to avoid slow object comparisons.
    sys.path.insert(0, str(args.repo / 'scripts'))
    from CellTypist_official_unicode_axis_v1 import official_annotate_unicode_axis
    a = ad.read_h5ad(args.shard, backed='r')
    total = a.n_obs
    idx = np.arange(total)
    if args.max_cells and total > args.max_cells:
        idx = np.sort(np.random.default_rng(20261103).choice(
            total, args.max_cells, replace=False))
    # Materialize X only, leaving raw-count layers and metadata on disk.
    b = ad.AnnData(X=a.X[idx, :].copy(),
                   obs=pd.DataFrame(index=a.obs_names[idx].copy()),
                   var=pd.DataFrame(index=a.var_names.copy()))
    a.file.close()
    if not b.obs_names.is_unique or not b.var_names.is_unique:
        raise ValueError('Duplicate cell or gene IDs')
    values = b.X.data if hasattr(b.X, 'indptr') else b.X
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Expected finite nonnegative log-normalized X')
    check = b.X[:min(100, b.n_obs)].copy()
    if hasattr(check, 'indptr'):
        check.data = np.expm1(check.data)
    else:
        check = np.expm1(check)
    library_sizes = np.asarray(check.sum(axis=1)).ravel()
    if not np.allclose(library_sizes, 10000, rtol=.01):
        raise ValueError('Input does not match log1p(normalize_total 1e4)')
    result = pd.DataFrame(index=b.obs_names)
    result.index.name = 'cell_id'
    summary = {'source': args.shard.stem, 'source_cells': total,
               'predicted_cells': b.n_obs, 'pilot': bool(args.max_cells),
               'celltypist_version': celltypist.__version__,
               'majority_voting': False, 'input': 'v2 X: log1p decontX, total 1e4',
               'interpretation': 'candidate model labels only; no accepted identity',
               'models': {}}
    for name in ('Immune_All_High', 'Immune_All_Low'):
        path = args.models / (name + '.pkl')
        model = celltypist.models.Model.load(str(path))
        overlap = int(b.var_names.isin(model.classifier.features).sum())
        pred = official_annotate_unicode_axis(b, model)
        labels = pred.predicted_labels['predicted_labels'].reindex(result.index)
        probabilities = pred.probability_matrix.reindex(result.index).to_numpy()
        top = np.sort(probabilities, axis=1)[:, -2:]
        result[name + '_candidate'] = labels
        result[name + '_max_probability'] = top[:, -1]
        result[name + '_probability_margin'] = top[:, -1] - top[:, -2]
        summary['models'][name] = {
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'features_shared': overlap,
            'features_total': len(model.classifier.features),
            'candidate_counts': labels.value_counts().to_dict()}
        print(name, 'finished', b.n_obs, 'cells', flush=True)
    # Refuse accidental replacement, including rerunning a pilot as a full source.
    with (args.out / 'candidate_labels.tsv.gz').open('xb') as f:
        result.to_csv(f, sep='\t', compression='gzip')
    with (args.out / 'summary.json').open('x') as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
