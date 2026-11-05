"""Marker-gated reference proposals; all original identities remain unchanged."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ.setdefault(key, '2')
for key, name in [('NUMBA_CACHE_DIR', 'numba'), ('MPLCONFIGDIR', 'matplotlib')]:
    os.environ.setdefault(key, f'/tmp/marker-routed-{name}-{os.getuid()}')
    os.makedirs(os.environ[key], exist_ok=True)

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite

sys.path.insert(0, str(Path(__file__).resolve().parent))
from marker_policy import (PANELS, IMMUNE_HIGH, EXPLICIT_NONIMMUNE,
                           MIN_SCORE, MIN_GAP, MIN_MARKERS,
                           routes, score_clusters, vascular_family, counterexamples)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for flag in ('shard', 'malignancy-dir', 'models', 'repo', 'rscript', 'r-cache', 'out'):
        p.add_argument('--' + flag, type=Path, required=True)
    p.add_argument('--max-cells', type=int, default=0)
    p.add_argument('--seed', type=int, default=20261103)
    args = p.parse_args()
    if args.max_cells < 0:
        p.error('--max-cells must be nonnegative')
    counterexamples()
    args.out.mkdir(parents=True, exist_ok=False)
    mpath = args.malignancy_dir / (args.shard.stem + '.tsv.gz')
    m = pd.read_csv(mpath, sep='\t', usecols=[
        'cell_id', 'unit_id', 'quick_cluster', 'lineage', 'malignancy_call'],
        dtype={'quick_cluster': str}).set_index('cell_id')
    a = ad.read_h5ad(args.shard, backed='r')
    if not m.index.is_unique or set(m.index) != set(a.obs_names):
        raise ValueError('Malignancy/source cell-ID mismatch')
    m = m.reindex(a.obs_names)
    if not np.array_equal(m.unit_id.astype(str), a.obs.unit_id.astype(str)):
        raise ValueError('Malignancy/source unit IDs differ')
    # JSON tuples avoid collisions between source-specific unit/cluster strings.
    groups = pd.Series([json.dumps([str(u), str(c)]) for u, c in
                        zip(m.unit_id, m.quick_cluster)], index=a.obs_names)
    genes = sorted({g for gs in PANELS.values() for g in gs if g in a.var_names})
    marker_x = a.X[:, a.var_names.get_indexer(genes)]
    markers = pd.DataFrame(marker_x.toarray() if sp.issparse(marker_x) else marker_x,
                           index=a.obs_names, columns=genes)
    clusters = score_clusters(markers, groups)
    clusters.to_csv(args.out / 'marker_cluster_evidence.tsv', sep='\t')
    cells = clusters.reindex(groups.to_numpy()).copy()
    cells.index = a.obs_names
    cells = cells[['marker_parent', 'marker_panel', 'marker_score', 'marker_gap',
                   'supporting_markers', 'marker_gate_status']]
    cells['marker_cluster_key'] = groups
    cells['malignancy_call'] = m.malignancy_call
    cells['malignancy_lineage'] = m.lineage
    cells['prediction_route'] = routes(cells.marker_parent, cells.malignancy_call)
    idx = np.arange(a.n_obs)
    if args.max_cells and len(idx) > args.max_cells:
        idx = np.sort(np.random.default_rng(args.seed).choice(
            len(idx), args.max_cells, replace=False))
    cells = cells.iloc[idx].copy()
    cells.index.name = 'cell_id'
    cells['reference_status'] = cells.prediction_route.map({
        'immune_celltypist': 'pending', 'vascular_celltypist': 'pending',
        'singler_hpca': 'pending', 'excluded_malignant': 'excluded_malignant',
        'unresolved_no_compatible_reference': 'unresolved_no_compatible_reference'})
    summary = {'source': args.shard.stem, 'source_cells': a.n_obs,
               'output_cells': len(cells), 'pilot': bool(args.max_cells),
               'marker_scope': 'all source cells, unit-specific quick clusters',
               'marker_thresholds': {'min_mean_logexpr': MIN_SCORE,
                   'min_parent_gap': MIN_GAP, 'min_markers_mean_ge_0.5': MIN_MARKERS},
               'marker_panels': PANELS, 'malignancy_sha256': sha(mpath),
               'script_sha256': sha(__file__),
               'policy_sha256': sha(Path(__file__).with_name('marker_policy.py')),
               'models': {}, 'accepted_fine_identity': False}
    sys.path.insert(0, str(args.repo / 'scripts'))
    import celltypist
    from CellTypist_official_unicode_axis_v1 import official_annotate_unicode_axis
    summary['celltypist_version'] = celltypist.__version__
    for route, names in [('immune_celltypist', ['Immune_All_High', 'Immune_All_Low']),
                         ('vascular_celltypist', ['Adult_Human_Vascular'])]:
        ids = cells.index[cells.prediction_route == route]
        if len(ids) == 0:
            continue
        pos = a.obs_names.get_indexer(ids)
        b = ad.AnnData(X=sp.csr_matrix(a.X[pos, :]),
                      obs=pd.DataFrame(index=ids), var=pd.DataFrame(index=a.var_names))
        values = b.X.data
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError('Invalid normalized expression')
        check = b.X[:min(100, len(ids))].copy(); check.data = np.expm1(check.data)
        if not np.allclose(np.asarray(check.sum(1)).ravel(), 10000, rtol=.01):
            raise ValueError('Expected log1p expression normalized to 10,000')
        for name in names:
            path = args.models / (name + '.pkl')
            model = celltypist.models.Model.load(str(path))
            pred = official_annotate_unicode_axis(b, model)
            lab = pred.predicted_labels.predicted_labels.reindex(ids).astype(str)
            prob = pred.probability_matrix.reindex(ids).to_numpy()
            top = np.sort(prob, axis=1)[:, -2:]
            cells.loc[ids, name + '_candidate'] = lab
            cells.loc[ids, name + '_max_probability'] = top[:, -1]
            cells.loc[ids, name + '_margin'] = top[:, -1] - top[:, -2]
            if name == 'Immune_All_High':
                ok = lab.isin(IMMUNE_HIGH)
                cells.loc[ids, 'reference_status'] = np.where(
                    ok, 'candidate_immune_parent_supported', 'reference_parent_conflict_or_unresolved')
            if name == 'Immune_All_Low':
                conflict = lab.isin(EXPLICIT_NONIMMUNE)
                cells.loc[ids, 'Immune_All_Low_explicit_parent_conflict'] = conflict
                cells.loc[ids[conflict.to_numpy()], 'reference_status'] = 'reference_parent_conflict_or_unresolved'
            if name == 'Adult_Human_Vascular':
                family = lab.map(vascular_family)
                cells.loc[ids, 'reference_status'] = np.where(
                    family == cells.loc[ids, 'marker_parent'],
                    'candidate_vascular_parent_supported', 'reference_parent_conflict_or_unresolved')
            summary['models'][name] = {'sha256': sha(path), 'input_cells': len(ids),
                'features_shared': int(b.var_names.isin(model.classifier.features).sum()),
                'candidate_counts': lab.value_counts().to_dict()}
            print(name, len(ids), 'finished', flush=True)
    ids = cells.index[cells.prediction_route == 'singler_hpca']
    if len(ids):
        qdir = args.out / 'singler_query'; qdir.mkdir()
        q = sp.csr_matrix(a.X[a.obs_names.get_indexer(ids), :])
        mmwrite(str(qdir / 'query.mtx'), q.T)
        (qdir / 'genes.txt').write_text('\n'.join(a.var_names) + '\n')
        (qdir / 'cells.txt').write_text('\n'.join(ids) + '\n')
        worker = Path(__file__).with_name('singler_hpca.R')
        command = [str(args.rscript), '--vanilla', str(worker), str(qdir),
                   str(args.r_cache), os.environ['OMP_NUM_THREADS']]
        with (args.out / 'singler.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
        pred = pd.read_csv(qdir / 'predictions.tsv', sep='\t').set_index('cell_id')
        if not pred.index.is_unique or set(pred.index) != set(ids):
            raise ValueError('SingleR output cell IDs differ')
        pred = pred.reindex(ids)
        for col in pred:
            cells.loc[ids, col] = pred[col]
        hpca_parent = {'Fibroblasts': 'fibroblast', 'Endothelial_cells': 'endothelial',
                       'Smooth_muscle_cells': 'mural', 'Epithelial_cells': 'epithelial',
                       'Hepatocytes': 'hepatocyte', 'Astrocyte': 'glial'}
        ok = pred.hpca_main_pruned.map(hpca_parent) == cells.loc[ids, 'marker_parent']
        cells.loc[ids, 'reference_status'] = np.where(
            ok, 'candidate_hpca_parent_supported', 'reference_parent_conflict_or_pruned')
        summary['models']['SingleR_HPCA'] = {'input_cells': len(ids),
            'matrix_sha256': sha(args.r_cache / '1d748702c9eed_3508'),
            'labels_sha256': sha(args.r_cache / '1d7483d3c0df_3509'),
            'worker_sha256': sha(worker),
            'reference_label_alignment': 'paired EH3492/EH3493 positional order; colData rownames absent',
            'version_log': (args.out / 'singler.log').read_text().strip()}
    a.file.close()
    # Routing is immutable: reference output must never rewrite marker parents.
    assert not (cells.reference_status == 'pending').any()
    cells.to_csv(args.out / 'candidate_labels.tsv.gz', sep='\t')
    summary['marker_parent_counts'] = cells.marker_parent.value_counts().to_dict()
    summary['route_counts'] = cells.prediction_route.value_counts().to_dict()
    summary['reference_status_counts'] = cells.reference_status.value_counts().to_dict()
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: summary[k] for k in ['source', 'output_cells', 'route_counts',
                                             'reference_status_counts']}), flush=True)


if __name__ == '__main__':
    main()
