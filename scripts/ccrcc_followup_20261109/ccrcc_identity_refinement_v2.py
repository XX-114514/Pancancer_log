"""Marker-guided subset review after the broad RNA-only diagnostic pools.

Selection is based on RNA identity programs, never on spatial effect/P values.
It does not alter the original 96-test family or freeze new cell identities.
"""
from pathlib import Path
from datetime import datetime, timezone
import sys
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.metrics import adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import ccrcc_followup_v1 as base
from resumable_public_zip_v1 import atomic_json, exclusive_lock

r = Path(sys.argv[1]).resolve()
o = r / 'identity_refinement_v2'


def main():
    o.mkdir(exist_ok=False)
    design = {'UTC_before_subset_effects': datetime.now(timezone.utc).isoformat(),
              'parent': 'immune pool seed0; identities from RNA markers/protein review',
              'B_ASC_parent_clusters': ['6', '7'], 'T_NK_parent_clusters': ['0', '1'],
              'selection_basis': 'MS4A1/CD79A/BANK1 vs MZB1/TENT5C/SDC1; CD3/T/cytotoxic programs',
              'spatial_effect_used_for_selection': False,
              'spatial_complete96_changed': False,
              'ACKR1_or_protein_or_position_clustering_features': False,
              'normalization': 'RNA excl ACKR1 CP10K/log1p, PCA30, neighbors30, Leiden0.5 seeds0/1/2',
              'scope': 'Exploratory identity diagnostics, not final labels or an independent replication'}
    atomic_json(o / 'design.json', design)
    ids, c, xy, xr, rn, xp, pn, scale, total, ng, a = base.load(r)
    pos, count, masks, conflict, depth, keep, label = base.identities(a, total, ng)
    parent = ad.read_h5ad(r / 'followup_v1/clusters_immune/RNA_only_cluster_object_PRIVATE.h5ad')
    original_index = parent.obs.source_index_PRIVATE.to_numpy(int)
    parent_group = parent.obs.seed_0.to_numpy().astype(str)
    assert len(set(original_index)) == len(original_index)
    base.PANELS['ASC_review'] = ['MZB1', 'SDC1', 'CD38', 'TENT5C', 'SLAMF7', 'IGHG2', 'IGKC', 'JCHAIN']
    base.PANELS['T_NK_review'] = ['CD2', 'CD4', 'CD8A', 'CD8B', 'KLRB1', 'FCGR3A', 'GZMB', 'GZMK', 'CCL5', 'TRDC', 'TRGC1', 'TRGC2']
    for g in base.PANELS['ASC_review'] + base.PANELS['T_NK_review']:
        if g in rn and g not in a:
            a[g] = xr[:, list(rn).index(g)].toarray().ravel()
    result = []
    coverage = []
    for name, parents in [('T_NK', ['0', '1']), ('B_ASC', ['6', '7'])]:
        ix = original_index[np.isin(parent_group, parents)]
        for typ in ['B', 'T', 'NK']:
            coverage.append({'review_pool': name, 'candidate': typ,
                             'main_candidate_total': int((keep & (label == typ)).sum()),
                             'main_candidate_in_review_pool': int((keep[ix] & (label[ix] == typ)).sum())})
        folder = o / name; folder.mkdir()
        features = rn != 'ACKR1'
        obj = ad.AnnData(xr[ix][:, features].astype(np.float32))
        obj.var_names = rn[features]
        obj.obs['source_index_PRIVATE'] = ix
        obj.obs['identity_rule_candidate'] = label[ix]
        obj.obs['conflicting_programs'] = conflict[ix]
        obj.obs['main_QC_keep'] = keep[ix]
        sc.pp.normalize_total(obj, target_sum=1e4); sc.pp.log1p(obj)
        variable = np.asarray(obj.X.power(2).mean(0)).ravel() - np.asarray(obj.X.mean(0)).ravel() ** 2 > 1e-9
        obj = obj[:, variable].copy()
        sc.pp.pca(obj, n_comps=30, random_state=0, svd_solver='arpack')
        sc.pp.neighbors(obj, n_neighbors=30, use_rep='X_pca', random_state=0)
        for seed in [0, 1, 2]:
            sc.tl.leiden(obj, resolution=.5, random_state=seed,
                         key_added='seed_' + str(seed), flavor='leidenalg', n_iterations=-1)
        stability = [{'seed_a': i, 'seed_b': j, 'ARI': float(adjusted_rand_score(obj.obs['seed_' + str(i)], obj.obs['seed_' + str(j)]))}
                     for i, j in [(0, 1), (0, 2), (1, 2)]]
        pd.DataFrame(stability).to_csv(folder / 'seed_stability.tsv', sep='\t', index=False)
        group = obj.obs.seed_0.to_numpy().astype(str)
        groups = sorted(set(group), key=int)
        marker = base.marker_table({g: v[ix] for g, v in a.items()}, group, groups)
        protein = base.protein_table(xp[ix], pn, group, groups)
        marker.to_csv(folder / 'marker_expected_identity.tsv', sep='\t', index=False, na_rep='NA')
        protein.to_csv(folder / 'protein_within_antibody.tsv', sep='\t', index=False, na_rep='NA')
        base.marker_plot(marker, folder / 'marker_expected_identity.png', name + ' subset RNA-only review; no frozen labels')
        base.protein_plot(protein, folder / 'protein_within_antibody.png')
        means = np.asarray(obj.X.mean(0)).ravel(); rows = []
        for g in groups:
            z = group == g; gm = np.asarray(obj.X[z].mean(0)).ravel()
            other = (means * len(group) - gm * z.sum()) / (len(group) - z.sum())
            top = obj.var_names[np.argsort(gm - other)[-10:][::-1]].tolist()
            row = {'cluster': g, 'n_cells': int(z.sum()), 'main_QC_keep_cells': int(keep[ix][z].sum()),
                   'conflicting_program_cells': int(conflict[ix][z].sum()),
                   'top_mean_logexpr_difference_markers': ';'.join(top)}
            for typ in ['B', 'T', 'NK', 'Blood_EC', 'Lymphatic', 'Mural', 'Other']:
                row['candidate_' + typ] = int((label[ix][z] == typ).sum())
            rows.append(row)
        pd.DataFrame(rows).to_csv(folder / 'cluster_review.tsv', sep='\t', index=False)
        obj.write_h5ad(folder / 'RNA_only_refinement_PRIVATE.h5ad', compression='gzip')
        rec = {'review_pool': name, 'cells': len(ix), 'clusters': len(groups),
               'ARI_min': min(s['ARI'] for s in stability), 'parent_cluster_ids': parents,
               'label_freeze_performed': False, 'spatial_complete96_changed': False}
        atomic_json(folder / 'complete.json', rec); result.append(rec)
    pd.DataFrame(coverage).to_csv(o / 'candidate_coverage.tsv', sep='\t', index=False)
    atomic_json(o / 'complete.json', {'UTC': datetime.now(timezone.utc).isoformat(),
                                    'status': 'subset_review_computed_not_frozen', 'pools': result,
                                    'spatial_models_refitted': False, 'source_and_V4_labels_preserved': True})


with exclusive_lock(r / 'identity_refinement_v2.lock'):
    main()
