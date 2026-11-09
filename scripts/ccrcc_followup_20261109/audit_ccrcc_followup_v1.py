"""Independent aggregate/row reconciliation; no private rows in the receipt."""
from pathlib import Path
import sys
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix
from scipy.spatial.distance import cdist
from statsmodels.stats.multitest import multipletests

r = Path(sys.argv[1]).resolve()
o = r / 'followup_v1'
c = pd.read_parquet(o / 'identity_candidates_PRIVATE.parquet')
summary = json.loads((o / 'identity_summary.json').read_text())
with h5py.File(r / 'inputs/member_000_PRIVATE.h5') as f:
    m = f['matrix']
    x = csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:])).T.tocsr()
    names = np.array([v.decode() for v in m['features/name'][:]])
    types = np.array([v.decode() for v in m['features/feature_type'][:]])
    ids = [v.decode() for v in m['barcodes'][:]]
    factor = float(f.attrs['protein_scaling_factor'])
assert list(c.index) == ids and len(c) == summary['input_cells'] == 690322
assert c.QC_keep.sum() == summary['QC_keep']
rna = types == 'Gene Expression'
xn = x[:, rna]
gn = list(names[rna])
assert np.array_equal(np.asarray(xn.sum(1)).ravel(), c.RNA_total)
assert np.array_equal(xn.getnnz(1), c.RNA_genes)
gate = sum((xn[:, gn.index(g)].toarray().ravel() >= 1).astype(int) for g in ['PECAM1', 'VWF', 'EMCN']) >= 2
ecmask = c.QC_keep & c.identity_candidate.eq('Blood_EC')
assert np.all(gate[ecmask])
geom = np.load(o / 'geometry_PRIVATE.npz')
assert np.array_equal(np.flatnonzero(ecmask), geom['ec_global'])
xy = c[['x_centroid', 'y_centroid']].to_numpy()
query, e = geom['query_global'], geom['ec_global']
sample = np.random.default_rng(312).choice(len(query), 128, replace=False)
dist = cdist(xy[query[sample]], xy[e]).min(1)
assert np.allclose(dist, geom['centroid_distance'][sample], atol=1e-8, rtol=0)
ack = xn[:, gn.index('ACKR1')].toarray().ravel()[e] >= 1
assert np.array_equal(ack, geom['ACKR1_detected'])
marker = pd.read_csv(o / 'marker_expected_identity.tsv', sep='\t')
groups = np.where(c.QC_keep, c.identity_candidate, 'Excluded_QC_or_conflict')
max_marker_error = 0.
for _, row in marker.iterrows():
    z = groups == row.group
    assert z.sum() == row.n_cells
    assert (row.marker in gn) == bool(row.available)
    if row.available and z.any():
        v = xn[z, gn.index(row.marker)].toarray().ravel()
        errors = [abs((v >= 1).mean() - row.fraction_positive),
                  abs(np.log1p(v).mean() - row.mean_log1p_RNA_UMI)]
        max_marker_error = max(max_marker_error, *errors)
assert max_marker_error < 1e-8
protein = pd.read_csv(o / 'protein_within_antibody.tsv', sep='\t')
xp = x[:, types == 'Protein Expression'].toarray() / factor
pg = list(names[types == 'Protein Expression'])
max_protein_error = 0.
for _, row in protein.iterrows():
    vals = np.quantile(xp[groups == row.group, pg.index(row.antibody)], [.25, .5, .75])
    expected = row[['q25_mean_intensity', 'median_mean_intensity', 'q75_mean_intensity']].to_numpy(float)
    max_protein_error = max(max_protein_error, float(np.max(abs(vals - expected))))
assert max_protein_error < 1e-8
table = pd.read_csv(o / 'spatial_complete96.tsv', sep='\t')
assert len(table) == 96
assert len(table[['geometry', 'radius_um', 'receiver', 'denominator', 'tile_um']].drop_duplicates()) == 96
assert set(table.geometry) == {'centroid', 'EC_cell_polygon'}
assert set(table.radius_um) == {25, 50, 100, 200}
assert set(table.receiver) == {'B', 'T', 'NK'}
assert set(table.tile_um) == {250, 500}
null = np.load(o / 'null_logodds_PRIVATE.npz')
weights = np.load(o / 'endpoint_weights_PRIVATE.npz')
max_spatial_error = 0.
label = c.identity_candidate.to_numpy()[query]
for tile in [250, 500]:
    tab = table[table.tile_um == tile]
    for endpoint, (_, row) in enumerate(tab.iterrows()):
        prefix = 'centroid' if row.geometry == 'centroid' else 'polygon'
        distance = geom[prefix + '_distance']
        owner = geom[prefix + '_owner']
        near = np.isfinite(distance) & (owner >= 0) & (distance <= row.radius_um)
        d = near if row.denominator == 'allQC_nonEC' else near & np.isin(label, ['B', 'T', 'NK'])
        y = near & (label == row.receiver)
        v = np.zeros(len(query), dtype=bool); v[owner >= 0] = ack[owner[owner >= 0]]
        values = [int((y & v).sum()), int((d & v).sum()), int((y & ~v).sum()), int((d & ~v).sum())]
        assert values == [row.target_y, row.target_denominator, row.control_y, row.control_denominator]
        yv = np.bincount(owner[y], minlength=len(e))
        dv = np.bincount(owner[d], minlength=len(e))
        assert np.array_equal(yv, weights['Y'][endpoint])
        assert np.array_equal(dv, weights['D'][endpoint])
        y1, d1, y0, d0 = values
        statistic = np.log((y1 + .5) / (d1 - y1 + .5) / ((y0 + .5) / (d0 - y0 + .5)))
        max_spatial_error = max(max_spatial_error, abs(statistic - row.log_odds))
        p = (1 + np.count_nonzero(null[str(tile)][:, endpoint] >= statistic - 1e-12)) / 1000
        # Numerical ties are compared with the exact producer statistic below.
        producer_obs = np.log((yv @ ack + .5) / (dv @ ack - yv @ ack + .5)) - np.log((yv @ ~ack + .5) / (dv @ ~ack - yv @ ~ack + .5))
        p = (1 + np.count_nonzero(null[str(tile)][:, endpoint] >= producer_obs)) / 1000
        assert abs(p - row.p_one_sided_conditional) < 1e-12
for method, column in [('fdr_bh', 'BH_q_complete96'), ('fdr_by', 'BY_q_complete96')]:
    expected = multipletests(table.p_one_sided_conditional, method=method)[1]
    max_spatial_error = max(max_spatial_error, float(np.max(abs(expected - table[column]))))
assert max_spatial_error < 1e-8
assert (table.p_one_sided_conditional >= .001).all()
cluster_checks = []
import anndata as ad
from sklearn.metrics import adjusted_rand_score
for folder in sorted(o.glob('clusters_*')):
    obj = ad.read_h5ad(folder / 'RNA_only_cluster_object_PRIVATE.h5ad')
    assert 'ACKR1' not in obj.var_names and set(obj.var_names).issubset(set(gn))
    assert len(set(obj.obs.source_index_PRIVATE)) == obj.n_obs
    for _, row in pd.read_csv(folder / 'seed_stability.tsv', sep='\t').iterrows():
        got = adjusted_rand_score(obj.obs['seed_' + str(int(row.seed_a))], obj.obs['seed_' + str(int(row.seed_b))])
        assert abs(got - row.ARI) < 1e-12
    cluster_checks.append({'pool': folder.name.removeprefix('clusters_'), 'cells': obj.n_obs,
                           'ACKR1_and_proteins_excluded': True, 'ARI_recomputed': True})
receipt = {'UTC': datetime.now(timezone.utc).isoformat(), 'audit_pass': True,
           'checks': ['690322 IDs/order and RNA-only depth', 'ACKR1-independent EC gate',
                      '128 independent brute-force nearest centroid distances',
                      'all marker aggregate values', 'all antibody quantiles / stored scale',
                      'complete96 endpoint cell counts and nearest-EC weights',
                      '999-null add-one P and independent statsmodels BH/BY',
                      'RNA-only clustering features and seed ARI'],
           'protein_scaling_factor': factor, 'max_marker_error': max_marker_error,
           'max_protein_quantile_error': max_protein_error, 'max_spatial_error': max_spatial_error,
           'cluster_checks': cluster_checks,
           'scope': 'Software/arithmetic, not final identity or patient-level biological validation'}
(o / 'independent_audit.json').write_text(json.dumps(receipt, indent=2))
print(json.dumps(receipt))
