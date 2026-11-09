"""Independently reconcile cluster intersections, cell counts and joint192 FDR."""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys
import hashlib
import numpy as np
import pandas as pd
import h5py
import anndata as ad
from scipy.sparse import csc_matrix
from scipy.spatial.distance import cdist
from statsmodels.stats.multitest import multipletests

r = Path(sys.argv[1]).resolve(); o = r / 'program_concordance_v3'
assert (o / 'complete.json').exists()
original = pd.read_parquet(r / 'followup_v1/identity_candidates_PRIVATE.parquet')
current = pd.read_parquet(o / 'identity_candidates_PRIVATE.parquet')
assert list(current.index) == list(original.index)
assert np.array_equal(current.QC_keep, original.QC_keep)
assert np.array_equal(current[['x_centroid', 'y_centroid']], original[['x_centroid', 'y_centroid']])
definitions = [('B', 'identity_refinement_v2/B_ASC/RNA_only_refinement_PRIVATE.h5ad', ['0']),
               ('T', 'identity_refinement_v2/T_NK/RNA_only_refinement_PRIVATE.h5ad', ['0', '1']),
               ('NK', 'identity_refinement_v2/T_NK/RNA_only_refinement_PRIVATE.h5ad', ['2']),
               ('Blood_EC', 'followup_v1/clusters_vascular_stromal/RNA_only_cluster_object_PRIVATE.h5ad', ['0'])]
counts = []
for typ, path, reference_clusters in definitions:
    obj = ad.read_h5ad(r / path, backed='r')
    indices = obj.obs.source_index_PRIVATE.to_numpy(int)
    reference = obj.obs.seed_0.astype(str)
    agree = reference.isin(reference_clusters).to_numpy()
    jaccards = []
    for seed in [1, 2]:
        target = obj.obs['seed_' + str(seed)].astype(str)
        cross = pd.crosstab(reference, target)
        choices = []
        for cluster in reference_clusters:
            sizes_a = (reference == cluster).sum()
            size_b = target.value_counts()
            candidates = [(float(cross.loc[cluster, b]) / (sizes_a + size_b[b] - cross.loc[cluster, b]), b)
                          for b in sorted(set(target), key=int)]
            best = sorted(candidates, key=lambda v: (-v[0], int(v[1])))[0]
            choices.append(best[1]); jaccards.append(best[0])
        agree &= target.isin(choices).to_numpy()
    expected = np.zeros(len(current), dtype=bool); expected[indices[agree]] = True
    expected &= original.QC_keep.to_numpy() & original.identity_candidate.eq(typ).to_numpy()
    observed = current.QC_keep & current.identity_candidate.eq(typ)
    assert np.array_equal(expected, observed), 'Consensus membership mismatch: ' + typ
    counts.append({'identity': typ, 'cells': int(observed.sum()), 'minimum_cluster_Jaccard': min(jaccards)})
    obj.file.close()
with h5py.File(r / 'inputs/member_000_PRIVATE.h5') as f:
    m = f['matrix']; names = [x.decode() for x in m['features/name'][:]]
    types = [x.decode() for x in m['features/feature_type'][:]]
    j = next(i for i, (name, kind) in enumerate(zip(names, types)) if name == 'ACKR1' and kind == 'Gene Expression')
    x = csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:])).T.tocsr()
    ack = x[:, j].toarray().ravel() >= 1
geom = np.load(o / 'geometry_PRIVATE.npz'); weights = np.load(o / 'endpoint_weights_PRIVATE.npz')
ec = np.flatnonzero(current.QC_keep & current.identity_candidate.eq('Blood_EC'))
assert np.array_equal(ec, geom['ec_global'])
assert np.array_equal(ack[ec], geom['ACKR1_detected'])
query = geom['query_global']; label = current.identity_candidate.to_numpy()[query]
xy = current[['x_centroid', 'y_centroid']].to_numpy()
pick = np.random.default_rng(331).choice(len(query), 128, replace=False)
assert np.allclose(cdist(xy[query[pick]], xy[ec]).min(1), geom['centroid_distance'][pick], atol=1e-8, rtol=0)
t = pd.read_csv(o / 'spatial_complete96.tsv', sep='\t'); null = np.load(o / 'null_logodds_PRIVATE.npz')
assert len(t) == 96
for tile in [250, 500]:
    z = t[t.tile_um == tile]
    for endpoint, (_, row) in enumerate(z.iterrows()):
        prefix = 'centroid' if row.geometry == 'centroid' else 'polygon'
        owner, distance = geom[prefix + '_owner'], geom[prefix + '_distance']
        near = (owner >= 0) & np.isfinite(distance) & (distance <= row.radius_um)
        d = near if row.denominator == 'allQC_nonEC' else near & np.isin(label, ['B', 'T', 'NK'])
        y = near & (label == row.receiver)
        yv = np.bincount(owner[y], minlength=len(ec)); dv = np.bincount(owner[d], minlength=len(ec))
        assert np.array_equal(yv, weights['Y'][endpoint]) and np.array_equal(dv, weights['D'][endpoint])
        y1, d1, y0, d0 = yv[ack[ec]].sum(), dv[ack[ec]].sum(), yv[~ack[ec]].sum(), dv[~ack[ec]].sum()
        assert [y1, d1, y0, d0] == [row.target_y, row.target_denominator, row.control_y, row.control_denominator]
        statistic = np.log((y1 + .5) / (d1 - y1 + .5)) - np.log((y0 + .5) / (d0 - y0 + .5))
        assert abs(statistic - row.log_odds) < 1e-10
        p = (1 + (null[str(tile)][:, endpoint] >= statistic).sum()) / 1000
        assert abs(p - row.p_one_sided_conditional) < 1e-12
both = pd.read_csv(o / 'complete192_two_views.tsv', sep='\t')
assert len(both) == 192
old = pd.read_csv(r / 'followup_v1/spatial_complete96.tsv', sep='\t')
pd.testing.assert_frame_equal(both.iloc[:96][old.columns].reset_index(drop=True), old, check_dtype=False, atol=1e-12, rtol=1e-12)
pd.testing.assert_frame_equal(both.iloc[96:][t.columns].reset_index(drop=True), t, check_dtype=False, atol=1e-12, rtol=1e-12)
error = 0.
for method, column in [('fdr_bh', 'BH_q_joint192'), ('fdr_by', 'BY_q_joint192')]:
    q = multipletests(both.p_one_sided_conditional, method=method)[1]
    error = max(error, float(np.max(abs(q - both[column]))))
assert error < 1e-12
receipt = {'UTC': datetime.now(timezone.utc).isoformat(), 'audit_pass': True,
           'checks': ['all3-seed intersections independently reconstructed with crosstabs',
                      'QC and coordinates unchanged', '128 brute-force nearest EC distances',
                      'all96 endpoint cell counts/weights/log-odds/add-one P',
                      'original96 preserved and joint192 independent statsmodels BH/BY'],
           'counts': counts, 'max_joint_FDR_error': error, 'joint_hypotheses': 192,
           'original96_SHA': hashlib.sha256((r / 'followup_v1/spatial_complete96.tsv').read_bytes()).hexdigest(),
           'scope': 'Exploratory identity sensitivity, not final labels or independent patient replication'}
(o / 'independent_audit.json').write_text(json.dumps(receipt, indent=2)); print(json.dumps(receipt))
