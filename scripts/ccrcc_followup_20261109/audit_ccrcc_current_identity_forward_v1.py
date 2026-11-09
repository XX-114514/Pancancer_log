"""Counts/geometry/P/FDR audit only; never recomputes or changes identities."""
from pathlib import Path
from datetime import datetime, timezone
import sys
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from statsmodels.stats.multitest import multipletests

r = Path(sys.argv[1]).resolve(); o = r / 'current_identity_forward_v1'
assert (o / 'complete.json').exists()
identity_path = r / 'spatial_pilot_v1/identity_candidates_PRIVATE.tsv.gz'
identity = pd.read_csv(identity_path, sep='\t', index_col=0)
design = json.loads((o / 'design.json').read_text())
assert hashlib.sha256(identity_path.read_bytes()).hexdigest() == design['identity_source_SHA']
keep = identity.QC_keep.to_numpy(bool); label = identity.identity_candidate.to_numpy()
assert keep.sum() == 358841
cells = pd.read_parquet(r / 'inputs/member_001_PRIVATE.parquet').set_index('cell_id').loc[identity.index]
xy = cells[['x_centroid', 'y_centroid']].to_numpy()
g = np.load(o / 'graph_PRIVATE.npz'); es, er, dd = g['sender'], g['receiver'], g['distance']
immune = g['immune_indices']
assert np.array_equal(er, immune[g['receiver_row']]) and np.all(es != er)
assert keep[es].all() and keep[er].all() and np.all(dd <= 100)
assert np.allclose(np.linalg.norm(xy[es] - xy[er], axis=1), dd, atol=1e-8, rtol=0)
sample = np.random.default_rng(940).choice(len(immune), 16, replace=False)
qc = np.flatnonzero(keep); maxerror = 0.
for row in sample:
    distance = cdist(xy[immune[row]:immune[row] + 1], xy[qc]).ravel()
    distance[qc == immune[row]] = np.inf
    brute = np.sort(np.partition(distance, 5)[:6]); brute = brute[brute <= 100]
    observed = np.sort(dd[g['receiver_row'] == row])
    assert len(brute) == len(observed)
    if len(brute):maxerror = max(maxerror, float(np.max(abs(brute - observed))))
assert maxerror < 1e-8
ec = np.flatnonzero(keep & (label == 'Blood_EC'))
mapping = np.full(len(identity), -1, dtype=int); mapping[ec] = np.arange(len(ec))
comp = pd.read_csv(o / 'CM03_existing_identity_spatial12.tsv', sep='\t')
comp_null = np.load(o / 'composition_null_PRIVATE.npz'); assert len(comp) == 12
for _, row in comp.iterrows():
    edges = (mapping[es] >= 0) & (dd <= row.radius_um)
    b = np.bincount(mapping[es[edges & (label[er] == 'B')]], minlength=len(ec))
    nk = np.bincount(mapping[es[edges & (label[er] == 'NK')]], minlength=len(ec))
    observed = {'EC_B_edges': int(b.sum()), 'EC_NK_edges': int(nk.sum()),
                'same_EC_B_and_NK': int(((b > 0) & (nk > 0)).sum())}[row.relation]
    assert observed == row.observed_count
    j = [50, 100].index(int(row.radius_um)) * 3 + ['EC_B_edges', 'EC_NK_edges', 'same_EC_B_and_NK'].index(row.relation)
    null = comp_null[str(row.tile_um)][:, j]
    assert abs((1 + (null >= observed).sum()) / 1000 - row.p_one_sided_conditional) < 1e-12
lr = pd.read_csv(o / 'spatial_expression_communication_complete108.tsv', sep='\t')
assert len(lr) == 108 and len(lr[['sender','receiver','ligand','receptor','radius_um','tile_um']].drop_duplicates()) == 108
inp = np.load(o / 'communication_inputs_PRIVATE.npz'); lr_null = np.load(o / 'communication_null_PRIVATE.npz')
genes = {'CCL19': 0, 'CXCL12': 1, 'CCR7': 2, 'CXCR4': 3}
for _, row in lr.iterrows():
    pair = ['Blood_EC', 'Mural', 'Lymphatic'].index(row.sender) * 3 + ['B', 'T', 'NK'].index(row.receiver)
    z = (inp['pair_code'] == pair) & (inp['distance'] <= row.radius_um)
    assert z.sum() == row.spatial_edges
    if row.ligand == 'CCL21':
        assert not row.components_measured and row.p_one_sided_conditional == 1 and np.isnan(row.observed_fraction)
        continue
    obs = int((z & inp['expression'][inp['es'], genes[row.ligand]] & inp['expression'][inp['er'], genes[row.receptor]]).sum())
    assert obs == row.observed_dualpositive_edges
    j = ([50, 100].index(int(row.radius_um)) * 2 + ['CCL19', 'CXCL12'].index(row.ligand)) * 9 + pair
    null = lr_null[str(row.tile_um)][:, j]
    p = (1 + (null >= obs).sum()) / 1000 if row.available else 1
    assert abs(p - row.p_one_sided_conditional) < 1e-12
max_fdr = 0.
for table, n in [(comp, 12), (lr, 108)]:
    assert len(table) == n
    for method, prefix in [('fdr_bh', 'BH'), ('fdr_by', 'BY')]:
        q = multipletests(table.p_one_sided_conditional, method=method)[1]
        max_fdr = max(max_fdr, float(np.max(abs(q - table[f'{prefix}_q_complete{n}']))))
assert max_fdr < 1e-12
metadata = pd.read_csv(o / 'imaging/ccRCC_current_identity/metadata_PRIVATE.tsv', sep='\t', index_col=0)
assert len(metadata) == 151833
assert np.array_equal(metadata.broad.to_numpy(), identity.loc[metadata.index].identity_candidate.to_numpy())
assert identity.loc[metadata.index].QC_keep.all()
result = {'UTC': datetime.now(timezone.utc).isoformat(), 'audit_pass': True,
          'identity_source_SHA_preserved': design['identity_source_SHA'], 'original_identity_changed': False,
          'checks': ['source identity SHA/QC/names unchanged', 'full graph Euclidean distances and 16 brute-force nearest6 sets',
                     '12 CM composition observed counts and exact plus-one P', '108 LR counts/missing slots/P',
                     'independent complete-family BH/BY', '151833 CellChat input labels match existing identities'],
          'max_graph_error_um': maxerror, 'max_FDR_error': max_fdr,
          'scope': 'Numerical audit, not reannotation or patient-level/causal validation'}
(o / 'independent_forward_audit.json').write_text(json.dumps(result, indent=2)); print(json.dumps(result))
