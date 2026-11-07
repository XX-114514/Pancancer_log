"""Independent neighbour, feature and saved-cluster checks, without CM effects."""
import json, sys, hashlib
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import shortest_path

r = Path(sys.argv[1])
key = sys.argv[2]
p = r / 'niches_retry2' / key
receipt = json.loads((p / 'complete.json').read_text())
spots = pd.read_csv(p / 'spot_niches_PRIVATE.tsv.gz', sep='\t', index_col=0)
eligibility = pd.read_csv(p / 'spot_niche_eligibility_PRIVATE.tsv.gz', sep='\t', index_col=0)
assert len(eligibility) == receipt['n_spots_before_niche_information_filter']
assert int((~eligibility.niche_eligible_nonzero_HVG).sum()) == receipt['n_spots_zero_in_selected_HVG_unassessable_for_niche']
assert spots.index.equals(eligibility.index[eligibility.niche_eligible_nonzero_HVG])
latent = np.load(p / 'scVI_latent_PRIVATE.npy')
features = np.load(p / 'CellCharter_features_PRIVATE.npy')
graph = sp.load_npz(p / 'spatial_6NN_PRIVATE.npz').tocsr()
assert len(spots) == receipt['n_spots'] and spots.index.is_unique
assert latent.shape == (len(spots), 10) and features.shape == (len(spots), 40)
assert np.isfinite(latent).all() and np.isfinite(features).all()
assert graph.shape == (len(spots), len(spots)) and np.all(graph.getnnz(axis=1) == 6)
assert graph.diagonal().sum() == 0
xy = spots[['x', 'y']].to_numpy()
distances, _ = cKDTree(xy).query(xy, k=7)
ii, jj = graph.nonzero()
observed = np.linalg.norm(xy[ii] - xy[jj], axis=1)
assert np.all(observed <= distances[:, -1][ii] + 1e-6)
# Independent shortest-path shells, permitting exact-distance ties in the 6NN graph.
errors = []
for start in range(0, len(spots), 250):
    indices = np.arange(start, min(start + 250, len(spots)))
    distances = shortest_path(graph, directed=True, unweighted=True, indices=indices)
    blocks = [latent[indices]]
    for hop in [1, 2, 3]:
        shell = sp.csr_matrix(distances == hop, dtype=float)
        denominator = np.asarray(shell.sum(axis=1)).ravel()
        block = np.asarray(shell @ latent) / np.maximum(denominator, 1)[:, None]
        blocks.append(block)
    expected = np.concatenate(blocks, axis=1)
    errors.append(float(np.max(np.abs(expected - features[indices]))))
assert max(errors) < 1e-4, max(errors)
history = np.load(p / 'scVI_ELBO_history_PRIVATE.npy')
assert len(history) == receipt['scVI_epochs_recorded'] and np.isfinite(history).all()
assert 1 <= len(history) <= 400
stability = pd.read_csv(p / 'K_stability_PRIVATE.tsv', sep='\t')
assert set(stability.k) == set(range(2, 11))
assert np.isfinite(stability.stability_mean).all()
selected = int(stability.loc[stability.stability_mean.idxmax(), 'k'])
assert selected == receipt['best_k']
assert spots.niche.nunique() == receipt['n_observed_niches']
assert hashlib.sha256((p / 'spot_niches_PRIVATE.tsv.gz').read_bytes()).hexdigest() == receipt['spot_labels_SHA256']
assert not receipt['uses_CM_labels_abundance_activity_or_clinical_outcomes_to_fit']
proof = {'status': 'PASS_independent_CellCharter_saved_features_6NN_distance_shells_and_K_selection',
    'UTC': datetime.now(timezone.utc).isoformat(), 'section_code': key, 'dataset': receipt['dataset'],
    'n_spots': len(spots), 'best_k': selected, 'n_observed_niches': spots.niche.nunique(),
    'max_feature_error': max(errors), 'scVI_epochs': len(history),
    'independent_identity_ground_truth_or_biological_validation_claim': False}
(p / 'independent_audit.json').open('x').write(json.dumps(proof, indent=2))
print(json.dumps(proof))
