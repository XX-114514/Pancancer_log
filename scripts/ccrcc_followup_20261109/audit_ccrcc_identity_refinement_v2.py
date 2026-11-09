"""Independent raw-marker/protein and subset-membership audit, no label freezing."""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix
import anndata as ad
from sklearn.metrics import adjusted_rand_score

r = Path(sys.argv[1]).resolve(); o = r / 'identity_refinement_v2'
assert (o / 'complete.json').exists()
parent = ad.read_h5ad(r / 'followup_v1/clusters_immune/RNA_only_cluster_object_PRIVATE.h5ad')
parent_ix = parent.obs.source_index_PRIVATE.to_numpy(int)
parent_label = parent.obs.seed_0.to_numpy().astype(str)
main = pd.read_parquet(r / 'followup_v1/identity_candidates_PRIVATE.parquet')
with h5py.File(r / 'inputs/member_000_PRIVATE.h5') as f:
    m = f['matrix']; x = csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:])).T.tocsr()
    names = np.array([v.decode() for v in m['features/name'][:]])
    types = np.array([v.decode() for v in m['features/feature_type'][:]])
    scale = float(f.attrs['protein_scaling_factor'])
rna_names = list(names[types == 'Gene Expression'])
protein_names = list(names[types == 'Protein Expression'])
reviews = []; max_error = 0.
for name, selected in [('T_NK', ['0', '1']), ('B_ASC', ['6', '7'])]:
    folder = o / name
    obj = ad.read_h5ad(folder / 'RNA_only_refinement_PRIVATE.h5ad')
    ix = obj.obs.source_index_PRIVATE.to_numpy(int)
    expected = parent_ix[np.isin(parent_label, selected)]
    assert np.array_equal(ix, expected), 'Subset must match marker-selected parent cells exactly'
    assert 'ACKR1' not in obj.var_names and set(obj.var_names).issubset(set(rna_names))
    assert np.array_equal(obj.obs.identity_rule_candidate.astype(str).to_numpy(), main.identity_candidate.to_numpy()[ix])
    assert np.array_equal(obj.obs.main_QC_keep.to_numpy(), main.QC_keep.to_numpy()[ix])
    group = obj.obs.seed_0.astype(str).to_numpy()
    raw = x[ix][:, types == 'Gene Expression'].toarray()
    prot = x[ix][:, types == 'Protein Expression'].toarray() / scale
    markers = pd.read_csv(folder / 'marker_expected_identity.tsv', sep='\t', dtype={'group': str})
    for _, row in markers.iterrows():
        z = group == row.group
        assert int(z.sum()) == row.n_cells
        assert bool(row.available) == (row.marker in rna_names)
        if row.available:
            v = raw[z, rna_names.index(row.marker)]
            max_error = max(max_error, abs((v >= 1).mean() - row.fraction_positive),
                            abs(np.log1p(v).mean() - row.mean_log1p_RNA_UMI))
    proteins = pd.read_csv(folder / 'protein_within_antibody.tsv', sep='\t', dtype={'group': str})
    for _, row in proteins.iterrows():
        vals = np.quantile(prot[group == row.group, protein_names.index(row.antibody)], [.25, .5, .75])
        target = row[['q25_mean_intensity', 'median_mean_intensity', 'q75_mean_intensity']].to_numpy(float)
        max_error = max(max_error, float(np.max(abs(vals - target))))
    for _, row in pd.read_csv(folder / 'seed_stability.tsv', sep='\t').iterrows():
        ari = adjusted_rand_score(obj.obs['seed_' + str(int(row.seed_a))], obj.obs['seed_' + str(int(row.seed_b))])
        max_error = max(max_error, abs(ari - row.ARI))
    review = pd.read_csv(folder / 'cluster_review.tsv', sep='\t', dtype={'cluster': str})
    for _, row in review.iterrows():
        z = group == row.cluster; assert z.sum() == row.n_cells
        for typ in ['B', 'T', 'NK', 'Blood_EC', 'Lymphatic', 'Mural', 'Other']:
            assert (obj.obs.identity_rule_candidate.astype(str).to_numpy()[z] == typ).sum() == row['candidate_' + typ]
    reviews.append({'pool': name, 'cells': obj.n_obs, 'clusters': len(review),
                    'parent_membership_exact': True, 'raw_marker_and_protein_quantiles_checked': True})
assert max_error < 1e-8
result = {'UTC': datetime.now(timezone.utc).isoformat(), 'audit_pass': True,
          'max_abs_error': float(max_error), 'protein_scaling_factor': scale, 'reviews': reviews,
          'checks': ['exact marker-selected parent membership', 'ACKR1/protein/position excluded from RNA features',
                     'original candidate/QC preservation', 'all measured/unmeasured marker aggregates',
                     'all unscaled within-antibody quantiles', 'seed ARI and cluster candidate counts'],
          'scope': 'Subset arithmetic/feature audit, not a frozen biological label acceptance'}
(o / 'independent_audit.json').write_text(json.dumps(result, indent=2))
print(json.dumps(result))
