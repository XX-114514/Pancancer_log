"""Continue ccRCC CM/communication with the EXISTING identity snapshot.

No reannotation, subset-cluster gate, or state prerequisite. All cell IDs,
coordinates, matrices and graph edges stay private on the research server.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.io import mmwrite
from scipy.sparse import diags
import ccrcc_followup_v1 as source
from resumable_public_zip_v1 import atomic_json, exclusive_lock, sha256

r = Path(sys.argv[1]).resolve(); o = r / 'current_identity_forward_v1'
SEND = ['Blood_EC', 'Mural', 'Lymphatic']; RECV = ['B', 'T', 'NK']
AXES = [('CCL19', 'CCR7'), ('CXCL12', 'CXCR4'), ('CCL21', 'CCR7')]
RADII = [50, 100]; TILES = [250, 500]


def strata(xy, total, tile, extra=None):
    keys = np.column_stack([np.floor(xy / tile).astype(int), np.floor(np.log2(total)).astype(int)])
    if extra is not None:
        keys = np.column_stack([keys, extra])
    _, inv = np.unique(keys, axis=0, return_inverse=True)
    order = np.argsort(inv, kind='stable')
    return np.split(order, np.flatnonzero(np.diff(inv[order])) + 1)


def main():
    o.mkdir(exist_ok=False)
    identity_path = r / 'spatial_pilot_v1/identity_candidates_PRIVATE.tsv.gz'
    design = {'UTC_before_forward_effects': datetime.now(timezone.utc).isoformat(),
              'identity': 'Existing original ccRCC spatial_pilot_v1; no new identity gates',
              'identity_source_SHA': sha256(identity_path),
              'user_steering': 'Record identity concerns; continue with existing identities',
              'identity_review_is_prerequisite': False, 'identity_or_V4_labels_changed': False,
              'graph': 'Each B/T/NK receiver: 6 nearest original-QC cells, exclude self; retain edges <=50/100um',
              'composition_family': 'EC-B, EC-NK, same EC has B AND NK; 2 radii x2 tile/depth nulls =12',
              'composition_null': '999 within-tile x floor(log2 RNA depth) immune-identity shuffles; EC/positions fixed',
              'communication_family': '3 senders x3 receivers x3 axes x2 radii x2 local grids =108',
              'communication_null': '999 whole measured 4-gene binary profile shuffles within identity x tile x RNA-depth; geometry fixed',
              'P': '(1+null>=observed)/1000; ties included; missing CCL21 retained NA/P1',
              'families': 'Complete12 composition and complete108 LR, each BH/BY; one donor conditional exploration',
              'CellChat': 'Original identity six groups, original 477RNA CP10K/log1p; native micrometre geometry; no protein-as-UMI',
              'no_causal_or_independent_patient_replication_claim': True}
    atomic_json(o / 'design.json', design)
    ids, c, xy, xr, rn, xp, pn, scale, total, ng, a = source.load(r)
    d = pd.read_csv(identity_path, sep='\t', index_col=0)
    assert list(d.index) == list(ids)
    assert np.allclose(d[['x_centroid', 'y_centroid']], xy, rtol=0, atol=1e-5)
    keep = d.QC_keep.to_numpy(bool); label = d.identity_candidate.astype(str).to_numpy()
    assert keep.sum() == 358841
    expected = {'B': 7376, 'T': 57230, 'NK': 1244, 'Blood_EC': 49600, 'Lymphatic': 1682, 'Mural': 34701}
    assert {k: int((keep & (label == k)).sum()) for k in expected} == expected
    atomic_json(o / 'identity_preservation.json', {'original_identity_SHA': design['identity_source_SHA'],
                 'QC_keep': int(keep.sum()), 'identity_counts': expected,
                 'new_identity_or_cluster_filters': 0, 'positions_rescaled': False})
    qc_indices = np.flatnonzero(keep); immune = np.flatnonzero(keep & np.isin(label, RECV))
    distance, neighbour = cKDTree(xy[qc_indices]).query(xy[immune], k=7)
    global_neighbour = qc_indices[neighbour]
    # Sort-preserving self removal, including possible coincident-centroid ties.
    rank = np.cumsum(global_neighbour != immune[:, None], axis=1)
    valid = (global_neighbour != immune[:, None]) & (rank <= 6) & (distance <= 100)
    receiver_rows, positions = np.where(valid)
    edge_sender = global_neighbour[receiver_rows, positions]
    edge_receiver = immune[receiver_rows]
    edge_distance = distance[receiver_rows, positions]
    np.savez_compressed(o / 'graph_PRIVATE.npz', sender=edge_sender, receiver=edge_receiver,
                        receiver_row=receiver_rows, distance=edge_distance, immune_indices=immune)
    # Composition: retain exactly the same graph; shuffle immune names locally.
    ec_indices = np.flatnonzero(keep & (label == 'Blood_EC'))
    ec_index = np.full(len(ids), -1, dtype=int); ec_index[ec_indices] = np.arange(len(ec_indices))
    ec_edges = ec_index[edge_sender] >= 0
    ec_owner = ec_index[edge_sender[ec_edges]]; ir = receiver_rows[ec_edges]; ed = edge_distance[ec_edges]
    immune_codes = pd.Categorical(label[immune], categories=RECV).codes.astype(np.int64)
    def composition(codes):
        answer = []
        for radius in RADII:
            z = ed <= radius
            nb = np.bincount(ec_owner[z & (codes[ir] == 0)], minlength=len(ec_indices))
            nn = np.bincount(ec_owner[z & (codes[ir] == 2)], minlength=len(ec_indices))
            answer += [int(nb.sum()), int(nn.sum()), int(((nb > 0) & (nn > 0)).sum())]
        return np.array(answer, dtype=np.int64)
    observed = composition(immune_codes); comp_rows = []; comp_null = {}
    for tile in TILES:
        groups = strata(xy[immune], total[immune], tile)
        variable = [g for g in groups if len(np.unique(immune_codes[g])) > 1]
        rng = np.random.default_rng(2026110900 + tile)
        null = np.zeros((999, 6), dtype=np.int64)
        for rep in range(999):
            shuffled = immune_codes.copy()
            for g in variable:
                shuffled[g] = rng.permutation(immune_codes[g])
            null[rep] = composition(shuffled)
        comp_null[str(tile)] = null
        for j, (radius, relation) in enumerate((rad, name) for rad in RADII for name in ['EC_B_edges', 'EC_NK_edges', 'same_EC_B_and_NK']):
            comp_rows.append({'radius_um': radius, 'tile_um': tile, 'relation': relation,
                              'observed_count': int(observed[j]), 'null_mean': float(null[:, j].mean()),
                              'difference': float(observed[j] - null[:, j].mean()),
                              'p_one_sided_conditional': float((1 + (null[:, j] >= observed[j]).sum()) / 1000),
                              'EC_cells': len(ec_indices), 'n_permutations': 999})
    comp = pd.DataFrame(comp_rows); assert len(comp) == 12
    comp['BH_q_complete12'] = source.fdr(comp.p_one_sided_conditional)
    comp['BY_q_complete12'] = source.fdr(comp.p_one_sided_conditional, by=True)
    comp.to_csv(o / 'CM03_existing_identity_spatial12.tsv', sep='\t', index=False)
    np.savez_compressed(o / 'composition_null_PRIVATE.npz', **comp_null)
    atomic_json(o / 'composition_complete.json', {'complete_hypotheses': 12,
                'BH_lt05': int((comp.BH_q_complete12 < .05).sum()),
                'BY_lt05': int((comp.BY_q_complete12 < .05).sum()), 'identity_changed': False})
    # Expression-supported communication on the same original identity snapshot.
    included = np.flatnonzero(keep & np.isin(label, SEND + RECV))
    mapping = np.full(len(ids), -1, dtype=int); mapping[included] = np.arange(len(included))
    gene_names = ['CCL19', 'CXCL12', 'CCR7', 'CXCR4']
    gene_index = {g: j for j, g in enumerate(gene_names)}
    values = np.column_stack([xr[included, list(rn).index(g)].toarray().ravel() >= 1 for g in gene_names])
    pair_good = np.isin(label[edge_sender], SEND)
    es, er, dd = mapping[edge_sender[pair_good]], mapping[edge_receiver[pair_good]], edge_distance[pair_good]
    s_code = pd.Categorical(label[edge_sender[pair_good]], categories=SEND).codes.astype(np.int64)
    r_code = pd.Categorical(label[edge_receiver[pair_good]], categories=RECV).codes.astype(np.int64)
    pair_code = s_code * 3 + r_code
    denominators = [np.bincount(pair_code[dd <= rad], minlength=9) for rad in RADII]
    def signal(profile):
        ans = []
        for radius in RADII:
            near = dd <= radius
            for ligand, receptor in AXES[:2]:
                positive = near & profile[es, gene_index[ligand]] & profile[er, gene_index[receptor]]
                ans.append(np.bincount(pair_code[positive], minlength=9))
        return np.array(ans, dtype=np.int64).reshape(-1)
    obs = signal(values); lr_rows = []; lr_null = {}
    identity_codes = pd.Categorical(label[included], categories=SEND + RECV).codes
    profile_code = values.astype(np.int64) @ (2 ** np.arange(4, dtype=np.int64))
    for tile in TILES:
        groups = strata(xy[included], total[included], tile, identity_codes)
        variable = [g for g in groups if len(np.unique(profile_code[g])) > 1]
        rng = np.random.default_rng(2026110910 + tile)
        null = np.zeros((999, len(obs)), dtype=np.int64)
        for rep in range(999):
            shuffled = values.copy()
            for g in variable:
                shuffled[g] = values[rng.permutation(g)]
            null[rep] = signal(shuffled)
        lr_null[str(tile)] = null
        for ri, radius in enumerate(RADII):
            for ai, (ligand, receptor) in enumerate(AXES):
                for si, sender in enumerate(SEND):
                    for ti, receiver in enumerate(RECV):
                        pair = si * 3 + ti; n_edges = int(denominators[ri][pair])
                        measured = ligand in rn and receptor in rn
                        available = measured and n_edges > 0
                        j = (ri * 2 + ai) * 9 + pair if ai < 2 else None
                        ob = int(obs[j]) if measured else None
                        expected_null = float(null[:, j].mean()) if measured else None
                        row = {'sender': sender, 'receiver': receiver, 'ligand': ligand, 'receptor': receptor,
                               'radius_um': radius, 'tile_um': tile, 'spatial_edges': n_edges,
                               'components_measured': measured, 'available': available,
                               'status': 'MEASURED_LOCAL_EXPRESSION_TEST' if available else ('MISSING_PANEL_COMPONENT' if not measured else 'NO_GRAPH_EDGES'),
                               'observed_dualpositive_edges': ob, 'null_mean_dualpositive_edges': expected_null,
                               'observed_fraction': ob / n_edges if available else None,
                               'null_mean_fraction': expected_null / n_edges if available else None,
                               'difference_pp': 100 * (ob - expected_null) / n_edges if available else None,
                               'p_one_sided_conditional': float((1 + (null[:, j] >= ob).sum()) / 1000) if available else 1.,
                               'n_permutations': 999}
                        lr_rows.append(row)
    lr = pd.DataFrame(lr_rows); assert len(lr) == 108
    lr['BH_q_complete108'] = source.fdr(lr.p_one_sided_conditional)
    lr['BY_q_complete108'] = source.fdr(lr.p_one_sided_conditional, by=True)
    lr.to_csv(o / 'spatial_expression_communication_complete108.tsv', sep='\t', index=False, na_rep='NA')
    np.savez_compressed(o / 'communication_null_PRIVATE.npz', **lr_null)
    np.savez_compressed(o / 'communication_inputs_PRIVATE.npz', included=included, expression=values,
                        es=es, er=er, distance=dd, pair_code=pair_code)
    # Reuse the established spatial CellChat input contract without classification.
    unit = o / 'imaging/ccRCC_current_identity'; unit.mkdir(parents=True)
    x = xr[included].astype(np.float64)
    x = (diags(1e4 / total[included]) @ x).tocsr(); x.data = np.log1p(x.data)
    mmwrite(unit / 'normalized_expression_PRIVATE.mtx', x.T)
    pd.Series(rn).to_csv(unit / 'genes.tsv', index=False, header=False, sep='\t')
    pd.DataFrame({'broad': label[included], 'source_FOV': 'one_tissue',
                  'x_um': xy[included, 0], 'y_um': xy[included, 1]}, index=ids[included]).to_csv(unit / 'metadata_PRIVATE.tsv', sep='\t')
    rec = {'UTC': datetime.now(timezone.utc).isoformat(), 'status': 'EXISTING_IDENTITY_FORWARD_COMPUTED_CELLCHAT_INPUT_READY',
           'identity': expected, 'original_QC_keep': int(keep.sum()), 'CellChat_cells': len(included),
           'composition_hypotheses': 12, 'composition_BH_lt05': int((comp.BH_q_complete12 < .05).sum()),
           'LR_hypotheses': 108, 'LR_measured_slots': int(lr.components_measured.sum()),
           'LR_BH_lt05': int((lr.BH_q_complete108 < .05).sum()), 'LR_BY_lt05': int((lr.BY_q_complete108 < .05).sum()),
           'gene_panel': 477, 'identity_changed': False, 'reannotation_prerequisite': False,
           'no_patient_causal_or_receptor_activation_claim': True}
    atomic_json(o / 'complete.json', rec); print(json.dumps(rec), flush=True)


with exclusive_lock(r / 'current_identity_forward_v1.lock'):
    main()
