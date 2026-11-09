"""Single-donor ccRCC: identity review and exposure-independent vascular tests.

All row-level objects stay under the supplied private run root. Public outputs
are aggregate tables/figures only. ACKR1 is never an identity/clustering feature.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys
import hashlib
import numpy as np
import pandas as pd
import h5py
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from resumable_public_zip_v1 import atomic_json, sha256, exclusive_lock


PANELS = {
    'B': ['CD19', 'MS4A1', 'CD79A'],
    'T': ['CD3D', 'CD3E', 'TRAC'],
    'NK_cytotoxic': ['KLRD1', 'NKG7', 'GNLY', 'NCAM1', 'PRF1'],
    'EC_identity': ['PECAM1', 'VWF', 'EMCN', 'ERG'],
    'EC_exposure_review': ['ACKR1', 'CA4', 'SOX17'],
    'Lymphatic': ['PROX1', 'LYVE1', 'PDPN'],
    'Mural': ['PDGFRB', 'ACTA2', 'MYH11', 'CSPG4', 'MCAM', 'NOTCH3'],
    'Fibroblast': ['PDGFRA', 'DCN', 'LUM', 'COL1A1', 'COL1A2', 'COL3A1'],
    'Myeloid': ['LST1', 'TYROBP', 'CD68', 'CD163', 'S100A8'],
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'CA9'],
    'Chemokine_review': ['CCL19', 'CCL21', 'CXCL12', 'CXCL13', 'CCR7', 'CXCR4'],
}
PROTEIN_EXPECTED = {'CD20': 'B', 'CD3E': 'T', 'CD4': 'T/myeloid', 'CD8A': 'T/NK',
                    'CD31': 'EC', 'alphaSMA': 'Mural', 'CD68': 'Myeloid',
                    'CD163': 'Myeloid', 'CD45': 'Immune', 'PanCK': 'Epithelial',
                    'GranzymeB': 'Cytotoxic', 'CD16': 'NK/myeloid'}
TYPES = ['B', 'T', 'NK']
RADII = [25, 50, 100, 200]


def stamp():
    return datetime.now(timezone.utc).isoformat()


def load(root):
    rec = json.loads((root / 'inputs/complete_PRIVATE.json').read_text())
    files = {}
    for item in rec['members']:
        p = root / 'inputs' / item['local_name']
        assert p.stat().st_size == item['bytes'] and sha256(p) == item['sha256']
        files[item['member']] = p
    with h5py.File(files['cell_feature_matrix.h5']) as f:
        m = f['matrix']
        x = csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]),
                       shape=tuple(m['shape'][:])).T.tocsr()
        names = np.array([v.decode() for v in m['features/name'][:]])
        ft = np.array([v.decode() for v in m['features/feature_type'][:]])
        ids = np.array([v.decode() for v in m['barcodes'][:]])
        scale = float(f.attrs['protein_scaling_factor'])
        metric = f.attrs['protein_metric']
        if isinstance(metric, bytes):
            metric = metric.decode()
    assert scale > 0 and metric == 'MEAN_PER_CELL_STAIN'
    rna = ft == 'Gene Expression'
    protein = ft == 'Protein Expression'
    assert rna.sum() == 477 and protein.sum() == 27
    rn, pn = names[rna], names[protein]
    assert len(set(rn)) == len(rn) and len(set(pn)) == len(pn)
    xr = x[:, rna]
    xp = x[:, protein].toarray().astype(float) / scale
    c = pd.read_parquet(files['cells.parquet'])
    assert c.cell_id.is_unique and set(c.cell_id) == set(ids)
    c = c.set_index('cell_id').loc[ids]
    xy = c[['x_centroid', 'y_centroid']].to_numpy()
    assert np.isfinite(xy).all()
    total = np.asarray(xr.sum(1)).ravel()
    ng = xr.getnnz(1)
    needed = sorted(set(sum(PANELS.values(), [])))
    a = {g: xr[:, list(rn).index(g)].toarray().ravel() for g in needed if g in rn}
    return ids, c, xy, xr, rn, xp, pn, scale, total, ng, a


def identities(a, total, ng):
    n = len(total)
    def pos(g):
        return a.get(g, np.zeros(n)) >= 1
    def count(gg):
        return np.sum([pos(g) for g in gg], axis=0)
    b = count(PANELS['B']) >= 2
    t = count(PANELS['T']) >= 2
    ec = count(['PECAM1', 'VWF', 'EMCN']) >= 2
    ly = count(PANELS['Lymphatic']) >= 2
    nk = pos('KLRD1') & (pos('NKG7') | pos('GNLY')) & (count(PANELS['T']) == 0) & ~b & ~ec
    mural = (count(['PDGFRB', 'ACTA2', 'MYH11']) >= 2) & ~ec & ~b & ~t
    masks = {'B': b, 'T': t, 'NK': nk, 'Blood_EC': ec & ~ly,
             'Lymphatic': ec & ly, 'Mural': mural}
    conflict = np.sum(list(masks.values()), axis=0) > 1
    depth = (total >= 50) & (ng >= 20)
    keep = depth & ~conflict
    label = np.full(n, 'Other', dtype=object)
    for k, v in masks.items():
        label[v & ~conflict] = k
    return pos, count, masks, conflict, depth, keep, label


def marker_table(a, groups, names):
    rows = []
    for group in names:
        z = groups == group
        for expected, gg in PANELS.items():
            for gene in gg:
                rows.append({'group': group, 'n_cells': int(z.sum()), 'marker': gene,
                             'expected_identity': expected, 'available': gene in a,
                             'fraction_positive': float((a[gene][z] >= 1).mean()) if gene in a and z.any() else np.nan,
                             'mean_log1p_RNA_UMI': float(np.log1p(a[gene][z]).mean()) if gene in a and z.any() else np.nan})
    return pd.DataFrame(rows)


def protein_table(xp, pn, groups, names):
    rows = []
    for group in names:
        z = groups == group
        if not z.any():
            continue
        for j, marker in enumerate(pn):
            q = np.quantile(xp[z, j], [.25, .5, .75])
            rows.append({'group': group, 'n_cells': int(z.sum()), 'antibody': marker,
                         'expected_identity': PROTEIN_EXPECTED.get(marker, 'review'),
                         'q25_mean_intensity': float(q[0]), 'median_mean_intensity': float(q[1]),
                         'q75_mean_intensity': float(q[2]), 'units': 'unscaled_arbitrary_mean_intensity'})
    return pd.DataFrame(rows)


def marker_plot(table, path, title):
    import matplotlib.pyplot as plt
    groups = list(dict.fromkeys(table.group))
    genes = list(dict.fromkeys(zip(table.marker, table.expected_identity)))
    fig, ax = plt.subplots(figsize=(max(13, len(genes) * .23), max(4, len(groups) * .32 + 3)))
    for y, group in enumerate(groups):
        for x, (gene, expected) in enumerate(genes):
            row = table[(table.group == group) & (table.marker == gene) & (table.expected_identity == expected)].iloc[0]
            if row.available and np.isfinite(row.fraction_positive):
                dots = ax.scatter(x, y, s=100 * row.fraction_positive,
                                  c=[row.mean_log1p_RNA_UMI], vmin=0, vmax=2.5, cmap='viridis')
            else:
                ax.text(x, y, '×', color='.6', ha='center', va='center', fontsize=7)
    ax.set_xticks(range(len(genes)), [g + ' | ' + e for g, e in genes], rotation=90, fontsize=6)
    ax.set_yticks(range(len(groups)), groups, fontsize=7)
    ax.set_title(title + '\nSize = RNA detection fraction; × = not measured')
    fig.colorbar(dots, ax=ax, label='Mean log1p native RNA UMI', fraction=.015, pad=.01)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def protein_plot(table, path):
    import matplotlib.pyplot as plt
    mat = table.pivot(index='antibody', columns='group', values='median_mean_intensity')
    log = np.log1p(mat.to_numpy())
    z = (log - log.mean(1, keepdims=True)) / np.maximum(log.std(1, keepdims=True), 1e-8)
    fig, ax = plt.subplots(figsize=(9, 10))
    im = ax.imshow(z, aspect='auto', cmap='coolwarm', vmin=-2, vmax=2)
    ax.set_yticks(range(len(mat)), [g + ' | ' + PROTEIN_EXPECTED.get(g, 'review') for g in mat.index], fontsize=8)
    ax.set_xticks(range(len(mat.columns)), mat.columns, rotation=55, ha='right')
    ax.set_title('RNA candidate identity / protein cross-review\nColour standardized within each antibody only')
    fig.colorbar(im, ax=ax, label='Within-antibody z of log1p group median')
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def fdr(p, by=False):
    p = np.array(p, float)
    order = np.argsort(p, kind='stable')
    adjust = len(p) / np.arange(1, len(p) + 1)
    if by:
        adjust *= np.sum(1 / np.arange(1, len(p) + 1))
    q = np.minimum.accumulate((p[order] * adjust)[::-1])[::-1]
    out = np.empty(len(p)); out[order] = np.minimum(q, 1)
    return out


def log_odds(y1, d1, y0, d0):
    return np.log((y1 + .5) / (d1 - y1 + .5)) - np.log((y0 + .5) / (d0 - y0 + .5))


def prepare_boundaries(root, ids, ec_global, xy, query_global, out):
    import shapely
    from shapely import STRtree
    receipt = json.loads((root / 'inputs/boundary_complete_PRIVATE.json').read_text())
    p = root / 'inputs' / receipt['local_name']
    assert p.stat().st_size == receipt['bytes'] and sha256(p) == receipt['sha256']
    b = pd.read_parquet(p)
    required = ['cell_id', 'vertex_x', 'vertex_y']
    assert all(col in b for col in required)
    target = pd.Index(ids[ec_global])
    b = b[b.cell_id.isin(target)]
    assert np.isfinite(b[['vertex_x', 'vertex_y']].to_numpy()).all()
    polygons, target_indices, invalid = [], [], 0
    for cell_id, g in b.groupby('cell_id', sort=False):
        polygon = shapely.Polygon(g[['vertex_x', 'vertex_y']].to_numpy())
        if not polygon.is_valid or polygon.is_empty or polygon.area <= 0:
            invalid += 1
            continue
        polygons.append(polygon)
        target_indices.append(target.get_loc(cell_id))
    polygons = np.array(polygons, dtype=object)
    target_indices = np.array(target_indices)
    assert len(polygons) > .95 * len(ec_global), 'Incomplete/invalid EC polygons; do not silently change geometry'
    ec_centers = shapely.points(xy[ec_global[target_indices]])
    center_offsets = shapely.distance(ec_centers, polygons)
    assert np.quantile(center_offsets, .99) < 10, 'Coordinate alignment requires review'
    tree = STRtree(polygons)
    points = shapely.points(xy[query_global])
    pair, distances = tree.query_nearest(points, all_matches=True, return_distance=True)
    hits = np.bincount(pair[0], minlength=len(points))
    single = hits[pair[0]] == 1
    # Exclude all ties and all zero-distance interiors, independent of ACKR1.
    accepted = single & (distances > 1e-8)
    q = pair[0, accepted]
    owner = np.full(len(points), -1, dtype=int)
    distance = np.full(len(points), np.nan)
    owner[q] = target_indices[pair[1, accepted]]
    distance[q] = distances[accepted]
    rng = np.random.default_rng(301)
    sample = rng.choice(q, min(64, len(q)), replace=False)
    brute_errors = []
    for i in sample:
        brute = shapely.distance(points[i], polygons).min()
        brute_errors.append(abs(float(brute) - distance[i]))
    assert max(brute_errors) < 1e-8
    audit = {'EC_polygons': len(polygons), 'EC_without_valid_polygon': len(ec_global) - len(polygons),
             'invalid_polygons': invalid, 'query_cells': len(points),
             'query_nonunique_nearest_excluded': int((hits > 1).sum()),
             'query_zero_distance_any_nearest': int(len(np.unique(pair[0, distances <= 1e-8]))),
             'query_accepted_unique_positive_distance': int(np.isfinite(distance).sum()),
             'EC_centroid_to_own_polygon_p99_um': float(np.quantile(center_offsets, .99)),
             'brute_force_checks': len(sample), 'max_abs_distance_error_um': max(brute_errors),
             'meaning': 'Distance to segmented EC cell region, not reconstructed complete vessel wall',
             'boundary_source_sha256': receipt['sha256']}
    atomic_json(out / 'boundary_audit.json', audit)
    return distance, owner, (polygons, target_indices), audit


def spatial_tests(root, xy, total, label, keep, pos, ids, out):
    ec_global = np.flatnonzero(keep & (label == 'Blood_EC'))
    query_global = np.flatnonzero(keep & ~np.isin(label, ['Blood_EC', 'Lymphatic']))
    e_xy, q_xy = xy[ec_global], xy[query_global]
    target = pos('ACKR1')[ec_global]
    assert target.any() and (~target).any()
    cd, co = cKDTree(e_xy).query(q_xy)
    sample = np.random.default_rng(302).choice(len(q_xy), min(128, len(q_xy)), replace=False)
    assert np.allclose(cd[sample], cdist(q_xy[sample], e_xy).min(1), atol=1e-8, rtol=0)
    pdist, po, polygons, ba = prepare_boundaries(root, ids, ec_global, xy, query_global, out)
    np.savez_compressed(out / 'geometry_PRIVATE.npz', ec_global=ec_global, query_global=query_global,
                        centroid_distance=cd, centroid_owner=co, polygon_distance=pdist, polygon_owner=po,
                        ACKR1_detected=target)
    qlabels = label[query_global]
    endpoints, y_vectors, d_vectors, descriptions = [], [], [], []
    for geometry, distance, owner in [('centroid', cd, co), ('EC_cell_polygon', pdist, po)]:
        for radius in RADII:
            near = np.isfinite(distance) & (distance <= radius) & (owner >= 0)
            for typ in TYPES:
                y = near & (qlabels == typ)
                for denominator in ['allQC_nonEC', 'B_T_NK_only']:
                    d = near if denominator == 'allQC_nonEC' else near & np.isin(qlabels, TYPES)
                    yv = np.bincount(owner[y], minlength=len(ec_global)).astype(float)
                    dv = np.bincount(owner[d], minlength=len(ec_global)).astype(float)
                    y1, d1, y0, d0 = yv[target].sum(), dv[target].sum(), yv[~target].sum(), dv[~target].sum()
                    available = min(d1, d0) > 0 and yv.sum() > 0
                    endpoints.append({'geometry': geometry, 'radius_um': radius, 'receiver': typ,
                                      'denominator': denominator, 'available': bool(available)})
                    y_vectors.append(yv); d_vectors.append(dv)
                    descriptions.append({'target_y': int(y1), 'target_denominator': int(d1),
                                         'control_y': int(y0), 'control_denominator': int(d0),
                                         'target_fraction': y1 / d1 if d1 else None,
                                         'control_fraction': y0 / d0 if d0 else None,
                                         'difference_pp': 100 * (y1 / d1 - y0 / d0) if min(d1, d0) else None,
                                         'log_odds': float(log_odds(y1, d1, y0, d0)) if available else None})
    Y, D = np.array(y_vectors), np.array(d_vectors)
    np.savez_compressed(out / 'endpoint_weights_PRIVATE.npz', Y=Y, D=D, target=target,
                        e_xy=e_xy, ec_rna_total=total[ec_global])
    ytotal, dtotal = Y.sum(1), D.sum(1)
    obs = log_odds(Y @ target, D @ target, ytotal - Y @ target, dtotal - D @ target)
    rows, strata_reviews = [], []
    null_arrays = {}
    for tile in [250, 500]:
        keys = np.column_stack([np.floor(e_xy / tile).astype(int),
                                np.floor(np.log2(total[ec_global])).astype(int)])
        _, inv = np.unique(keys, axis=0, return_inverse=True)
        order = np.argsort(inv, kind='stable')
        groups = np.split(order, np.flatnonzero(np.diff(inv[order])) + 1)
        movable = [g for g in groups if target[g].any() and (~target[g]).any()]
        strata_reviews.append({'tile_um': tile, 'depth_bin': 'floor_log2_native_EC_RNA_UMI',
                                'strata': len(groups), 'mixed_strata': len(movable),
                                'EC_in_mixed_strata': int(sum(map(len, movable))),
                                'EC_total': len(target), 'target_EC': int(target.sum())})
        rng = np.random.default_rng(20261109 + tile)
        permutations = np.tile(target, (999, 1)).astype(float)
        for g in movable:
            for k in range(999):
                permutations[k, g] = rng.permutation(target[g])
        # All endpoints share the same permutation realization, not independent tests.
        assert np.all(permutations.sum(1) == target.sum())
        for g in groups:
            assert np.all(permutations[:, g].sum(1) == target[g].sum())
        y1, d1 = permutations @ Y.T, permutations @ D.T
        null = log_odds(y1, d1, ytotal - y1, dtotal - d1)
        null_arrays[str(tile)] = null
        for i, (endpoint, desc) in enumerate(zip(endpoints, descriptions)):
            p = float((1 + (null[:, i] >= obs[i]).sum()) / 1000) if endpoint['available'] else 1.
            rows.append({**endpoint, **desc, 'tile_um': tile, 'n_permutations': 999,
                         'null_mean_log_odds': float(null[:, i].mean()),
                         'excess_over_null_log_odds': float(obs[i] - null[:, i].mean()),
                         'p_one_sided_conditional': p})
        atomic_json(out / ('progress_null_' + str(tile) + '.json'), {'UTC': stamp(), 'complete': True})
    assert len(rows) == 96
    table = pd.DataFrame(rows)
    table['BH_q_complete96'] = fdr(table.p_one_sided_conditional)
    table['BY_q_complete96'] = fdr(table.p_one_sided_conditional, by=True)
    table.to_csv(out / 'spatial_complete96.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(strata_reviews).to_csv(out / 'local_null_strata_review.tsv', sep='\t', index=False)
    np.savez_compressed(out / 'null_logodds_PRIVATE.npz', **null_arrays)
    depth = []
    for v in [True, False]:
        z = ec_global[target == v]
        depth.append({'ACKR1_detected': v, 'EC_cells': len(z),
                      'median_RNA_UMI': float(np.median(total[z])),
                      'RNA_UMI_q25': float(np.quantile(total[z], .25)),
                      'RNA_UMI_q75': float(np.quantile(total[z], .75))})
    pd.DataFrame(depth).to_csv(out / 'ACKR1_EC_depth_review.tsv', sep='\t', index=False)
    # Spatial display selected without exposure or effect: whole tissue, then most
    # EC-rich and most B/T/NK-rich 500um tiles (not independent statistical units).
    import matplotlib.pyplot as plt
    keys = np.floor(xy[keep] / 500).astype(int)
    tile = pd.DataFrame(keys, columns=['x', 'y'])
    tile['EC'] = (label[keep] == 'Blood_EC').astype(int)
    tile['lymph'] = np.isin(label[keep], TYPES).astype(int)
    bytile = tile.groupby(['x', 'y'])[['EC', 'lymph']].sum()
    rois = [tuple(bytile.EC.idxmax()), tuple(bytile.lymph.idxmax())]
    fig, axes = plt.subplots(1, 3, figsize=(17, 6))
    colors = {'B': '#2878b5', 'T': '#6d6d6d', 'NK': '#2ca25f'}
    for k, ax in enumerate(axes):
        region = keep.copy()
        if k:
            rx, ry = rois[k - 1]
            region &= (xy[:, 0] >= rx * 500) & (xy[:, 0] < (rx + 1) * 500) & (xy[:, 1] >= ry * 500) & (xy[:, 1] < (ry + 1) * 500)
        back = np.flatnonzero(region)[::(15 if k == 0 else 1)]
        ax.scatter(xy[back, 0], xy[back, 1], s=.25, c='#dddddd', rasterized=True)
        for exposure, color in [(False, '#e6ac65'), (True, '#c52c35')]:
            z = region & (label == 'Blood_EC') & (pos('ACKR1') == exposure)
            ax.scatter(xy[z, 0], xy[z, 1], s=.6 if k == 0 else 6, c=color,
                       label='EC ACKR1 ' + ('detected' if exposure else 'undetected'), rasterized=True)
        for typ, color in colors.items():
            z = region & (label == typ)
            ax.scatter(xy[z, 0], xy[z, 1], s=.5 if k == 0 else 5, c=color, label=typ + ' candidate', rasterized=True)
        if k:
            import shapely
            pol, ix = polygons
            from matplotlib.collections import LineCollection
            local = region[ec_global[ix]]
            ax.add_collection(LineCollection([np.asarray(g.exterior.coords) for g in pol[local]], colors='#992222', linewidths=.35))
        ax.set_aspect('equal'); ax.invert_yaxis()
        ax.set_xlabel('Native x (µm)'); ax.set_ylabel('Native y (µm)')
        ax.set_title(['Whole source', 'EC-rich tile; exposure-independent selection', 'Lymphocyte-rich tile; exposure-independent selection'][k], fontsize=9)
    axes[0].legend(fontsize=6, markerscale=4)
    fig.suptitle('ccRCC one donor: RNA candidates, not vessel/TLS ground truth')
    fig.tight_layout(); fig.savefig(out / 'spatial_identity_context.png', dpi=180); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, denominator in zip(axes, ['allQC_nonEC', 'B_T_NK_only']):
        for typ, color in colors.items():
            for geom, style in [('centroid', '-'), ('EC_cell_polygon', '--')]:
                z = table[(table.denominator == denominator) & (table.receiver == typ) & (table.geometry == geom) & (table.tile_um == 250)]
                ax.plot(z.radius_um, z.difference_pp, style, color=color, marker='o', label=typ + ' / ' + geom)
        ax.axhline(0, color='.7', lw=.8); ax.set_xscale('log', base=2)
        ax.set_xticks(RADII, RADII); ax.set_xlabel('Radius (µm)'); ax.set_title(denominator)
    axes[0].set_ylabel('ACKR1 detected − undetected EC neighborhood (pp)')
    axes[1].legend(fontsize=6)
    fig.suptitle('Descriptive composition differences; conditional tests in complete96 table')
    fig.tight_layout(); fig.savefig(out / 'geometry_radius_composition.png', dpi=180); plt.close(fig)
    return {'EC_cells': len(ec_global), 'ACKR1_detected_EC': int(target.sum()),
            'complete_hypotheses': 96, 'BH_q_lt05': int((table.BH_q_complete96 < .05).sum()),
            'BY_q_lt05': int((table.BY_q_complete96 < .05).sum()),
            'centroid_bruteforce_checks': len(sample), 'boundary': ba}


def diagnostic_clusters(xr, rn, xp, pn, a, count, depth, conflict, keep, label, out):
    import scanpy as sc
    import anndata as ad
    from sklearn.metrics import adjusted_rand_score
    immune = (count(PANELS['B'] + PANELS['T'] + PANELS['NK_cytotoxic'] + PANELS['Myeloid']) >= 2)
    vascular = (count(['PECAM1', 'VWF', 'EMCN'] + PANELS['Lymphatic'] + PANELS['Mural'] + PANELS['Fibroblast']) >= 2)
    features = rn != 'ACKR1'
    results = []
    for pool_name, pool in [('immune', immune), ('vascular_stromal', vascular)]:
        ix = np.flatnonzero(depth & pool)
        if len(ix) < 100:
            continue
        folder = out / ('clusters_' + pool_name)
        folder.mkdir(exist_ok=False)
        obj = ad.AnnData(xr[ix][:, features].astype(np.float32))
        obj.var_names = rn[features]
        obj.obs['identity_rule_candidate'] = label[ix]
        obj.obs['conflicting_programs'] = conflict[ix]
        obj.obs['main_QC_keep'] = keep[ix]
        sc.pp.normalize_total(obj, target_sum=1e4)
        sc.pp.log1p(obj)
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
        # Only diagnostic numeric cluster IDs, never automatic publication labels.
        markers = marker_table({g: v[ix] for g, v in a.items()}, group, groups)
        markers.to_csv(folder / 'marker_expected_identity.tsv', sep='\t', index=False, na_rep='NA')
        proteins = protein_table(xp[ix], pn, group, groups)
        proteins.to_csv(folder / 'protein_within_antibody.tsv', sep='\t', index=False, na_rep='NA')
        marker_plot(markers, folder / 'marker_expected_identity.png', pool_name + ' RNA-only clusters (seed 0)')
        protein_plot(proteins, folder / 'protein_within_antibody.png')
        means = np.asarray(obj.X.mean(0)).ravel()
        review = []
        for g in groups:
            z = group == g
            gm = np.asarray(obj.X[z].mean(0)).ravel()
            # Mean log-normalized difference for review, not DE significance.
            rm = (means * len(group) - gm * z.sum()) / (len(group) - z.sum())
            top = obj.var_names[np.argsort(gm - rm)[-8:][::-1]].tolist()
            row = {'cluster': g, 'n_cells': int(z.sum()), 'conflicting_program_cells': int(conflict[ix][z].sum()),
                   'main_QC_keep_cells': int(keep[ix][z].sum()), 'top_mean_logexpr_difference_markers': ';'.join(top)}
            for name in ['B', 'T', 'NK', 'Blood_EC', 'Lymphatic', 'Mural', 'Other']:
                row['candidate_' + name] = int((label[ix][z] == name).sum())
            review.append(row)
        pd.DataFrame(review).to_csv(folder / 'cluster_review.tsv', sep='\t', index=False)
        obj.obs['source_index_PRIVATE'] = ix
        obj.write_h5ad(folder / 'RNA_only_cluster_object_PRIVATE.h5ad', compression='gzip')
        results.append({'pool': pool_name, 'cells': len(ix), 'RNA_features': obj.n_vars,
                        'seed0_clusters': len(groups), 'ARI_min': min(s['ARI'] for s in stability),
                        'seed_stability': stability, 'ACKR1_used_for_clustering': False,
                        'protein_or_position_used_for_clustering': False})
        atomic_json(folder / 'complete.json', results[-1])
    return results


def main(root):
    out = root / 'followup_v1'
    out.mkdir(exist_ok=False)
    atomic_json(out / 'started.json', {'UTC': stamp(), 'status': 'running'})
    ids, c, xy, xr, rn, xp, pn, scale, total, ng, a = load(root)
    pos, count, masks, conflict, depth, keep, label = identities(a, total, ng)
    groups = np.where(keep, label, 'Excluded_QC_or_conflict')
    names = ['B', 'T', 'NK', 'Blood_EC', 'Lymphatic', 'Mural', 'Other', 'Excluded_QC_or_conflict']
    mt = marker_table(a, groups, names)
    pt = protein_table(xp, pn, groups, names)
    mt.to_csv(out / 'marker_expected_identity.tsv', sep='\t', index=False, na_rep='NA')
    pt.to_csv(out / 'protein_within_antibody.tsv', sep='\t', index=False, na_rep='NA')
    import matplotlib
    matplotlib.use('Agg')
    marker_plot(mt, out / 'marker_expected_identity.png', 'ccRCC ACKR1-independent identity review')
    protein_plot(pt, out / 'protein_within_antibody.png')
    counts = {'input_cells': len(ids), 'RNA_depth_QC_pass_before_conflicts': int(depth.sum()),
              'depth_QC_conflicts': int((depth & conflict).sum()), 'QC_keep': int(keep.sum()),
              'candidate_counts': {k: int((keep & (label == k)).sum()) for k in names[:-1]},
              'protein_scaling_factor_removed': scale, 'RNA_features': len(rn), 'protein_features': len(pn),
              'EC_gate': '>=2 of PECAM1/VWF/EMCN; ACKR1 excluded',
              'RNA_gene_names_unique': len(set(rn)) == len(rn),
              'RNA_protein_shared_names': sorted(set(rn) & set(pn)),
              'native_centroid_units': 'micrometre; no image-pixel scaling applied'}
    atomic_json(out / 'identity_summary.json', counts)
    c['identity_candidate'] = label; c['QC_keep'] = keep; c['RNA_total'] = total; c['RNA_genes'] = ng
    c['conflicting_programs'] = conflict
    c.to_parquet(out / 'identity_candidates_PRIVATE.parquet')
    spatial = spatial_tests(root, xy, total, label, keep, pos, ids, out)
    atomic_json(out / 'spatial_complete.json', {'UTC': stamp(), **spatial})
    clusters = diagnostic_clusters(xr, rn, xp, pn, a, count, depth, conflict, keep, label, out)
    receipt = {'UTC': stamp(), 'status': 'computed_not_frozen_identity_or_patient_validation',
               'identity': counts, 'spatial': spatial, 'diagnostic_clusters': clusters,
               'limits': ['one donor/section', 'RNA rule candidates, not final fine subtypes',
                          '477-gene panel; IGHD/XCL1/NCR1/CCL21 not measured',
                          'no independent pathology/TLS truth', 'segmented EC region, not complete vessel wall',
                          'conditional permutation significance is not patient replication'],
               'communications_run': False, 'clinical_models_run': False}
    atomic_json(out / 'complete.json', receipt)
    print(json.dumps(receipt), flush=True)


if __name__ == '__main__':
    root = Path(sys.argv[1]).resolve()
    with exclusive_lock(root / 'followup_writer.lock'):
        main(root)
