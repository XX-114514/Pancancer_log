#!/usr/bin/env python3
"""Sample the CT75 annotation scheme from the qualified TME atlas as a spatial reference.

Run inside wehi_oldpath_sandbox.sh (frozen JSONs name old-server paths).
Every node of the CT_identity_cluster scheme (all ten parents, immune and
non-immune, including reference-unresolved fallback nodes) inside the 439 CT75
work units is sampled (immune cells from the CT75 immune analysis domain,
non-immune cells from the same units): up to --n-ref cells for the reference and a disjoint
up to --n-heldout cells for annotation validation. Genes are the 12,783 genes
measured in every source. X holds the adopted raw integer counts (as in a
Visium/Xenium count matrix); layers['analysis'] subtracts the adopted ambient
estimate, the layer the atlas CellTypist labels were computed on.
obs['total_native'] is the full-gene analysis RNA total used for the atlas
CP10K, kept so held-out cells can be normalised either way.

Usage: build_atlas_reference_v1.py --v5 OLD_V5_PATH --output DIR [--n-ref 1000 --n-heldout 200 --seed 20261102]
"""
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad


def load(p):
    return json.loads(Path(p).read_text())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--v5', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--n-ref', type=int, default=1000)
    ap.add_argument('--n-heldout', type=int, default=200)
    ap.add_argument('--min-node-cells', type=int, default=30)
    ap.add_argument('--seed', type=int, default=20261102)
    a = ap.parse_args()
    V5 = Path(a.v5); runs = V5 / 'runs'; out = Path(a.output); out.mkdir(parents=True, exist_ok=False)
    J = runs / '20261026_joined_native_atlas_delivery_v1'
    sys.path.insert(0, str(J / 'scripts'))
    from native_archive_collection_v1 import NativeArchiveCollection
    from native_annotation_review_v1 import native_measured_genes
    from native_panel_aggregate_v1 import gather_csr_rows
    sys.path.insert(0, str(runs / '20261101_TME_CT75_native_program_LR_review_v2/scripts'))
    from audit_TME_standard_Harmony_v2 import gather_nonnegative_float_delta

    fit = runs / '20261101_TME_immune_CM_main_CT_identity_own_v1'
    acc = load(load(fit / 'config/spec_v1.json')['input_acceptance']['path'])
    m = load(acc['annotation_acceptance']['path'])
    meta = pd.read_csv(acc['metadata']['path'], sep='\t')
    fq = pd.read_csv(acc['frequency']['path'], sep='\t', index_col=0)
    prep = load(runs / '20261101_TME_immune_CM_preparation_v1/config/spec_v1.json')
    H = Path(prep['model'])
    axes = np.load(H / 'outputs/private_cell_axes_v1.npz')
    scope = load(H / 'outputs/scope_and_common_genes_v1.json')
    sources = scope['sources']; genes = list(scope['genes'])
    codes = np.load(m['node_codes']['path']); parent = np.load(m['parent_codes']['path'])
    unit = np.load(load(Path(prep['unit_closure']) / 'outputs/accepted_fit_qualification_unit_bridge_v3.json')['unit_codes']['path'])
    domain = np.load(Path(prep['design']) / 'outputs/private_technical_material_candidate_mask_v1.npy')
    totals = np.load(H / 'outputs/private_analysis_RNA_totals_v1.npy', mmap_mode='r')
    lineage = {x['source']: x for x in load(H / 'outputs/source_input_lineage_v1.json')}
    nodes = m['nodes']; parents = {r['code']: r['name'] for r in m['parents']}
    ct75 = set(fq.index)

    in_units = np.isin(unit, meta.qualification_unit_code.to_numpy()) & (codes >= 0)
    immune = np.isin(parent, [c for c, n in parents.items() if n in ('T_NK', 'B_Plasma', 'Myeloid', 'Mast')])
    # Immune cells: the CT75 immune analysis domain. Non-immune cells of the same
    # work units are added so spot deconvolution has tumour/stromal signatures.
    mask = (domain & in_units) | (in_units & ~immune)
    rng = np.random.default_rng(a.seed)
    take_ref, take_held = [], []
    node_rows = []
    for code in np.unique(codes[mask]):
        ix = np.flatnonzero(mask & (codes == code)); n = nodes[int(code)]
        row = dict(code=int(code), node=n['name'], parent=parents[n['parent_code']], numeric_cluster=n['numeric_cluster'],
                   identity=n['reference_label'], discovery_eligible=n['discovery_eligible'], in_CT75=n['name'] in ct75,
                   n_domain_cells=len(ix))
        if len(ix) >= a.min_node_cells:
            perm = rng.permutation(ix)
            r = perm[:a.n_ref]; h = perm[a.n_ref:a.n_ref + a.n_heldout]
            take_ref.append(r); take_held.append(h); row.update(n_ref=len(r), n_heldout=len(h))
        else:
            row.update(n_ref=0, n_heldout=0)
        node_rows.append(row)
    pd.DataFrame(node_rows).to_csv(out / 'reference_node_table_v1.tsv', sep='\t', index=False)
    ref = np.sort(np.concatenate(take_ref)); held = np.sort(np.concatenate(take_held))
    split = np.zeros(len(codes), np.int8); split[ref] = 1; split[held] = 2
    cells = np.sort(np.concatenate([ref, held]))
    print('nodes', len(node_rows), 'reference', len(ref), 'heldout', len(held), flush=True)

    native = NativeArchiveCollection(load(J / 'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'])
    blocks_raw, blocks_an, order = [], [], []
    for sc, source in enumerate(sources):
        ids = cells[axes['source_code'][cells] == sc]
        if not len(ids):
            continue
        names, physical = native_measured_genes(native, source); lookup = dict(zip(names, physical))
        cols = np.array([int(lookup[g]) for g in genes])
        z = native.store(source)
        for lo in range(0, len(ids), 1024):
            ix = ids[lo:lo + 1024]; nr = axes['native_row'][ix]
            x = gather_csr_rows(z[lineage[source]['counts_key']], nr, max_span=8192).astype(np.float64).tocsr()
            if 'ambient_removed_counts' in z['layers']:
                d = gather_nonnegative_float_delta(z['layers/ambient_removed_counts'], nr)
                xa = (x - d).tocsr(); xa.data[xa.data < 0] = 0; xa.eliminate_zeros()
            else:
                xa = x
            np.testing.assert_allclose(np.asarray(xa.sum(axis=1)).ravel(), totals[ix], rtol=1e-8, atol=1e-6)
            blocks_raw.append(x[:, cols].astype(np.float32)); blocks_an.append(xa[:, cols].astype(np.float32)); order.append(ix)
        native._fields.clear()
        print('read', source, len(ids), flush=True)
    order = np.concatenate(order); X = sp.vstack(blocks_raw).tocsr(); XA = sp.vstack(blocks_an).tocsr()
    um = meta.set_index('qualification_unit_code')
    obs = pd.DataFrame({
        'node': [nodes[int(c)]['name'] for c in codes[order]],
        'parent': [parents[int(p)] for p in parent[order]],
        'numeric_cluster': [nodes[int(c)]['numeric_cluster'] for c in codes[order]],
        'identity': [nodes[int(c)]['reference_label'] for c in codes[order]],
        'source': np.asarray(sources)[axes['source_code'][order]],
        'unit': unit[order],
        'cancer_type': um.loc[unit[order], 'cancer_type'].to_numpy(),
        'split': np.where(split[order] == 1, 'reference', 'heldout'),
        'total_native': np.asarray(totals[order], dtype=np.float64),
        'fit_row': order,
    }, index=['c' + str(i) for i in order])
    obs['parent_cluster'] = obs.parent + '|C' + obs.numeric_cluster.astype(str)
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes), layers={'analysis': XA})
    adata.uns['scheme'] = 'CT_identity_cluster'; adata.uns['seed'] = a.seed
    adata.write_h5ad(out / 'atlas_reference_v1.h5ad', compression='gzip')
    json.dump(dict(n_reference=int(len(ref)), n_heldout=int(len(held)), n_genes=len(genes), n_nodes_sampled=len(take_ref),
                   min_node_cells=a.min_node_cells, n_ref=a.n_ref, n_heldout_per_node=a.n_heldout, seed=a.seed,
                   counts='X adopted raw counts; layers.analysis ambient-subtracted'),
              open(out / 'reference_summary_v1.json', 'x'), indent=2)


if __name__ == '__main__':
    main()
