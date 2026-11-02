#!/usr/bin/env python3
"""cell2location deconvolution of GSE273952 ccRCC Visium with the CT75-scheme atlas reference.

Follows the CoVarNet paper (Shi et al., Nature 2025, Methods "Analysis of spatial
transcriptomics data"): reference subsampled to <=1,000 cells per subset,
N_cells_per_location=5, detection_alpha=20, otherwise cell2location defaults
(reference regression 250 epochs; spatial model 30,000 epochs, full batch).
Reference subsets are all nodes of the CT_identity_cluster scheme with >=30
cells in the CT75 work units, batch = atlas source; Visium batch = section.
With --section-index each section gets its own full-batch spatial model.

Usage: visium_cell2location_v1.py --reference REF.h5ad --visium-root RAW_DIR --output DIR [--stage ref|spatial|all]
"""
import argparse, glob, gzip, json
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


def read_visium(root, only=None):
    parts = []
    for gsm in sorted(Path(root).glob('GSM*')):
        if only is not None and gsm.name != only:
            continue
        h5 = glob.glob(str(gsm / '*_filtered_feature_bc_matrix.h5'))[0]
        sample = Path(h5).name.split('_filtered')[0]
        a = sc.read_10x_h5(h5); a.var_names_make_unique()
        a.var['gene_ids'] = a.var['gene_ids'].astype(str)
        pos = pd.read_csv(glob.glob(str(gsm / '*_tissue_positions_list.csv.gz'))[0], header=None,
                          names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row', 'pxl_col'])
        if str(pos.iloc[0, 0]) == 'barcode':
            pos = pos.iloc[1:]
        pos = pos.set_index('barcode').loc[a.obs_names]
        a = a[pos.in_tissue.astype(int).to_numpy() == 1].copy(); pos = pos.loc[a.obs_names]
        a.obs['sample'] = sample; a.obs['gsm'] = gsm.name
        for c in ['array_row', 'array_col', 'pxl_row', 'pxl_col']:
            a.obs[c] = pos[c].astype(float).to_numpy()
        a.obsm['spatial'] = a.obs[['pxl_col', 'pxl_row']].to_numpy()
        a.obs_names = sample + ':' + a.obs_names
        parts.append(a)
    return ad.concat(parts, join='inner', merge='same')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reference', required=True)
    ap.add_argument('--visium-root', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--stage', default='all', choices=['ref', 'spatial', 'all'])
    ap.add_argument('--seed', type=int, default=20261102)
    ap.add_argument('--section-index', type=int, default=None,
                    help='fit one section (0-based, sorted GSM order) with its own full-batch model; '
                         'used because all 19 sections in one full batch exceed a 40 GB GPU')
    a = ap.parse_args()
    import cell2location, scvi, torch
    from cell2location.models import RegressionModel
    from cell2location.utils.filtering import filter_genes
    scvi.settings.seed = a.seed
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    acc = 'gpu' if torch.cuda.is_available() else 'cpu'
    inf_path = out / 'reference_signatures_v1.tsv'
    if a.stage in ('ref', 'all'):
        ref = ad.read_h5ad(a.reference)
        ref = ref[ref.obs.split == 'reference'].copy()
        ref.X = ref.X.astype(np.float32)
        selected = filter_genes(ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
        ref = ref[:, selected].copy()
        RegressionModel.setup_anndata(ref, batch_key='source', labels_key='node')
        mod = RegressionModel(ref)
        mod.train(max_epochs=250, accelerator=acc)
        ref = mod.export_posterior(ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500})
        mod.save(str(out / 'reference_model'), overwrite=True)
        cols = [f'means_per_cluster_mu_fg_{i}' for i in ref.uns['mod']['factor_names']]
        inf = ref.varm['means_per_cluster_mu_fg'][cols].copy(); inf.columns = ref.uns['mod']['factor_names']
        inf.to_csv(inf_path, sep='\t')
        np.save(out / 'reference_training_history_v1.npy', np.asarray(mod.history['elbo_train']).ravel())
        print('reference done', inf.shape, flush=True)
    if a.stage in ('spatial', 'all'):
        inf = pd.read_csv(inf_path, sep='\t', index_col=0)
        only = None
        if a.section_index is not None:
            only = sorted(p.name for p in Path(a.visium_root).glob('GSM*'))[a.section_index]
            out = out / 'sections' / only; out.mkdir(parents=True, exist_ok=True)
        vis = read_visium(a.visium_root, only)
        vis.var['MT'] = vis.var_names.str.startswith('MT-')
        vis = vis[:, ~vis.var.MT.to_numpy()].copy()
        shared = np.intersect1d(vis.var_names, inf.index)
        vis = vis[:, shared].copy(); inf = inf.loc[shared]
        vis.X = vis.X.astype(np.float32)
        cell2location.models.Cell2location.setup_anndata(vis, batch_key='sample' if only is None else None)
        mod = cell2location.models.Cell2location(vis, cell_state_df=inf, N_cells_per_location=5, detection_alpha=20)
        mod.train(max_epochs=30000, batch_size=None, train_size=1, accelerator=acc)
        vis = mod.export_posterior(vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs})
        np.save(out / 'spatial_training_history_v1.npy', np.asarray(mod.history['elbo_train']).ravel())
        q05 = vis.obsm['q05_cell_abundance_w_sf'].copy(); q05.columns = vis.uns['mod']['factor_names']
        means = vis.obsm['means_cell_abundance_w_sf'].copy(); means.columns = vis.uns['mod']['factor_names']
        spots = vis.obs[['sample', 'gsm', 'array_row', 'array_col', 'pxl_row', 'pxl_col']].copy()
        spots['n_counts'] = np.asarray(vis.X.sum(axis=1)).ravel()
        spots.to_csv(out / 'visium_spots_v1.tsv.gz', sep='\t')
        q05.to_csv(out / 'visium_q05_abundance_v1.tsv.gz', sep='\t')
        means.to_csv(out / 'visium_mean_abundance_v1.tsv.gz', sep='\t')
        json.dump(dict(n_spots=int(vis.n_obs), n_genes=int(vis.n_vars), n_sections=int(vis.obs['sample'].nunique()),
                       n_subsets=int(inf.shape[1]), N_cells_per_location=5, detection_alpha=20, max_epochs=30000,
                       cell2location=cell2location.__version__, scvi_tools=scvi.__version__, seed=a.seed),
                  open(out / 'visium_cell2location_summary_v1.json', 'w'), indent=2)


if __name__ == '__main__':
    main()
