# Harmonise QC-passed supplementary shards (v2/v3 batches) to the obs schema of
# pancancer_analysis_ready_v1 (same rules as the supplement branch of
# scripts/source_ingest/build_analysis_zarr_v1.py) so they can go through decontX and integration v2.
# Usage: harmonise_supplement_v1.py REPO OUT_DIR BATCH_LABEL QC_DIR [QC_DIR ...]
import sys, os, glob, re
import numpy as np, pandas as pd, anndata as ad

repo, out, batch = sys.argv[1:4]
dirs = sys.argv[4:]
os.makedirs(out, exist_ok=True)
cmap = pd.concat([pd.read_csv(f'{repo}/inventories/cancer_type/{f}', sep='\t') for f in
                  ('cancer_type_map_v1.tsv', 'cancer_type_map_supplement_v1.tsv', 'cancer_type_map_supplement_v2.tsv')]).drop_duplicates('source', keep='last').set_index('source')
ONCO_MAIN = {'Hepatocellular carcinoma': 'HCC', 'Intrahepatic cholangiocarcinoma': 'IHCH', 'Endometrial cancer': 'UCEC', 'Ovarian cancer': 'OVARY'}
SITE = {'primary_organ_site': 'primary_organ_site', 'metastatic_site': 'metastatic_site', 'tumour_thrombus': 'tumour_thrombus'}
QCFLAG = {'GSE277524': 'floor_cell_calling_on_unfiltered_matrix'}
COLS = ['cell_id', 'source', 'data_batch', 'unit_id', 'sample_id', 'patient_id', 'patient_id_rule',
        'analysis_main', 'analysis_fine', 'oncotree_code', 'cancer_label_resolution',
        'sampling_site', 'material_detail', 'treatment_status', 'treatment_detail', 'treatment_rule',
        'fresh_frozen', 'fresh_frozen_detail', 'technology', 'suspension', 'capture_scope',
        'origin_store', 'origin_counts_key', 'origin_native_row', 'origin_cell_name',
        'ambient_correction', 'doublet_method', 'author_cell_type', 'qc_flag', 'total_counts', 'n_genes', 'pct_mt']


def fresh(v):
    v = str(v)
    f, c = bool(re.search(r'fresh', v, re.I)), bool(re.search(r'cryopreserv|frozen|thaw', v, re.I))
    return 'fresh_or_cryopreserved' if f and c else 'cryopreserved_or_frozen' if c else 'fresh' if f else 'unknown'


for d in dirs:
    for f in sorted(glob.glob(f'{d}/*.h5ad')):
        src = os.path.basename(f)[:-5]
        if src not in cmap.index:
            print('skip (not in cancer map / excluded):', src); continue
        a = ad.read_h5ad(f); o = a.obs; c = cmap.loc[src]
        h = pd.DataFrame(index=o.index)
        h['cell_id'] = o['cell_id'].astype(str); h['source'] = src; h['data_batch'] = batch
        h['unit_id'] = src + '|' + o['sample_id'].astype(str); h['sample_id'] = o['sample_id'].astype(str)
        h['patient_id'] = o['patient_id'].astype(str); h['patient_id_rule'] = 'GEO_sample_record'
        main = o['analysis_main'].astype(str) if 'analysis_main' in o else pd.Series(c['analysis_main'], index=o.index)
        h['analysis_main'] = main
        h['analysis_fine'] = o['analysis_fine'].astype(str) if 'analysis_fine' in o else main.map(
            lambda m: {'Hepatocellular carcinoma': 'HCC', 'Intrahepatic cholangiocarcinoma': 'ICC'}.get(m, c['analysis_fine']))
        h['oncotree_code'] = main.map(lambda m: ONCO_MAIN.get(m, c['oncotree_code'])) if c['resolution'] == 'per_sample' else c['oncotree_code']
        h['cancer_label_resolution'] = c['resolution']
        h['sampling_site'] = o['sampling_site'].astype(str).map(lambda s: SITE.get(s, 'unspecified'))
        h['material_detail'] = o['material'].astype(str)
        h['treatment_status'] = o['treatment_status'].astype(str); h['treatment_detail'] = o['treatment_detail'].astype(str)
        h['treatment_rule'] = 'GEO_sample_record_or_paper'
        h['fresh_frozen_detail'] = o['fresh_frozen'].astype(str); h['fresh_frozen'] = h['fresh_frozen_detail'].map(fresh)
        h['technology'] = o['technology'].astype(str); h['suspension'] = o['suspension'].astype(str); h['capture_scope'] = o['capture_scope'].astype(str)
        h['origin_store'] = f'GEO:{src.split("_")[0]}'; h['origin_counts_key'] = 'supplementary_file'; h['origin_native_row'] = -1
        h['origin_cell_name'] = o['origin_barcode'].astype(str)
        h['ambient_correction'] = o['ambient_correction'].astype(str); h['doublet_method'] = o['doublet_method'].astype(str)
        h['author_cell_type'] = o['author_cell_type'].astype(str) if 'author_cell_type' in o else ''
        h['qc_flag'] = QCFLAG.get(src, 'qc_v1_pass')
        X = a.X.tocsr(); mt = a.var_names.str.upper().str.startswith('MT-')
        h['total_counts'] = np.asarray(X.sum(1)).ravel().astype(np.int64); h['n_genes'] = np.diff(X.indptr).astype(np.int32)
        h['pct_mt'] = (np.asarray(X[:, mt].sum(1)).ravel() / np.maximum(h['total_counts'].values, 1) * 100).astype(np.float32)
        h = h[COLS]; h.index = h['cell_id'].values
        ad.AnnData(X=X, obs=h, var=pd.DataFrame(index=a.var_names.astype(str))).write_h5ad(f'{out}/{src}.h5ad')
        print('harmonised', src, a.shape, flush=True)
