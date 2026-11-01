# Build the combined analysis-ready zarr (pancancer_analysis_ready_v1.zarr) from
#   - existing analysis set v2 shards (extract_existing_set_v2.py; 40 sources), and
#   - supplementary sources after QC (ingest_supplement_v1.py + qc_supplement_v1.py; 13 sources).
# Steps: (1) rewrite each shard's obs to one harmonised schema and recompute QC metrics from X;
# (2) concatenate all shards on disk (outer gene join, X = raw integer counts, unmeasured genes are
# structural zeros); (3) add the gene x source measured mask, gene coverage columns and provenance;
# (4) write the composition tables used by the data description.
# Usage: build_analysis_zarr_v1.py EXISTING_DIR SUPPLEMENT_QC_DIR REPO OUT_ZARR REPORT_DIR
import sys, os, json, glob, hashlib, re, subprocess, datetime
import numpy as np, pandas as pd, h5py, zarr, anndata as ad
from anndata.io import read_elem, write_elem
from anndata.experimental import concat_on_disk

exist_dir, supp_dir, repo, out_zarr, rep_dir = sys.argv[1:6]
os.makedirs(rep_dir, exist_ok=True)
work = out_zarr + '.shards'
os.makedirs(work, exist_ok=True)
cmap = pd.concat([pd.read_csv(f'{repo}/inventories/cancer_type/cancer_type_map_v1.tsv', sep='\t'),
                  pd.read_csv(f'{repo}/inventories/cancer_type/cancer_type_map_supplement_v1.tsv', sep='\t')]).set_index('source')
ONCO_MAIN = {'Hepatocellular carcinoma': 'HCC', 'Intrahepatic cholangiocarcinoma': 'IHCH'}
COLS = ['cell_id', 'source', 'data_batch', 'unit_id', 'sample_id', 'patient_id', 'patient_id_rule',
        'analysis_main', 'analysis_fine', 'oncotree_code', 'cancer_label_resolution',
        'sampling_site', 'material_detail', 'treatment_status', 'treatment_detail', 'treatment_rule',
        'fresh_frozen', 'fresh_frozen_detail', 'technology', 'suspension', 'capture_scope',
        'origin_store', 'origin_counts_key', 'origin_native_row', 'origin_cell_name',
        'ambient_correction', 'doublet_method', 'author_cell_type', 'qc_flag',
        'total_counts', 'n_genes', 'pct_mt']
SITE = {'primary_tumor': 'primary_organ_site', 'metastatic_tumor': 'metastatic_site', 'involved_lymph_node': 'metastatic_site',
        'primary_organ_site': 'primary_organ_site', 'metastatic_site': 'metastatic_site', 'tumour_thrombus': 'tumour_thrombus'}


def fresh(v):
    v = str(v)
    if re.search(r'冻存|cryopreserv|frozen|thaw', v, re.I) and re.search(r'新鲜|fresh', v, re.I):
        return 'fresh_or_cryopreserved'
    if re.search(r'冻存|cryopreserv|frozen|thaw', v, re.I):
        return 'cryopreserved_or_frozen'
    if re.search(r'新鲜|fresh', v, re.I):
        return 'fresh'
    return 'unknown'


def harmonise(o, src, batch):
    c = cmap.loc[src]
    h = pd.DataFrame(index=o.index)
    h['cell_id'] = o['cell_id'].astype(str)
    h['source'] = src
    h['data_batch'] = batch
    if batch == 'existing_analysis_set_v2':
        h['unit_id'] = o['qualification_unit_id'].astype(str)
        h['sample_id'] = o['qualification_unit_id'].astype(str)
        pk = o['patient_key'].astype(str)
        # patient keys are source-scoped; unknown patient -> the unit stands in for the patient
        h['patient_id'] = np.where(pk == '', src + ':unit:' + h['unit_id'], np.where(pk.str.startswith(src), pk, src + ':' + pk))
        h['patient_id_rule'] = o['patient_rule'].astype(str) if 'patient_rule' in o else o.get('gap2_rule', 'recorded').astype(str)
        h['analysis_main'] = c['analysis_main']
        h['analysis_fine'] = c['analysis_fine']
        h['oncotree_code'] = c['oncotree_code']
        h['sampling_site'] = o['sampling_site'].astype(str).where(o['sampling_site'].astype(str) != '', o['material_class'].astype(str)).map(lambda s: SITE.get(s, 'unspecified'))
        h['material_detail'] = o['material'].astype(str)
        h['treatment_rule'] = o['treatment_rule'].astype(str) if 'treatment_rule' in o else o['gap1_rule'].astype(str)
        h['fresh_frozen_detail'] = o['fresh_frozen'].astype(str)
        h['technology'] = 'scRNA (see source registry)'
        h['suspension'] = 'cell'
        h['capture_scope'] = 'source_gate_v7_qualified'
        h['origin_store'] = o['origin_store'].astype(str)
        h['origin_counts_key'] = o['origin_counts_key'].astype(str)
        h['origin_native_row'] = o['origin_native_row'].astype(np.int64)
        h['origin_cell_name'] = o['origin_obs_name'].astype(str)
        h['doublet_method'] = 'inherited_from_origin_store'
        h['author_cell_type'] = ''
        h['qc_flag'] = 'inherited_fit_cell'
    else:
        h['unit_id'] = src + '|' + o['sample_id'].astype(str)
        h['sample_id'] = o['sample_id'].astype(str)
        h['patient_id'] = o['patient_id'].astype(str)
        h['patient_id_rule'] = 'GEO_sample_record'
        main = o['analysis_main'].astype(str) if 'analysis_main' in o else pd.Series(c['analysis_main'], index=o.index)
        h['analysis_main'] = main
        h['analysis_fine'] = main.map(lambda m: {'Hepatocellular carcinoma': 'HCC', 'Intrahepatic cholangiocarcinoma': 'ICC'}.get(m, c['analysis_fine']))
        h['oncotree_code'] = main.map(lambda m: ONCO_MAIN.get(m, c['oncotree_code'])) if src in ('GSE125449', 'GSE189903') else c['oncotree_code']
        if src == 'GSE277524':
            h.loc[o['sample_id'].astype(str) == 'P03_R1', 'oncotree_code'] = 'UTUC'
        h['sampling_site'] = o['sampling_site'].astype(str).map(lambda s: SITE.get(s, 'unspecified'))
        h['material_detail'] = o['material'].astype(str)
        h['treatment_rule'] = 'GEO_sample_record_or_paper'
        h['fresh_frozen_detail'] = o['fresh_frozen'].astype(str)
        h['technology'] = o['technology'].astype(str)
        h['suspension'] = o['suspension'].astype(str)
        h['capture_scope'] = o['capture_scope'].astype(str)
        h['origin_store'] = f'GEO:{src}'
        h['origin_counts_key'] = 'supplementary_file'
        h['origin_native_row'] = -1
        h['origin_cell_name'] = o['origin_barcode'].astype(str)
        h['doublet_method'] = o['doublet_method'].astype(str)
        h['author_cell_type'] = o['author_cell_type'].astype(str) if 'author_cell_type' in o else ''
        h['qc_flag'] = 'low_depth_platform' if src == 'GSE163203' else 'qc_v1_pass'
    h['cancer_label_resolution'] = c['resolution']
    h['treatment_status'] = o['treatment_status'].astype(str)
    h['treatment_detail'] = o['treatment_detail'].astype(str)
    h['fresh_frozen'] = h['fresh_frozen_detail'].map(fresh)
    h['ambient_correction'] = o['ambient_correction'].astype(str)
    return h


# (1) harmonise shards (copies in work dir; original shards untouched)
shards, gene_sets = {}, {}
for batch, d in (('existing_analysis_set_v2', exist_dir), ('supplement_v1', supp_dir)):
    for f in sorted(glob.glob(os.path.join(d, '*.h5ad'))):
        src = os.path.basename(f)[:-5]
        a = ad.read_h5ad(f)
        X = a.X.tocsr()
        o = harmonise(a.obs, src, batch)
        mt = a.var_names.str.upper().str.startswith('MT-')
        o['total_counts'] = np.asarray(X.sum(1)).ravel().astype(np.int64)
        o['n_genes'] = np.diff(X.indptr).astype(np.int32)
        o['pct_mt'] = (np.asarray(X[:, mt].sum(1)).ravel() / np.maximum(o['total_counts'].values, 1) * 100).astype(np.float32)
        o = o[COLS]
        o.index = o['cell_id'].values
        b = ad.AnnData(X=X, obs=o, var=pd.DataFrame(index=a.var_names.astype(str)))
        p = os.path.join(work, f'{src}.h5ad')
        b.write_h5ad(p)
        shards[src] = p
        gene_sets[src] = set(b.var_names)
        print('harmonised', src, b.shape, flush=True)

# (2) on-disk concatenation
order = list(shards)
concat_on_disk([shards[s] for s in order], out_zarr, join='outer', max_loaded_elems=200_000_000)

# (3) var coverage, measured mask, obs categoricals, provenance
g = zarr.open(out_zarr, 'r+')
var = read_elem(g['var'])
genes = var.index.astype(str)
mask = np.zeros((len(genes), len(order)), dtype=bool)
gi = pd.Index(genes)
for j, s in enumerate(order):
    mask[gi.get_indexer(sorted(gene_sets[s])), j] = True
var['n_sources_measured'] = mask.sum(1).astype(np.int32)
var['measured_in_all_sources'] = var['n_sources_measured'] == len(order)
existing = [i for i, s in enumerate(order) if os.path.exists(os.path.join(exist_dir, s + '.h5ad'))]
var['measured_in_all_existing_sources'] = mask[:, existing].all(1)
del g['var']; write_elem(g, 'var', var)
if 'varm' not in g:
    g.create_group('varm')
write_elem(g['varm'], 'measured_in_source', mask)
obs = read_elem(g['obs'])
for col in obs.columns:
    if obs[col].dtype == object:
        obs[col] = obs[col].astype('category')
del g['obs']; write_elem(g, 'obs', obs)
commit = subprocess.run(['git', '-C', repo, 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
uns = {
    'schema': 'pancancer_analysis_ready_v1',
    'created_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'code_commit': commit,
    'X': 'raw integer UMI/read counts per source (int32); genes keyed by symbol, duplicate symbols summed; non-RNA features removed',
    'unmeasured_genes': 'structural zeros: use varm/measured_in_source; never treat as biological zero',
    'source_order': order,
    'existing_cells': 'qualified fit cells of the TME standard integration (20261029_TME_standard_HVG_PCA_Harmony_v1) minus excluded sources (inventories/data_gaps/source_exclusions_v1.tsv); counts from the origin store key in obs origin_counts_key (no ambient subtraction applied here)',
    'supplement_cells': 'GEO supplementary counts, tumour samples only, QC v1 (>=200 genes, >=500 UMI, per-sample MT cap min(25,max(15,median+3MAD)), Scrublet per sample, samples >=100 cells); ambient not estimated',
    'not_done': 'no normalisation, batch correction, embedding, cell-type annotation or cross-source de-duplication of patients beyond the fingerprint screen',
}
write_elem(g, 'uns', {k: json.dumps(v) if not isinstance(v, str) else v for k, v in uns.items()})
zarr.consolidate_metadata(out_zarr)

# (4) composition tables
obs = obs.reset_index(drop=True)
def comp(by):
    return obs.groupby(by, observed=True).agg(cells=('cell_id', 'size'), units=('unit_id', 'nunique'), patients=('patient_id', 'nunique'),
                                              sources=('source', 'nunique')).reset_index()
comp(['data_batch']).to_csv(f'{rep_dir}/composition_by_batch.tsv', sep='\t', index=False)
comp(['analysis_main']).sort_values('cells', ascending=False).to_csv(f'{rep_dir}/composition_by_cancer.tsv', sep='\t', index=False)
comp(['source', 'data_batch', 'analysis_main']).to_csv(f'{rep_dir}/composition_by_source.tsv', sep='\t', index=False)
for col in ['sampling_site', 'treatment_status', 'fresh_frozen', 'technology', 'capture_scope', 'qc_flag']:
    comp([col]).to_csv(f'{rep_dir}/composition_by_{col}.tsv', sep='\t', index=False)
summ = {'cells': int(len(obs)), 'genes': int(len(genes)), 'sources': len(order), 'units': int(obs.unit_id.nunique()),
        'patients': int(obs.patient_id.nunique()), 'genes_measured_in_all_sources': int(var['measured_in_all_sources'].sum()),
        'cancers': int(obs.analysis_main.nunique())}
json.dump(summ, open(f'{rep_dir}/summary.json', 'w'), indent=1)
print(json.dumps(summ))
