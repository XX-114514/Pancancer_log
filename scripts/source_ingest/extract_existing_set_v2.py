# Extract analysis set v2 (qualified fit cells of the 40 retained sources) from the three native
# stores into one h5ad shard per source: raw integer counts, genes keyed by symbol (duplicates
# summed, as in the core store), only genes measured by the source, non-RNA features dropped, plus per-cell unit / patient /
# treatment / sampling-site / cancer-type metadata. Cell set = the cell axis of the current TME
# standard integration run (20261029_TME_standard_HVG_PCA_Harmony_v1) minus excluded sources.
# Usage: extract_existing_set_v2.py PC_ROOT REPO OUT_DIR [SOURCE ...]
import sys, os, json, re
import numpy as np, pandas as pd, scipy.sparse as sp, zarr, anndata as ad
from anndata.io import read_elem

pc, repo, out = sys.argv[1:4]
only = set(sys.argv[4:])
runs = f'{pc}/Project_v5/runs'
A = f'{pc}/analysis/data_gaps'
os.makedirs(out, exist_ok=True)
pm = [l.rstrip('\n').split('\t') for l in open(f'{pc}/PATH_MAP.tsv')][1:]
pm = sorted([(a, b) for a, b, *_ in pm if b.startswith('/')], key=lambda x: -len(x[0]))
local = lambda p: next((b + p[len(a):] for a, b in pm if p.startswith(a)), p)

lineage = json.load(open(f'{runs}/20261029_TME_standard_HVG_PCA_Harmony_v1/outputs/source_input_lineage_v1.json'))
axes = np.load(f'{runs}/20261029_TME_standard_HVG_PCA_Harmony_v1/outputs/private_cell_axes_v1.npz')
excluded = set(pd.read_csv(f'{repo}/inventories/data_gaps/source_exclusions_v1.tsv', sep='\t').source)
nonrna = pd.read_csv(f'{repo}/inventories/feature_exclusions/non_RNA_features_v1.tsv', sep='\t')
nonrna = set(nonrna.feature.str.upper())
cmap = pd.read_csv(f'{repo}/inventories/cancer_type/cancer_type_map_v1.tsv', sep='\t').set_index('source')

# unit-level metadata (server-private tables)
mat = pd.read_csv(f'{runs}/20261030_TME_CM_technical_material_design_v1/outputs/private_analysis_unit_context_and_coverage_v1.tsv', sep='\t')
qual = pd.read_csv(f'{runs}/20261027_TME_module_scope_registry_v2/outputs/TME_source_qualified_units_core_v7.tsv', sep='\t')
qual['parents'] = qual.native_parent_unit_sha256s.apply(json.loads)
k = 'qualification_unit_key_SHA256'
meta = mat.merge(qual[['qualification_unit_key_sha256', 'qualification_unit_id', 'parents']],
                 left_on=k, right_on='qualification_unit_key_sha256')
for f, cols in (('gap1_treatment_resolution_private_v1.tsv', ['gap1_rule', 'gap1_status']),
                ('gap2_patient_keys_private_v1.tsv', ['gap2_rule', 'patient_key']),
                ('gap5_sampling_site_private_v1.tsv', ['gap5_rule', 'sampling_site'])):
    meta = meta.merge(pd.read_csv(f'{A}/{f}', sep='\t')[[k] + cols], on=k, how='left')
ff = pd.read_csv(f'{repo}/inventories/data_gaps/fresh_frozen_evidence_v1.tsv', sep='\t').set_index('source').status


def symbol_var(store, g):
    v = read_elem(g['var'])
    sym = v['gene_symbol'].astype(str).values if 'gene_symbol' in v else v.index.astype(str).values
    return np.asarray(sym)


def treatment_harmonised(r):
    ph, st = str(r.treatment_phase), str(r.gap1_status)
    if r.gap1_rule == 'recorded':
        txt = ph
    elif r.gap1_rule == 'unresolved':
        return 'unknown', ''
    else:
        txt = st
    if re.search(r'naive|untreated|未治疗|治疗前|治疗/干预前|preRT|基线|newly_diagnosed', txt, re.I):
        return 'untreated_or_pretreatment', txt
    if re.search(r'prior|treated|recurrent|治疗中|进展后|onRT', txt, re.I):
        return 'on_or_after_treatment', txt
    return 'other', txt


summary = []
for code, L in enumerate(lineage):
    src = L['source']
    if src in excluded or (only and src not in only):
        continue
    sel = axes['source_code'] == code
    rows = axes['native_row'][sel]
    assert len(rows) == L['n_fit'], (src, len(rows), L['n_fit'])
    store = local(L['store'])
    g = zarr.open(store, 'r')
    m = g[L['counts_key']]
    indptr = m['indptr'][:]
    lo, hi = int(rows.min()), int(rows.max()) + 1
    a, b = int(indptr[lo]), int(indptr[hi])
    data = m['data'][a:b]; ind = m['indices'][a:b]
    blk = sp.csr_matrix((data, ind, indptr[lo:hi + 1] - a), shape=(hi - lo, m.attrs['shape'][1]))
    X = blk[rows - lo]
    if not np.allclose(X.data[:100000], np.round(X.data[:100000])):
        raise SystemExit(f'{src}: non-integer counts in {L["counts_key"]}')
    sym = symbol_var(store, g)
    keep_cols = ~pd.Series(sym).str.upper().isin(nonrna).values
    # keep only genes the source actually measured (stores zero-pad unmeasured genes)
    if 'measured_in_source' in g['varm']:
        order = re.findall(r"'([^']+)'", str(read_elem(g['uns']['source_order'])))
        keep_cols &= np.asarray(read_elem(g['varm']['measured_in_source']))[:, order.index(src)]
    else:
        per = ad.read_h5ad(f'{pc}/Project_v3/annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse/{src}_symbol_counts.h5ad', backed='r')
        keep_cols &= np.isin(sym, per.var_names.astype(str))
    # sum duplicate symbols
    u, inv = np.unique(sym[keep_cols], return_inverse=True)
    agg = sp.csr_matrix((np.ones(keep_cols.sum()), (np.flatnonzero(keep_cols), inv)), shape=(len(sym), len(u)))
    Xs = (X @ agg).tocsr()
    Xs.data = np.rint(Xs.data).astype(np.int32); Xs = Xs.astype(np.int32)
    # per-cell native identifiers and unit
    obs_names = np.asarray(read_elem(g['obs'][g['obs'].attrs['_index']]))[rows].astype(str)
    uv = np.load(f'{runs}/20261026_archive_expression_unit_views_v1/outputs/{src}_archive_unit_views_v1.npz')
    pos = pd.Series(np.arange(len(uv['native_row'])), index=uv['native_row'])
    unit_code = uv['native_expression_unit_code'][pos.loc[rows].values]
    native_unit = uv['expression_unit_dictionary'][unit_code]
    msrc = meta[meta.source == src].copy()
    if src == 'GSE166555':
        br = np.load(f'{runs}/20261027_GSE166555_exact_author_cell_bridge_v1/outputs/author_metadata_bridge_v1.npz')
        reg = json.load(open(f'{runs}/20261027_GSE166555_exact_author_cell_bridge_v1/outputs/author_sample_registry_v1.json'))
        sd = reg['sample_dictionary']
        bpos = pd.Series(np.arange(len(br['native_row'])), index=br['native_row'])
        samp = np.asarray(sd)[br['author_sample_code'][bpos.loc[rows].values]]
        gsm = {r['author_sample_id']: r['GSM'] for r in reg['records']}
        unit_id = np.array([next((q for q in msrc.qualification_unit_id if q.endswith(f'::{s}::{gsm.get(s, "")}')), '') for s in samp])
    else:
        unit_id = native_unit
    msrc = msrc.set_index('qualification_unit_id')
    missing = set(unit_id) - set(msrc.index)
    keep = np.isin(unit_id, msrc.index)
    th = msrc.apply(lambda r: pd.Series(treatment_harmonised(r), index=['treatment_status', 'treatment_detail']), axis=1)
    msrc = msrc.join(th)
    cols = {'unit_id': 'qualification_unit_id', 'patient_key': 'patient_key', 'patient_rule': 'gap2_rule',
            'treatment_status': 'treatment_status', 'treatment_detail': 'treatment_detail', 'treatment_rule': 'gap1_rule',
            'material_class': 'material_class', 'material': 'material', 'sampling_site': 'sampling_site',
            'sampling_site_rule': 'gap5_rule', 'fresh_frozen_recorded': 'fresh_frozen', 'RNA_modality': 'RNA_modality'}
    o = pd.DataFrame({'cell_id': [f'{src}|{x}' for x in obs_names], 'source': src, 'origin_store': os.path.basename(store),
                      'origin_counts_key': L['counts_key'], 'origin_native_row': rows, 'origin_obs_name': obs_names,
                      'qualification_unit_id': unit_id})
    o = o.join(msrc.drop(columns=[c for c in msrc.columns if c not in cols.values()]), on='qualification_unit_id')
    o['patient_key'] = o['patient_key'].fillna('')
    o['fresh_frozen'] = np.where(o['fresh_frozen'].astype(str).str.contains('待核实'),
                                 'fresh' if ff.get(src) == 'fresh' else 'unknown', o['fresh_frozen'].astype(str))
    for c in ['analysis_main', 'analysis_fine', 'oncotree_code', 'resolution']:
        o['cancer_' + c if c == 'resolution' else c] = cmap.loc[src, c]
    o['ambient_correction'] = 'available_in_origin_store' if L.get('has_adopted_ambient_delta') else 'none'
    o['data_batch'] = 'existing_analysis_set_v2'
    o = o[keep].reset_index(drop=True)
    Xs = Xs[keep]
    o.index = o['cell_id'].values
    adata = ad.AnnData(X=Xs, obs=o, var=pd.DataFrame(index=u))
    adata.write_h5ad(os.path.join(out, f'{src}.h5ad'), compression='gzip')
    summary.append({'source': src, 'cells': int(adata.n_obs), 'genes': int(adata.n_vars), 'units': int(o.qualification_unit_id.nunique()),
                    'unmatched_unit_cells': int((~keep).sum()), 'missing_units': sorted(missing)[:3]})
    print(json.dumps(summary[-1]), flush=True)
json.dump(summary, open(os.path.join(out, f'extract_summary{"_" + "_".join(sorted(only)) if only else ""}.json'), 'w'), indent=1)
