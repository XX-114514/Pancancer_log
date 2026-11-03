# Data-consistency report for integration v2 -> malignancy v1 -> TME v1: one row per source (aggregates
# only, no per-cell or per-patient rows), plus per-cancer totals and a JSON summary.
# Columns trace each source through the pipeline: cells before the decontX filter, removed for
# contamination >= 0.5, in integration v2, malignant / uncertain / non-malignant, kept in TME v1; and the
# per-source metadata the analysis relies on (cancer, OncoTree, technology, units, patients, depth,
# genes, treatment, sampling site, fresh/frozen, QC/doublet/ambient provenance, malignancy-unit status).
# Usage: report_dataset_detail_v1.py V2_SHARD_DIR V2_INTEGRATION_SUMMARY MALIG_DIR TME_REPORT_DIR OUT_DIR
import sys, os, glob, json
import numpy as np, pandas as pd, anndata as ad

sdir, v2sum, mdir, tdir, out = sys.argv[1:6]
os.makedirs(out, exist_ok=True)
filt = json.load(open(v2sum))['filter']
EXCLUDE = {'GSE179994': 'excluded_lineage_enriched_all_immune'}  # mirrors integrate_tme_v1.py
dist = lambda s: '; '.join(f'{k}:{v}' for k, v in s.value_counts().items())
rows = []
for f in sorted(glob.glob(f'{sdir}/*.h5ad')):
    s = os.path.basename(f)[:-5]
    a = ad.read_h5ad(f, backed='r'); o = a.obs
    m = pd.read_csv(f'{mdir}/{s}.tsv.gz', sep='\t', usecols=['cell_id', 'malignancy_call', 'malignancy_rule', 'lineage'])
    assert len(m) == len(o) and set(m.cell_id) == set(o.index.astype(str)), s
    units = [json.loads(l) for l in open(f'{mdir}/{s}.units.jsonl')]
    mc = m.malignancy_call.value_counts()
    r = {'source': s, 'data_batch': o.data_batch.iloc[0], 'analysis_main': dist(o.analysis_main), 'oncotree_code': dist(o.oncotree_code),
         'technology': dist(o.technology), 'units': o.unit_id.nunique(), 'patients': o.patient_id.nunique(), 'genes_measured': a.n_vars,
         'cells_pre_decontX_filter': filt[s]['cells_in'], 'removed_contamination_ge_0.5': filt[s]['removed_high_contamination'],
         'cells_integration_v2': len(o), 'malignant': int(mc.get('malignant', 0)), 'uncertain': int(mc.get('uncertain', 0)),
         'non_malignant': int(mc.get('non_malignant', 0))}
    r['tme_v1_status'] = EXCLUDE.get(s, 'included')
    r['cells_tme_v1'] = 0 if s in EXCLUDE else r['uncertain'] + r['non_malignant']
    r['pct_malignant'] = round(100 * r['malignant'] / len(o), 1)
    r['pct_uncertain_of_tme'] = round(100 * r['uncertain'] / max(r['cells_tme_v1'], 1), 1)
    rule = m.malignancy_rule.value_counts()
    r['uncertain_methods_disagree'] = int(rule.get('methods_disagree', 0)); r['uncertain_cnv_not_evaluable'] = int(rule.get('CNV_not_evaluable', 0))
    r['units_scatomic_failed'] = sum(str(u.get('scatomic', '')).startswith('failed') for u in units)
    r['units_infercnv_not_done'] = sum(u.get('infercnv') != 'done' for u in units)
    for k, c in (('umi', 'total_counts'), ('genes', 'n_genes'), ('pct_mt', 'pct_mt'), ('contamination', 'decontX_contamination')):
        r[f'median_{k}'] = round(float(np.nanmedian(o[c].astype(float))), 3)
    for c in ('treatment_status', 'sampling_site', 'fresh_frozen', 'suspension', 'capture_scope', 'qc_flag', 'doublet_method',
              'ambient_correction', 'cancer_label_resolution', 'patient_id_rule'):
        r[c] = dist(o[c].astype(str))
    r['ambient_correction'] = 'decontX per unit (integration v2, all sources); origin store: ' + r['ambient_correction']
    rows.append(r)
    print(s, r['cells_integration_v2'], r['malignant'], r['uncertain'], flush=True)
D = pd.DataFrame(rows)
D.to_csv(f'{out}/dataset_detail_by_source.tsv', sep='\t', index=False)
num = ['cells_pre_decontX_filter', 'removed_contamination_ge_0.5', 'cells_integration_v2', 'malignant', 'uncertain', 'non_malignant', 'cells_tme_v1']
# per cancer (sources with several cancers are split by their cell counts in the v2 shards)
cc = []
for f in sorted(glob.glob(f'{sdir}/*.h5ad')):
    s = os.path.basename(f)[:-5]
    o = ad.read_h5ad(f, backed='r').obs[['analysis_main', 'oncotree_code', 'unit_id', 'patient_id']].copy()
    m = pd.read_csv(f'{mdir}/{s}.tsv.gz', sep='\t', usecols=['cell_id', 'malignancy_call']).set_index('cell_id').reindex(o.index.astype(str))
    o['call'] = m.malignancy_call.values; o['source'] = s
    cc.append(o)
C = pd.concat(cc)
P = C.groupby('analysis_main', observed=True).agg(sources=('source', 'nunique'), units=('unit_id', 'nunique'), patients=('patient_id', 'nunique'),
                                                   cells_integration_v2=('call', 'size'), malignant=('call', lambda x: (x == 'malignant').sum()),
                                                   uncertain=('call', lambda x: (x == 'uncertain').sum()), non_malignant=('call', lambda x: (x == 'non_malignant').sum()))
C['tme'] = (C['call'] != 'malignant') & ~C['source'].isin(list(EXCLUDE))
P['cells_tme_v1'] = C.groupby('analysis_main', observed=True)['tme'].sum()
P['oncotree_codes'] = C.groupby('analysis_main', observed=True)['oncotree_code'].agg(lambda x: ','.join(sorted(set(map(str, x)))))
P.sort_values('cells_integration_v2', ascending=False).to_csv(f'{out}/dataset_detail_by_cancer.tsv', sep='\t')
tsum = json.load(open(f'{tdir}/tme_summary.json')) if os.path.exists(f'{tdir}/tme_summary.json') else {}
summ = {'sources': len(D), 'cancers': int(P.shape[0]), 'units': int(D.units.sum()), 'patients_v2': int(C.patient_id.nunique()),
        **{k: int(D[k].sum()) for k in num}, 'pct_malignant': round(100 * D.malignant.sum() / D.cells_integration_v2.sum(), 2),
        'units_scatomic_failed': int(D.units_scatomic_failed.sum()), 'units_infercnv_not_done': int(D.units_infercnv_not_done.sum()),
        'tme_v1_check': {'cells_in_zarr': tsum.get('cells'), 'matches_sum': tsum.get('cells') == int(D.cells_tme_v1.sum()) if tsum else None}}
json.dump(summ, open(f'{out}/dataset_detail_summary.json', 'w'), indent=1)
print(json.dumps(summ, indent=1))
