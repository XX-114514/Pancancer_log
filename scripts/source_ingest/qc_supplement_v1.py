# Per-sample QC for ingested supplementary sources, using the project's native-branch rules
# (methods/project-v3-unified-workflow.md; clinical expansion native branches):
#   floor: >=200 detected genes and >=500 UMI;
#   mitochondrial: pct_mt <= min(25, max(15, median + 3*MAD)) per sample, on floor-passing cells;
#   doublets: Scrublet per sample (scanpy defaults, seed 0) on floor/mt-passing cells, predicted doublets removed;
#   samples with <100 cells after QC are dropped.
# Ambient RNA is not estimated (recorded as not_estimated). Writes the QC-passed h5ad and a per-sample report.
# Usage: qc_supplement_v1.py IN_DIR OUT_DIR [SOURCE ...]
import sys, os, json, glob
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc

inp, out = sys.argv[1:3]
only = set(sys.argv[3:])
os.makedirs(out, exist_ok=True)
rep_all = []
for f in sorted(glob.glob(os.path.join(inp, '*.h5ad'))):
    src = os.path.basename(f)[:-5]
    if only and src not in only:
        continue
    a = ad.read_h5ad(f)
    X = a.X.tocsr()
    mt = a.var_names.str.upper().str.startswith('MT-')
    tot = np.asarray(X.sum(1)).ravel()
    ng = np.diff(X.indptr)
    pmt = np.asarray(X[:, mt].sum(1)).ravel() / np.maximum(tot, 1) * 100
    a.obs['total_counts'] = tot; a.obs['n_genes'] = ng; a.obs['pct_mt'] = pmt
    a.obs['qc_floor_pass'] = (ng >= 200) & (tot >= 500)
    a.obs['qc_mt_pass'] = False; a.obs['doublet_score'] = np.nan; a.obs['predicted_doublet'] = False
    a.obs['mt_threshold'] = np.nan
    for s, idx in a.obs.groupby('sample_id', observed=True).indices.items():
        idx = np.asarray(idx)
        fl = idx[a.obs['qc_floor_pass'].values[idx]]
        if len(fl) == 0:
            continue
        v = pmt[fl]; med = np.median(v); mad = np.median(np.abs(v - med))
        thr = min(25.0, max(15.0, med + 3 * mad)) if mt.any() else np.inf
        a.obs.iloc[idx, a.obs.columns.get_loc('mt_threshold')] = thr
        ok = fl[pmt[fl] <= thr]
        a.obs.iloc[ok, a.obs.columns.get_loc('qc_mt_pass')] = True
        if len(ok) >= 100:
            sub = ad.AnnData(X=X[ok].astype(np.float32))
            try:
                sc.pp.scrublet(sub, random_state=0, verbose=False)
                a.obs.iloc[ok, a.obs.columns.get_loc('doublet_score')] = sub.obs['doublet_score'].values
                a.obs.iloc[ok, a.obs.columns.get_loc('predicted_doublet')] = sub.obs['predicted_doublet'].values.astype(bool)
            except Exception as e:  # too few cells / degenerate sample: keep unscored, record
                print(src, s, 'scrublet_failed', repr(e)[:80], flush=True)
    keep = a.obs['qc_floor_pass'] & a.obs['qc_mt_pass'] & ~a.obs['predicted_doublet'].astype(bool)
    per = a.obs.assign(keep=keep).groupby('sample_id', observed=True).agg(
        cells_in=('keep', 'size'), floor_pass=('qc_floor_pass', 'sum'), mt_pass=('qc_mt_pass', 'sum'),
        doublets=('predicted_doublet', 'sum'), cells_qc=('keep', 'sum'), mt_threshold=('mt_threshold', 'first'),
        median_umi=('total_counts', 'median'), patient_id=('patient_id', 'first'))
    per['sample_kept'] = per['cells_qc'] >= 100
    keep &= a.obs['sample_id'].map(per['sample_kept']).astype(bool)
    b = a[keep.values].copy()
    b.obs['ambient_correction'] = 'not_estimated'
    b.obs['doublet_method'] = 'scrublet_per_sample'
    b.write_h5ad(os.path.join(out, f'{src}.h5ad'), compression='gzip')
    per.insert(0, 'source', src)
    rep_all.append(per.reset_index())
    print(json.dumps({'source': src, 'cells_in': int(a.n_obs), 'cells_qc': int(b.n_obs), 'samples_in': int(len(per)),
                      'samples_kept': int(per['sample_kept'].sum()), 'patients_kept': int(b.obs['patient_id'].nunique())}), flush=True)
if rep_all:
    pd.concat(rep_all).to_csv(os.path.join(out, f'qc_report{"_" + "_".join(sorted(only)) if only else ""}.tsv'), sep='\t', index=False)
