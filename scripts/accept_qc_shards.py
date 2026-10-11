#!/usr/bin/env python3
"""Independent full-cell QC acceptance audit; never edits counts or grants source AC.

Run in the data environment (numpy/pandas/scipy/h5py/zarr/anndata). Outputs
aggregate evidence only. Count matrices are scanned in bounded CSR row blocks.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import zarr
from anndata._io.specs import read_elem


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024**2), b''):
            value.update(block)
    return value.hexdigest()


def strings(group, field):
    return np.asarray(read_elem(group[field])).astype(str)


def index(group):
    return strings(group, group.attrs['_index'])


def block(group, pointers, lo, hi, n_genes):
    start, end = int(pointers[lo]), int(pointers[hi])
    return sparse.csr_matrix((np.asarray(group['data'][start:end]),
                              np.asarray(group['indices'][start:end]),
                              pointers[lo:hi+1] - start), shape=(hi-lo, n_genes))


def one(run, proof_path, global_ids, db):
    proof = json.loads(proof_path.read_text()); name = proof['gse_id']
    print(json.dumps({'dataset': name, 'phase': 'ledger_and_policy'}), flush=True)
    paths = {k: run / 'manifests' / (name + suffix) for k, suffix in {
        'ledger': '_cell_qc_and_filter_reason.tsv.gz',
        'thresholds': '_sample_qc_thresholds.tsv',
        'registry': '_sample_registry.tsv'}.items()}
    ledger = pd.read_csv(paths['ledger'], sep='\t', keep_default_na=False)
    thresholds = pd.read_csv(paths['thresholds'], sep='\t', keep_default_na=False)
    registry = pd.read_csv(paths['registry'], sep='\t', keep_default_na=False)
    keep = ledger['keep'].astype(str).str.lower().isin(['true', '1']).to_numpy()
    selected = np.flatnonzero(keep)
    checks = {'raw_ledger_rows': len(ledger) == proof['n_raw'],
              'retained_ledger_rows': len(selected) == proof['n_qc_nonduplicate'],
              'threshold_library_unique': thresholds.library_id.is_unique,
              'registry_library_unique': registry.library_id.is_unique,
              'registry_and_threshold_libraries_equal': set(registry.library_id) == set(thresholds.library_id),
              'library_policy_recomputed': True, 'library_raw_sizes': True,
              'all_library_doublet_computed': bool(thresholds.doublet_status.eq('computed').all()),
              'all_library_mt_measured': bool(thresholds.mt_metric_status.eq('measured').all())}
    limits = dict(zip(thresholds.library_id, thresholds.mt_threshold.astype(float)))
    for lib, frame in ledger.groupby('library_id', sort=False):
        t = thresholds.loc[thresholds.library_id.eq(lib)]
        checks['library_raw_sizes'] &= len(t) == 1 and int(t.iloc[0].n_raw) == len(frame)
        floor = frame.n_genes.ge(200) & frame.total_counts.ge(500)
        mt = frame.loc[floor, 'pct_counts_mt'].astype(float).to_numpy()
        expected = min(25, max(15, float(np.median(mt) + 3*np.median(np.abs(mt-np.median(mt))))))
        checks['library_policy_recomputed'] &= bool(np.isclose(limits[lib], expected, atol=1e-10))
    kept = ledger.loc[keep]
    checks['retained_pass_floor_and_mt'] = bool((kept.n_genes.ge(200) & kept.total_counts.ge(500) &
        kept.pct_counts_mt.le(kept.library_id.map(limits) + 1e-10)).all())
    checks['retained_singlets_finite_scores'] = bool(kept.doublet_status.eq('predicted_singlet').all() and
        np.isfinite(pd.to_numeric(kept.doublet_score, errors='coerce')).all())
    output = Path(proof['output']); output.relative_to(run)
    z = zarr.open_group(str(output), mode='r'); obs = z['obs']; ids = index(obs)
    originals = strings(obs, 'original_cell_id'); libraries = strings(obs, 'original_library_id')
    checks['stable_ID_unique_within_and_across_shards'] = len(set(ids)) == len(ids) and not bool(set(ids) & global_ids)
    global_ids.update(ids)
    checks['stable_ID_exact_source_position'] = all(x == name+'::'+str(i)+'::'+str(c)
        for x, i, c in zip(ids, selected, kept.original_cell_id))
    checks['original_IDs_match_ledger'] = np.array_equal(originals, kept.original_cell_id.astype(str))
    checks['original_libraries_match_ledger'] = np.array_equal(libraries, kept.library_id.astype(str))
    lookup = registry.set_index('library_id')
    for field in ['sample_uid', 'patient_uid', 'clinical_tier', 'patient_identity_status', 'response',
                  'treatment', 'timepoint', 'tissue', 'sample_type']:
        checks['registry_obs_'+field] = np.array_equal(strings(obs, field), libraries_to_field(libraries, lookup, field))
    genes = index(z['var']); n_genes = len(genes)
    checks['gene_symbols_unique_nonempty'] = len(set(genes)) == n_genes and all(g and g != 'nan' for g in genes)
    source = Path(str(read_elem(z['uns']['raw_counts_source'])))
    source.relative_to(run)
    checks['source_sha_matches_registered_proof'] = digest(source) == proof['source_sha256']
    pre_path = run / 'manifests' / ('preflight_'+name+'.json')
    pre = json.loads(pre_path.read_text())
    official = run / pre['official_metadata_path']
    checks['official_metadata_sha_matches_proof'] = digest(official) == proof['official_source_metadata_sha256']
    checks['db_registered_cells_equal_proof'] = db.execute('SELECT COUNT(*) FROM new_cells WHERE gse=?', (name,)).fetchone()[0] == len(selected)
    numeric_checks = ['all_counts_finite_nonnegative_integer', 'X_equals_counts_layer',
                      'retained_counts_equal_original_source', 'recomputed_QC_metrics_match_ledger',
                      'stored_obs_QC_metrics_match_ledger']
    checks.update({key: True for key in numeric_checks})
    stored_metrics = {field: np.asarray(read_elem(obs[field])).astype(float)
                      for field in ['total_counts', 'n_genes', 'pct_counts_mt', 'doublet_score']}
    x_ptr = np.asarray(z['X']['indptr']); c_ptr = np.asarray(z['layers']['counts']['indptr'])
    checks['CSR_pointer_integrity'] = bool(x_ptr[0] == 0 and x_ptr[-1] == len(z['X']['data']) and
        len(x_ptr) == len(selected)+1 and (np.diff(x_ptr) >= 0).all() and np.array_equal(x_ptr, c_ptr))
    total_sum = total_nnz = 0
    mt_mask = np.char.startswith(np.char.upper(genes), 'MT-')
    with h5py.File(source, 'r') as h:
        hobs = h['obs']; source_ids = index(hobs); source_genes = index(h['var'])
        checks['ledger_original_ID_order_equals_source'] = np.array_equal(source_ids, ledger.original_cell_id.astype(str))
        checks['source_gene_order_equals_QC'] = np.array_equal(source_genes, genes)
        sg = h['layers/counts'] if 'layers/counts' in h else h['X']
        if sg.attrs['encoding-type'] != 'csr_matrix':
            raise ValueError('Source CSR contract requires separate review')
        src_ptr = np.asarray(sg['indptr'])
        for start in range(0, len(selected), 5000):
            end = min(start+5000, len(selected))
            x = block(z['X'], x_ptr, start, end, n_genes)
            c = block(z['layers']['counts'], c_ptr, start, end, n_genes)
            checks['all_counts_finite_nonnegative_integer'] &= bool(np.isfinite(x.data).all() and
                (x.data >= 0).all() and (np.abs(x.data-np.rint(x.data)) <= 1e-5).all())
            checks['X_equals_counts_layer'] &= bool(np.array_equal(x.data, c.data) and np.array_equal(x.indices, c.indices))
            rows = selected[start:end]; lo, hi = int(rows[0]), int(rows[-1])+1
            src = block(sg, src_ptr, lo, hi, n_genes)[rows-lo]
            checks['retained_counts_equal_original_source'] &= (src-x).nnz == 0
            totals = np.asarray(x.sum(axis=1)).ravel(); ng = x.getnnz(axis=1)
            pct = np.asarray(x[:, mt_mask].sum(axis=1)).ravel()*100/totals
            frame = kept.iloc[start:end]
            checks['recomputed_QC_metrics_match_ledger'] &= bool(np.array_equal(totals, frame.total_counts) and
                np.array_equal(ng, frame.n_genes) and np.allclose(pct, frame.pct_counts_mt, rtol=0, atol=1e-9))
            for field in ['total_counts', 'n_genes', 'pct_counts_mt', 'doublet_score']:
                checks['stored_obs_QC_metrics_match_ledger'] &= bool(np.allclose(
                    stored_metrics[field][start:end], frame[field].astype(float), rtol=0, atol=1e-9))
            total_sum += int(totals.sum()); total_nnz += int(x.nnz)
            if start % 50000 == 0:
                print(json.dumps({'dataset': name, 'phase': 'full_counts_comparison', 'rows_verified': end,
                                  'rows_total': len(selected)}), flush=True)
    removed = Counter(ledger.loc[~keep, 'filter_reason'])
    result = {'dataset': name, 'raw_cells': len(ledger), 'QC_cells': len(selected),
              'QC_retention_percent': round(100*len(selected)/len(ledger), 4), 'n_genes': n_genes,
              'n_libraries': len(registry), 'original_modality': pre['original_modality'],
              'filtered_reasons': dict(removed), 'raw_count_sum_QC': total_sum, 'QC_nonzeros': total_nnz,
              'library_clinical_match_status': dict(Counter(registry.clinical_match_status)),
              'library_patient_identity_status': dict(Counter(registry.patient_identity_status)),
              'tier_cell_counts': dict(Counter(strings(obs, 'clinical_tier'))),
              'technical_checks': {k: bool(v) for k, v in checks.items()},
              'source_acceptance_status': 'pending_original_study_alias_and_library_overlap_review',
              'final_accepted_new_cells': 0,
              'evidence_sha256': {p.name: digest(p) for p in [proof_path, pre_path, official]+list(paths.values())},
              'technical_audit_status': 'PASS' if all(checks.values()) else 'FAIL',
              'limits': ['Scrublet calls were verified, not rerun; sensitivity is not established.',
                         'Clinical mapping missingness remains explicit; no outcome inferred.',
                         'Exact full-count/barcode nonmatch alone does not close source overlap.']}
    return result


def libraries_to_field(libraries, lookup, field):
    return lookup[field].reindex(libraries).astype(str).to_numpy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True); parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); run = args.run.resolve()
    if args.output.exists(): raise FileExistsError('Preserve earlier acceptance evidence')
    report = {'schema': 'independent_QC_acceptance_v1', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'audit_scope': 'all_pending_QC_shards_full_cell_and_count_comparison', 'shards': [],
              'source_final_acceptance_granted': False}
    db = sqlite3.connect('file:'+str(run/'manifests/count_fingerprints.sqlite')+'?mode=ro', uri=True)
    global_ids = set()
    for path in sorted((run/'manifests').glob('accepted_*.json')):
        report['shards'].append(one(run, path, global_ids, db))
    if not report['shards']: raise ValueError('No pending QC shards')
    report['finished_utc'] = datetime.now(timezone.utc).isoformat()
    report['technical_QC_pass_cells'] = sum(s['QC_cells'] for s in report['shards'] if s['technical_audit_status']=='PASS')
    report['all_technical_checks_pass'] = all(s['technical_audit_status']=='PASS' for s in report['shards'])
    with args.output.open('x') as handle: json.dump(report, handle, indent=2, ensure_ascii=False); handle.write('\n')
    print(json.dumps({'phase':'audit_finished','all_technical_checks_pass':report['all_technical_checks_pass'],
                      'technical_QC_pass_cells':report['technical_QC_pass_cells']}), flush=True)
    raise SystemExit(0 if report['all_technical_checks_pass'] else 1)


if __name__ == '__main__': main()
