"""Independently reconstruct native ligand counts and means using CSC columns."""
import hashlib
import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from diagnose_frozen_MCP_two_end_functions_v1 import hnscc_inputs


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    project = Path(os.environ['PANCANCER_PROJECT_ROOT'])
    out = root/'HNSCC_all_sender_ligand_measurements_20261016_v1'
    audit_path = out/'private_independent_audit.json'
    assert not audit_path.exists()
    terminal = json.loads((out/'private_terminal.json').read_text())
    expression = pd.read_csv(out/'private_all_patient_group_gene_expression.tsv', sep='\t')
    meta = pd.read_csv(out/'private_all_patient_group_capture_quality.tsv', sep='\t')
    coverage = pd.read_csv(out/'all400_private_sender_ligand_coverage.tsv', sep='\t')
    receptors = pd.read_csv(out/'private_all_patient_receptor_complex_RNA.tsv', sep='\t')
    for name, key in [('private_all_patient_group_gene_expression.tsv', 'expression_sha256'),
                      ('private_all_patient_group_capture_quality.tsv', 'quality_sha256'),
                      ('all400_private_sender_ligand_coverage.tsv', 'coverage_sha256')]:
        assert hashlib.sha256((out/name).read_bytes()).hexdigest() == terminal[key]
    assert not expression.duplicated(['sender_group', 'patient_id', 'gene']).any()
    assert not meta.duplicated(['sender_group', 'patient_id']).any()
    g, source_rows, cells, var, eligible, measured = hnscc_inputs(root, project)
    genes = sorted(expression.gene.unique())
    matches = defaultdict(set)
    wanted = set(genes)
    for i, (name, row) in enumerate(var.iterrows()):
        if not measured[i] or str(name).startswith('AMBIGUOUS::'):
            continue
        known = {str(name), str(row.get('gene_symbol', ''))}
        if row.get('feature_identity_status') == 'stable_Ensembl_ID':
            known.update(str(row.get('source_gene_aliases', '')).split('|'))
        for gene in known & wanted:
            matches[gene].add(i)
    columns = {gene: next(iter(matches[gene])) if len(matches[gene]) == 1 else None for gene in genes}
    keys = list(zip(meta.sender_group, meta.patient_id))
    key_index = {key: i for i, key in enumerate(keys)}
    patients = sorted(meta.patient_id.unique())
    base = cells.patient_id.isin(patients) & cells.timepoint.eq('pre') & cells.library_id.isin(eligible.library_id)
    positions = np.flatnonzero(base.to_numpy())
    groups = cells.iloc[positions].major.astype(str).to_numpy().astype(object)
    is_receiver = cells.iloc[positions].major.eq('Myeloid').to_numpy() & cells.iloc[positions].subtype.eq('Macro_C1QC').to_numpy() & cells.iloc[positions].evidence.eq('multigene_subtype_candidate').to_numpy()
    groups[groups == 'Myeloid'] = 'Myeloid_excluding_C1QC'
    groups[is_receiver] = 'Receiver_C1QC'
    ids = cells.iloc[positions].patient_id.to_numpy()
    codes = np.array([key_index.get((group, patient), -1) for group, patient in zip(groups, ids)])
    counts = np.bincount(codes[codes >= 0], minlength=len(keys))
    assert np.array_equal(counts, meta.native_cells)
    totalTME = cells.iloc[positions].groupby('patient_id').size()
    totalMyeloid = cells.iloc[positions][cells.iloc[positions].major.eq('Myeloid')].groupby('patient_id').size()
    for row in meta.itertuples():
        assert row.retained_TME_cells == totalTME[row.patient_id]
        assert row.Myeloid_cells == totalMyeloid.get(row.patient_id, 0)
        assert np.isclose(row.retained_TME_fraction, row.native_cells/row.retained_TME_cells)
        if row.sender_group == 'Receiver_C1QC':
            assert np.isclose(row.receiver_fraction_within_Myeloid, row.native_cells/row.Myeloid_cells)
    physical = sorted({c for c in columns.values() if c is not None})
    column_index = {c: i for i, c in enumerate(physical)}
    hit = np.zeros((len(keys), len(physical)), np.int64)
    amount = np.zeros((len(keys), len(physical)))
    qRNA = np.zeros(len(keys)); qgenes = np.zeros(len(keys))
    receptor_positive = {c: set() for c in physical}
    matrix = ad.io.sparse_dataset(g['layers/counts'])
    for start in range(0, len(positions), 1500):
        pos = positions[start:start+1500]
        raw = matrix[source_rows[pos], :].tocsr()
        raw.sum_duplicates(); raw.eliminate_zeros()
        assert (raw.data > 0).all() and np.isfinite(raw.data).all()
        totals = np.asarray(raw.sum(axis=1)).ravel()
        number = np.asarray((raw > 0).sum(axis=1)).ravel()
        z = codes[start:start+len(pos)]
        valid = z >= 0
        qRNA += np.bincount(z[valid], weights=np.log1p(totals[valid]), minlength=len(keys))
        qgenes += np.bincount(z[valid], weights=np.log1p(number[valid]), minlength=len(keys))
        csc = raw[:, physical].tocsc()
        for j, c in enumerate(physical):
            lo, hi = csc.indptr[j:j+2]
            rows = csc.indices[lo:hi]; vals = csc.data[lo:hi]
            use = z[rows] >= 0
            hit[:, j] += np.bincount(z[rows[use]], minlength=len(keys))
            amount[:, j] += np.bincount(z[rows[use]], weights=np.log1p(10000*vals[use]/totals[rows[use]]), minlength=len(keys))
            rec = is_receiver[start+rows]
            receptor_positive[c].update((start+rows[rec]).tolist())
    use = counts > 0
    assert np.allclose(qRNA[use]/counts[use], meta.loc[use, 'mean_logRNA'], atol=1e-10, rtol=1e-10)
    assert np.allclose(qgenes[use]/counts[use], meta.loc[use, 'mean_loggenes'], atol=1e-10, rtol=1e-10)
    zero_values = 0; missing_values = 0
    for row in expression.itertuples():
        k = key_index[row.sender_group, row.patient_id]; c = columns[row.gene]
        assert bool(row.unique_measured_column) == (c is not None)
        if c is None or counts[k] == 0:
            assert pd.isna(row.mean_log1p_CPTT) and pd.isna(row.positive_cells)
            missing_values += 1
        else:
            j = column_index[c]
            assert row.native_column == c and row.positive_cells == hit[k, j]
            assert np.isclose(row.mean_log1p_CPTT, amount[k, j]/counts[k], atol=1e-11, rtol=1e-10)
            assert np.isclose(row.positive_fraction, hit[k, j]/counts[k], atol=1e-12)
            zero_values += int(hit[k, j] == 0)
    joint_sets = {}
    for receptor in receptors.receptor.unique():
        chain = [columns[gene] for gene in receptor.split('_')]
        joint_sets[receptor] = set.intersection(*(receptor_positive[c] for c in chain)) if all(c is not None for c in chain) and len(set(chain)) == len(chain) else None
    for row in receptors.itertuples():
        positive = joint_sets[row.receptor]
        if positive is None:
            assert pd.isna(row.same_cell_positive_cells)
        else:
            n = sum(ids[p] == row.patient_id for p in positive)
            assert row.same_cell_positive_cells == n
            assert np.isclose(row.same_cell_positive_fraction, n/row.native_receiver_cells, atol=1e-12)
    previous = pd.read_csv(root/'HNSCC_receiver_prior_ligand_ranking_20261016_v1/private_all_sender_patient_capture.tsv', sep='\t').set_index(['sender_group', 'patient_id'])
    for row in meta[meta.sender_group.ne('Receiver_C1QC')].itertuples():
        assert row.native_cells == previous.loc[(row.sender_group, row.patient_id), 'n_cells']
    for row in coverage.itertuples():
        f = expression[expression.sender_group.eq(row.sender_group) & expression.gene.eq(row.ligand)]
        u = f[f.native_cells.ge(20) & f.status.eq('measured')]
        assert len(u) == row.measured_matched_patients
        assert int(u.positive_cells.gt(0).sum()) == row.patients_with_any_positive_ligand_cells
        assert int((u.positive_cells*100 >= u.native_cells).sum()) == row.patients_with_atleast1pct_positive_ligand_cells
        assert int((u.positive_cells*10 >= u.native_cells).sum()) == row.patients_with_atleast10pct_positive_ligand_cells
        if len(u):
            assert np.isclose(u.positive_cells.sum()/u.native_cells.sum(), row.pooled_sender_positive_fraction)
        else:
            assert pd.isna(row.pooled_sender_positive_fraction)
    assert len(coverage) == 400
    report = dict(status='PASS_ALL400_LIGAND_SLOTS_CSC_COUNTS_FULLRNA_MEANS_COMPLEX_PATIENT_COUNTS_AND_MISSINGNESS',
        finished_UTC=datetime.now(timezone.utc).isoformat(), expression_rows=len(expression), receptor_patient_rows=len(receptors),
        measured_true_zero_rows=zero_values, unavailable_or_no_cell_rows=missing_values,
        coverage_sha256=terminal['coverage_sha256'], expression_sha256=terminal['expression_sha256'],
        quality_sha256=terminal['quality_sha256'], no_effect_tests_or_candidate_reselection=True)
    audit_path.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
