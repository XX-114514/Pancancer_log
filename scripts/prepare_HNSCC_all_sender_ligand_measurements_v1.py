"""Measure the frozen ligand universe in every prospectively eligible sender.

All native baseline cells are used. No association tests, clinical selection,
replacement ligands, inferred zero for unavailable genes, or new downloads.
"""
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


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    project = Path(os.environ['PANCANCER_PROJECT_ROOT'])
    prior = root/'HNSCC_receiver_prior_ligand_ranking_20261016_v1'
    qualification = root/'HNSCC_prior_receptor_qualification_20261016_v1'
    prior_audit = json.loads((prior/'private_independent_audit.json').read_text())
    receptor_audit = json.loads((qualification/'private_independent_audit.json').read_text())
    ranking_path = prior/'private_all_ligand_prior_rankings.tsv'
    qualified_path = qualification/'private_all_frozen_ligand_qualification.tsv'
    assert sha(ranking_path) == prior_audit['prior_ranking_sha256']
    assert sha(qualified_path) == receptor_audit['qualified_table_sha256']
    ranks = pd.read_csv(ranking_path, sep='\t')
    selected = ranks[ranks.prospective_top20].sort_values(['programme_sign', 'ligand'])
    assert len(selected) == 40 and selected.ligand.is_unique
    senders = pd.read_csv(prior/'private_sender_capture_summary.tsv', sep='\t')
    senders = senders[senders.prospective_sender_eligible].sort_values('sender_group')
    assert len(senders) == 10
    receptor_edges = pd.read_csv(qualification/'private_all_shortlist_receptor_edges.tsv', sep='\t')
    qualified = pd.read_csv(qualified_path, sep='\t').set_index('ligand')
    genes = sorted(set(selected.ligand) | {g for r in receptor_edges.receptor for g in r.split('_')})
    out = root/'HNSCC_all_sender_ligand_measurements_20261016_v1'
    assert not out.exists()
    out.mkdir()
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(Path(__file__)),
        ranking_sha256=sha(ranking_path), receptor_qualification_sha256=sha(qualified_path),
        universe='Frozen forty ligand/sign entries x all ten originally capture-eligible major sender groups; no lower-rank replacements.',
        receiver='All native pre Macro_C1QC multigene candidate cells in the original twenty-six receiver patients.',
        sender='All native pre cells in eligible libraries and matched receiver patients; remove all C1QC receiver cells from Myeloid.',
        expression='Mean per-cell log1p(10000 * gene count / full native cell RNA count), plus positive cell counts and fractions.',
        measurement='Only unique source-measured physical columns; unavailable/ambiguous genes remain missing. Measured true zeros retained.',
        coverage='All original patients, including absent sender groups; absent groups have zero cells and missing expression. Report >=20-cell coverage and any/1pct/10pct detection without selecting a favorable threshold.',
        receptor='All constituent genes and same-cell joint counts per patient, all literal reference edges; constituent RNA is not protein or binding evidence.',
        quality='Mean native log1p(full RNA), mean log1p(detected physical genes), full retained baseline TME denominator, Myeloid denominator.',
        no_effect_or_P_BH=True, no_clinical_DEG_or_signature_refit=True, no_new_download=True,
        statistical_model='Not estimated here; required controls and interface-disjoint response must be frozen before effects. The original full MCP endpoint remains primary.',
        limitations='Receiver/genes follow an already observed clinical lead, so follow-up is exploratory. No inference of sender identity from prior ranks.')
    (out/'plan_before_native_measurement.json').write_text(json.dumps(plan, indent=2)+'\n')
    g, source_rows, cells, var, eligible, measured = hnscc_inputs(root, project)
    original = pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_patient_subtype_scores.tsv', sep='\t')
    original = original[original.subtype_candidate.eq('Macro_C1QC')].set_index('patient_id').sort_index()
    assert original.index.is_unique and len(original) == 26
    baseline = cells.timepoint.eq('pre') & cells.patient_id.isin(original.index) & cells.library_id.isin(eligible.library_id)
    receiver = baseline & cells.major.eq('Myeloid') & cells.subtype.eq('Macro_C1QC') & cells.evidence.eq('multigene_subtype_candidate')
    assert not cells[baseline].groupby('patient_id').library_id.nunique().gt(1).any()
    aliases = defaultdict(set)
    for c, (name, row) in enumerate(var.iterrows()):
        if not measured[c] or str(name).startswith('AMBIGUOUS::'):
            continue
        names = {str(name), str(row.get('gene_symbol', ''))}
        if row.get('feature_identity_status') == 'stable_Ensembl_ID':
            names |= set(str(row.get('source_gene_aliases', '')).split('|'))
        for name in names - {'', 'nan', 'None'}:
            aliases[name].add(c)
    gene_columns = {gene: next(iter(aliases[gene])) if len(aliases[gene]) == 1 else None for gene in genes}
    physical = sorted({c for c in gene_columns.values() if c is not None})
    local = {c: i for i, c in enumerate(physical)}
    groups = senders.sender_group.astype(str).tolist() + ['Receiver_C1QC']
    patients = original.index.astype(str).tolist()
    keys = [(group, patient) for group in groups for patient in patients]
    key_index = {key: i for i, key in enumerate(keys)}
    n = np.zeros(len(keys), np.int64)
    sums = np.zeros((len(keys), len(physical)))
    hits = np.zeros_like(sums, dtype=np.int64)
    logRNA = np.zeros(len(keys)); loggenes = np.zeros(len(keys))
    edge_columns = []
    for r in receptor_edges.itertuples():
        chain = [gene_columns[v] for v in r.receptor.split('_')]
        edge_columns.append(chain if all(c is not None for c in chain) and len(set(chain)) == len(chain) else None)
    joint = np.zeros((len(patients), len(receptor_edges)), np.int64)
    positions = np.flatnonzero(baseline.to_numpy())
    matrix = ad.io.sparse_dataset(g['layers/counts'])
    for start in range(0, len(positions), 2000):
        pos = positions[start:start+2000]
        md = cells.iloc[pos]
        raw = matrix[source_rows[pos], :].tocsr().astype(float)
        raw.sum_duplicates(); raw.eliminate_zeros()
        assert np.isfinite(raw.data).all() and (raw.data > 0).all()
        totals = np.asarray(raw.sum(axis=1)).ravel(); detected = np.diff(raw.indptr)
        assert (totals > 0).all()
        x = raw[:, physical].toarray()
        normalized = np.log1p(x*10000/totals[:, None])
        assigned = md.major.astype(str).to_numpy().astype(object)
        rx = receiver.iloc[pos].to_numpy()
        assigned[(assigned == 'Myeloid') & ~rx] = 'Myeloid_excluding_C1QC'
        assigned[rx] = 'Receiver_C1QC'
        ids = md.patient_id.astype(str).to_numpy()
        for group in groups:
            for patient in sorted(set(ids[assigned == group])):
                use = (assigned == group) & (ids == patient)
                k = key_index[group, patient]
                n[k] += int(use.sum()); sums[k] += normalized[use].sum(axis=0)
                hits[k] += (x[use] > 0).sum(axis=0)
                logRNA[k] += np.log1p(totals[use]).sum(); loggenes[k] += np.log1p(detected[use]).sum()
                if group == 'Receiver_C1QC':
                    p = patients.index(patient)
                    for e, columns in enumerate(edge_columns):
                        if columns is not None:
                            joint[p, e] += np.all(x[use][:, [local[c] for c in columns]] > 0, axis=1).sum()
    records = []; metadata = []
    whole = cells[baseline].groupby('patient_id').size().reindex(patients)
    major = cells[baseline & cells.major.eq('Myeloid')].groupby('patient_id').size().reindex(patients, fill_value=0)
    for k, (group, patient) in enumerate(keys):
        metadata.append(dict(sender_group=group, patient_id=patient, native_cells=int(n[k]),
            mean_logRNA=logRNA[k]/n[k] if n[k] else np.nan,
            mean_loggenes=loggenes[k]/n[k] if n[k] else np.nan,
            retained_TME_cells=int(whole[patient]), Myeloid_cells=int(major[patient]),
            retained_TME_fraction=n[k]/whole[patient],
            receiver_fraction_within_Myeloid=n[k]/major[patient] if group == 'Receiver_C1QC' and major[patient] else np.nan))
        for gene in genes:
            c = gene_columns[gene]
            exists = c is not None
            records.append(dict(sender_group=group, patient_id=patient, gene=gene,
                unique_measured_column=exists, native_column=c, native_cells=int(n[k]),
                positive_cells=int(hits[k, local[c]]) if exists and n[k] else np.nan,
                positive_fraction=hits[k, local[c]]/n[k] if exists and n[k] else np.nan,
                mean_log1p_CPTT=sums[k, local[c]]/n[k] if exists and n[k] else np.nan,
                status='measured' if exists and n[k] else 'no_cells' if exists else 'unavailable_or_ambiguous_gene'))
    expression = pd.DataFrame(records); meta = pd.DataFrame(metadata)
    receptor_records = []
    for p, patient in enumerate(patients):
        count = int(n[key_index['Receiver_C1QC', patient]])
        assert count == original.loc[patient, 'native_subtype_cells']
        for e, edge in enumerate(receptor_edges.itertuples()):
            exists = edge_columns[e] is not None
            receptor_records.append(dict(patient_id=patient, ligand=edge.ligand, receptor=edge.receptor,
                native_receiver_cells=count, all_components_measured=exists,
                same_cell_positive_cells=int(joint[p, e]) if exists else np.nan,
                same_cell_positive_fraction=joint[p, e]/count if exists else np.nan,
                original_pooled_reference_and_RNA_qualified=bool(edge.prior_reference_and_RNA_qualified)))
    summaries = []
    for candidate in selected.itertuples():
        for group in senders.sender_group:
            frame = expression[expression.sender_group.eq(group) & expression.gene.eq(candidate.ligand)]
            assert len(frame) == len(patients)
            usable = frame[frame.native_cells.ge(20) & frame.status.eq('measured')]
            total = int(usable.native_cells.sum())
            summaries.append(dict(programme_sign=candidate.programme_sign, ligand=candidate.ligand, sender_group=group,
                receiver_prior_qualified=bool(qualified.loc[candidate.ligand, 'prospective_receiver_qualified']),
                source_gene_uniquely_measured=gene_columns[candidate.ligand] is not None,
                captured_matched_patients=int(frame.native_cells.ge(20).sum()), measured_matched_patients=len(usable),
                patients_with_any_positive_ligand_cells=int(usable.positive_cells.gt(0).sum()),
                patients_with_atleast1pct_positive_ligand_cells=int((usable.positive_cells*100 >= usable.native_cells).sum()),
                patients_with_atleast10pct_positive_ligand_cells=int((usable.positive_cells*10 >= usable.native_cells).sum()),
                pooled_sender_positive_fraction=float(usable.positive_cells.sum()/total) if total else np.nan,
                distinct_mean_expression_values=int(usable.mean_log1p_CPTT.nunique()),
                expression_association_not_tested=True))
    summary = pd.DataFrame(summaries)
    assert len(summary) == 400 and not summary.duplicated(['programme_sign', 'ligand', 'sender_group']).any()
    expression.to_csv(out/'private_all_patient_group_gene_expression.tsv', sep='\t', index=False, na_rep='NA')
    meta.to_csv(out/'private_all_patient_group_capture_quality.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(receptor_records).to_csv(out/'private_all_patient_receptor_complex_RNA.tsv', sep='\t', index=False, na_rep='NA')
    summary.to_csv(out/'all400_private_sender_ligand_coverage.tsv', sep='\t', index=False, na_rep='NA')
    terminal = dict(status='ALL400_FROZEN_SENDER_LIGANDS_NATIVE_BASELINE_MEASURED_NO_EFFECTS',
        finished_UTC=datetime.now(timezone.utc).isoformat(), all_native_baseline_cells=len(positions),
        native_receiver_cells=int(n[[key_index['Receiver_C1QC', p] for p in patients]].sum()),
        planned_sender_ligand_slots=len(summary), measured_gene_aliases=sum(v is not None for v in gene_columns.values()),
        unique_physical_gene_columns=len(physical), metadata_rows=len(meta), patient_gene_rows=len(expression),
        receptor_patient_records=len(receptor_records),
        coverage_sha256=sha(out/'all400_private_sender_ligand_coverage.tsv'),
        expression_sha256=sha(out/'private_all_patient_group_gene_expression.tsv'),
        quality_sha256=sha(out/'private_all_patient_group_capture_quality.tsv'),
        no_effect_P_or_threshold_selection=True, direct_axis_verified=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')
    print(json.dumps(terminal))
    print(summary.groupby('sender_group').agg(slots=('ligand', 'size'),
        nonzero_in_any_patient=('patients_with_any_positive_ligand_cells', lambda s: int(s.gt(0).sum())),
        nonzero_in_20_patients=('patients_with_any_positive_ligand_cells', lambda s: int(s.ge(20).sum())),
        receiver_qualified=('receiver_prior_qualified', 'sum')).to_json(orient='index'))


if __name__ == '__main__':
    main()
