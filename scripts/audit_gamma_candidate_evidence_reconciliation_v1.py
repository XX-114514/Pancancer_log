"""Complete metadata/flag audit of reference-preserving candidate reconciliation.

Reuses previously accepted native axes; does not claim independent RNA or TCR
truth. The scientific decision remains separate from exact metadata transfer.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from external_functional_coverage_inputs_v1 import read_inputs
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 1024**2)
    src = root / 'gamma_candidate_evidence_nomenclature_reconciliation_20261019_v1'
    out = root / 'gamma_candidate_evidence_reconciliation_audit_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    terminal = wait(root, out, 'reconcile_gamma_candidate_evidence_and_nomenclature_v1', src / 'private_terminal.json', {'COMPLETE_PRIVATE_AUTHOR_REFERENCE_CANDIDATE_EVIDENCE_AND_VERSIONED_NOMENCLATURE_CORRECTION_NO_IDENTITY_ADOPTION'})
    for name, info in terminal['files'].items():
        assert sha(src / name) == info['sha256']
    z = dict(np.load(src / 'private_union_reference_and_candidate_identity_evidence.npz', allow_pickle=False))
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    u = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False).set_index('unit_index')
    unit = np.load(native / 'private_cell_unit_index.npy', allow_pickle=False)
    types = u.sub_cell_type.to_numpy(str)[unit]
    reference = np.select([types == 'T_gdT_TRDV1', types == 'T_gdT_TRDV2'], [1, 2], default=0).astype(np.uint8)
    h = dict(np.load(root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1/private_candidate_identity_sidecar.npz', allow_pickle=False))
    p = dict(np.load(root / 'NSCLC_gamma_scGate_fixed_pilot_20261019_v1/private_complete_pilot_candidate_classification_sidecar.npz', allow_pickle=False))
    rows = z['native_row']
    expected_rows = sorted(set(np.where(reference != 0)[0]) | set(h['native_row']) | set(p['native_row']))
    assert rows.tolist() == expected_rows
    np.testing.assert_array_equal(z['author_reference_code'], reference[rows])
    expected_bits = np.zeros(len(unit), np.uint8)
    expected_bits[h['native_row']] = h['rule_bits']
    np.testing.assert_array_equal(z['raw_rule_bits'], expected_bits[rows])
    for name in ['TRDV1_detected', 'TRDV2_detected']:
        expected = np.zeros(len(unit), bool)
        expected[h['native_row']] = h[name]
        np.testing.assert_array_equal(z[name + '_in_raw_candidate'], expected[rows])
    decoder = {1: 'Pure', 2: 'Impure', 3: 'UNMEASURABLE'}
    for j, model in enumerate(p['models']):
        got = z['pilot_model_codes'][np.searchsorted(rows, p['native_row']), j]
        assert np.array([decoder[int(v)] for v in got]).tolist() == p['classes'][:, j].tolist()
    not_pilot = ~np.isin(rows, p['native_row'])
    assert not z['pilot_model_codes'][not_pilot].any()
    axes = read_inputs(root)
    fixed_plan = json.loads((root / 'all61_fixed_gamma_scGate_native_opportunities_20261019_v1/plan_before_all_source_model_axis_qualification.json').read_text())
    genes = fixed_plan['models']['HiTME_full']
    assert 'TRDV5' in genes and 'TRAV29DV5' not in genes
    own = pd.read_csv(root / 'GSE243013_genes.csv.gz', keep_default_na=False).geneSymbol.str.strip().str.upper()
    own_axis = {v: {i} for i, v in enumerate(own)}
    source_maps = {k: v['lookup'] for k, v in axes['source_axes'].items()}
    source_maps[('original_native', 'GSE243013')] = own_axis
    tab = pd.read_csv(src / 'private_all62_full_model_symbol_normalization_qualification.tsv', sep='\t', keep_default_na=False)
    assert len(tab) == len(source_maps) == 62 and set(zip(tab.arm, tab.dataset)) == set(source_maps)
    for row in tab.itertuples():
        axis = source_maps[(row.arm, row.dataset)]
        resolved = []
        missing = []
        ambiguous = []
        for gene in genes:
            names = [gene, 'TRAV29DV5'] if gene == 'TRDV5' else [gene]
            hits = {position for name in names for position in axis.get(name, set())}
            if not hits:
                missing.append(gene)
            elif len(hits) > 1:
                ambiguous.append(gene)
            else:
                resolved.append(next(iter(hits)))
        collisions = len(resolved) - len(set(resolved))
        assert row.fixed_features == len(genes)
        assert row.missing_features == len(missing) and row.ambiguous_features == len(ambiguous) and row.column_collisions == collisions
        assert row.all_fixed_features_uniquely_measurable == (not missing and not ambiguous and not collisions)
        assert json.loads(row.remaining_unmeasurable_features) == sorted(missing + ambiguous)
        assert row.source_native_stable_ID_not_newly_verified and not row.RNA_scored
    # A complete source-only descriptive record, not a selected clinical result.
    groups = ['not_author_gamma', 'author_Vdelta1', 'author_Vdelta2']
    summary = []
    for code, name in enumerate(groups):
        fullmask = reference == code
        selected_mask = reference[rows] == code
        summary.append(dict(source='GSE243013', author_reference_group=name, native_cells=int(fullmask.sum()), evidence_union_cells=int(selected_mask.sum()), any_raw_rule_candidate_cells=int((fullmask & (expected_bits > 0)).sum()), pilot_assessed_cells=int(np.isin(rows[selected_mask], p['native_row']).sum()), independent_TCR_truth=False))
    pd.DataFrame(summary).to_csv(out / 'private_complete_source_reference_candidate_counts.tsv', sep='\t', index=False)
    receipt = dict(status='ACCEPTED_ALL_REFERENCE_RAW_AND_PILOT_FLAGS_ALL62_SYMBOL_MAPS_METADATA_ONLY_NO_FORMAL_IDENTITY', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), all_candidate_metadata_rows_exact=True, all62_mapping_positions_exact=True, reference_authority_mutated=False, full_RNA_or_TCR_replayed=False, scientific_summary_public=False, resource_guard=guard, complete_source_summary_sha256=sha(out / 'private_complete_source_reference_candidate_counts.tsv'))
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    assert sum(f.stat().st_size for f in out.iterdir() if f.is_file()) < 1024**2
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
