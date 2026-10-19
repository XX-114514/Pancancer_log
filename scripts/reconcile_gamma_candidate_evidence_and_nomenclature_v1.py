"""Private identity evidence sidecar and explicit public-symbol mapping correction.

No RNA, clinical effect or formal reannotation. Original frozen pilot and literal
model qualification remain untouched; corrected measurement is a separate version.
"""
import io
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from external_functional_coverage_inputs_v1 import read_inputs
from audit_all61_fixed_gamma_scGate_native_opportunities_v1 import resolve
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 2 * 1024**2)
    fixed = root / 'NSCLC_gamma_scGate_fixed_pilot_20261019_v1'
    review = root / 'NSCLC_fixed_gamma_scGate_pilot_review_20261019_v1'
    hard = root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1'
    required = {
        'NSCLC_gamma_delta_RNA_candidate_audit_20261019_v1/private_terminal.json': 'e4d6c57c1727d33dd4d123178929125787706afa64964520c7c0ab1ff7429502',
        'NSCLC_gamma_scGate_fixed_pilot_audit_20261019_v1/private_terminal.json': '519a8671d37cc311d45c799acbd6ce70a2c76fafe956a70f209ee7fa62175cbf',
        'all61_fixed_gamma_scGate_native_opportunity_audit_20261019_v1/private_terminal.json': 'db230a445874a26c99fea20e9482b3ec0aa95df8edf540e7d2f9ab7cb3388899',
    }
    for name, expected in required.items():
        assert sha(root / name) == expected
    t = json.loads((review / 'private_terminal.json').read_text())
    assert t['status'] == 'COMPLETE_ALL_SOURCE_REFERENCE_CALIBRATION_AND_JOINT_TIERS_FIGURES_VISUAL_AND_SCIENTIFIC_REVIEW_PENDING'
    for name, info in t['files'].items():
        assert sha(review / name) == info['sha256']
    out = root / 'gamma_candidate_evidence_nomenclature_reconciliation_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    correction = dict(model_token='TRDV5', accepted_native_symbol='TRAV29DV5', IMGT_name='TRAV29/DV5', public_gene_id='28653', public_HGNC='HGNC:12127', public_Ensembl='ENSG00000211810', references=['https://www.imgt.org/IMGTrepertoire/LocusGenes/genetable/autotable.php?group=TRDV&species=Homo+sapiens', 'https://www.ncbi.nlm.nih.gov/gene/28653/'], shared_alpha_delta_variable_segment_not_gamma_specific=True, no_fuzzy_aliasing=True)
    save('plan_before_metadata_reconciliation.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), inputs=required, correction=correction, resource_guard=guard, new_RNA_read=False, no_original_artifact_mutation=True, no_new_identity_or_clinical_adoption=True))
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    unit = np.load(native / 'private_cell_unit_index.npy', allow_pickle=False)
    u = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False).sort_values('unit_index')
    assert len(unit) == 1254749 and u.unit_index.tolist() == list(range(11701))
    reference = np.zeros(len(unit), np.uint8)
    labels = u.sub_cell_type.to_numpy(str)[unit]
    reference[labels == 'T_gdT_TRDV1'] = 1
    reference[labels == 'T_gdT_TRDV2'] = 2
    h = dict(np.load(hard / 'private_candidate_identity_sidecar.npz', allow_pickle=False))
    p = dict(np.load(fixed / 'private_complete_pilot_candidate_classification_sidecar.npz', allow_pickle=False))
    bits = np.zeros(len(unit), np.uint8)
    bits[h['native_row']] = h['rule_bits']
    selected = np.union1d(np.union1d(np.flatnonzero(reference), h['native_row']), p['native_row']).astype(np.int32)
    assert np.array_equal(unit[h['native_row']], h['unit_index'])
    # Missing pilot assessment is distinct from a measured Impure classification.
    calls = np.zeros((len(selected), 4), np.uint8)
    class_codes = dict(Pure=1, Impure=2, UNMEASURABLE=3)
    assert p['models'].tolist() == ['HiTME_full', 'PBMC_exact', 'HiTME_native_subset', 'NK_HiTME_exact']
    for j in range(4):
        for text in np.unique(p['classes'][:, j]):
            assert text in class_codes
            calls[np.searchsorted(selected, p['native_row'][p['classes'][:, j] == text]), j] = class_codes[text]
    vd1 = np.zeros(len(unit), bool)
    vd2 = vd1.copy()
    vd1[h['native_row']] = h['TRDV1_detected']
    vd2[h['native_row']] = h['TRDV2_detected']
    payload = io.BytesIO()
    np.savez_compressed(payload, native_row=selected, author_reference_code=reference[selected], raw_rule_bits=bits[selected], pilot_model_codes=calls, TRDV1_detected_in_raw_candidate=vd1[selected], TRDV2_detected_in_raw_candidate=vd2[selected])
    assert len(payload.getvalue()) < 1024**2
    sidecar = out / 'private_union_reference_and_candidate_identity_evidence.npz'
    sidecar.write_bytes(payload.getvalue())
    reread = dict(np.load(sidecar, allow_pickle=False))
    np.testing.assert_array_equal(reread['native_row'], selected)
    np.testing.assert_array_equal(reread['author_reference_code'], reference[selected])
    np.testing.assert_array_equal(reread['raw_rule_bits'][np.searchsorted(selected, h['native_row'])], h['rule_bits'])
    decoder = {0: 'NOT_ASSESSED', 1: 'Pure', 2: 'Impure', 3: 'UNMEASURABLE'}
    decoded = np.array([[decoder[int(v)] for v in row] for row in calls[np.searchsorted(selected, p['native_row'])]])
    np.testing.assert_array_equal(decoded, p['classes'])
    save('private_identity_evidence_contract.json', dict(source='GSE243013', source_native_cells=len(unit), union_rows=len(selected), native_row_scope='Frozen zero-based original native row, not patient or cell barcode', authoritative_cell_unit_axis_sha256=sha(native / 'private_cell_unit_index.npy'), author_unit_metadata_sha256=sha(native / 'private_all_native_units_and_quality.tsv'), author_reference_codes={'0': 'not_original_author_gamma_reference', '1': 'T_gdT_TRDV1', '2': 'T_gdT_TRDV2'}, pilot_model_order=p['models'].tolist(), pilot_class_codes=decoder, raw_rule_bit_definition_sha256=sha(hard / 'plan_before_any_gamma_delta_RNA_counts.json'), TRDV_detection_flags_scope='Only hard-rule candidate rows are assessed; false outside those rows is not a negative gene test', rows_omitted='Outside author gamma reference, raw candidate union and uniform pilot; no pilot assessment or formal identity inference', original_author_identities_preserved=True, final_identity_accepted=False, all_selected_metadata_roundtrips_exact=True))
    # Correct one explicitly documented nomenclature alias, not a gene deletion.
    runtime = root / 'gamma_scGate_official_runtime_qualification_20261019_v1'
    folder = runtime / 'scGate_models-41a45cd3f8bb5f5a7daf21ec276f6a726f6ee0d4/human/HiTME'
    full, _ = resolve(folder, 'gdT')
    assert 'TRDV5' in full and 'TRAV29DV5' not in full
    axes = read_inputs(root)
    own = pd.read_csv(root / 'GSE243013_genes.csv.gz', keep_default_na=False).geneSymbol.str.strip().str.upper()
    assert own.eq('TRAV29DV5').sum() == 1 and own.eq('TRDV5').sum() == 0
    own_lookup = {v: {i} for i, v in enumerate(own)}
    mappings = list(axes['source_axes'].items()) + [(('original_native', 'GSE243013'), {'lookup': own_lookup})]
    records = []
    for (arm, dataset), axis in mappings:
        match = {}
        for gene in full:
            match[gene] = set(axis['lookup'].get(gene, set()))
            if gene == 'TRDV5':
                match[gene] |= set(axis['lookup'].get('TRAV29DV5', set()))
        columns = [next(iter(v)) for v in match.values() if len(v) == 1]
        missing = [k for k, v in match.items() if not v]
        ambiguous = [k for k, v in match.items() if len(v) > 1]
        collisions = len(columns) - len(set(columns))
        records.append(dict(arm=arm, dataset=dataset, model='HiTME_full_one_public_symbol_normalization_v2', fixed_features=len(full), missing_features=len(missing), ambiguous_features=len(ambiguous), column_collisions=collisions, all_fixed_features_uniquely_measurable=not missing and not ambiguous and not collisions, remaining_unmeasurable_features=json.dumps(sorted(missing + ambiguous)), source_native_stable_ID_not_newly_verified=True, RNA_scored=False))
    mapped = pd.DataFrame(records)
    assert len(mapped) == 62 and not mapped.duplicated(['arm', 'dataset']).any()
    mapped.to_csv(out / 'private_all62_full_model_symbol_normalization_qualification.tsv', sep='\t', index=False)
    save('private_pilot_scientific_review_and_visual_acceptance.json', dict(review_UTC=datetime.now(timezone.utc).isoformat(), source='GSE243013', fixed_all_slot_calibration_sha256=sha(review / 'private_complete_source_reference_calibration.tsv'), fixed_all_joint_tier_sha256=sha(review / 'private_complete_source_joint_candidate_tiers.tsv'), figures_actual_visual_reviewed=True, figure_sha256={name: sha(review / name) for name in ['source_complete_reference_calibration.png', 'source_complete_joint_candidate_tiers.png']}, visual_scope='All source-only counts/denominators, targets, missing gray slot, labels/legend/footnote visible; gray slot means original literal mapping failed, not physical absence after nomenclature correction', decision='PRESERVE_ALL_AUTHOR_REFERENCE_AND_ALL_CANDIDATE_FLAGS_NO_FULL_SOURCE_FINAL_IDENTITY_ADOPTION', reasoning=['Both integer co-detection and fixed scGate pilot recover a small fraction of author gamma references; reference labels are not TCR truth.', 'Two gamma models disagree, and NK conflict remains; do not select their small intersection as definitive truth or call noncaptured author cells false.', 'Do not retune thresholds on clinical effects. Literal TRDV5 mapping failed; canonical TRAV29DV5 is present. This is a separate versioned measurement correction, not demonstrated accuracy improvement.', 'Other-source RNA and clinical adoption remain pending; existing TCR/protein metadata and cluster-level lineage evidence are next orthogonal options.'], clinical_tests_performed=False, formal_identity_or_state_accepted=False))
    total = sum(f.stat().st_size for f in out.iterdir() if f.is_file())
    assert total < 2 * 1024**2
    receipt = dict(status='COMPLETE_PRIVATE_AUTHOR_REFERENCE_CANDIDATE_EVIDENCE_AND_VERSIONED_NOMENCLATURE_CORRECTION_NO_IDENTITY_ADOPTION', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), resource_guard=guard, source_axis_positions=62, pilot_classification_and_raw_candidate_metadata_exact=True, full_expression_read=False, public_scientific_outputs=False, files={f.name: dict(sha256=sha(f), bytes=f.stat().st_size) for f in out.iterdir() if f.is_file()})
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
