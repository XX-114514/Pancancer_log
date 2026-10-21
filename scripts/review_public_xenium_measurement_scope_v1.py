"""Source-level scientific acceptance of joint spatial measurement limits."""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verified(path):
    t = json.loads((path / 'private_terminal.json').read_text())
    for name, value in t['files'].items():
        assert sha(path / name) == value['SHA256']
    return t


def main():
    joint = Path(os.environ['PANCANCER_JOINT_AUDIT'])
    units = Path(os.environ['PANCANCER_SOURCE_UNIT_AUDIT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    j, u = verified(joint), verified(units)
    assert j['all18_Counter_keys_and_finite_coordinates_replayed']
    assert u['independent_XML_parser_replayed_source_names_titles_case_runs_and_archive_keys']
    rows = json.loads((joint / 'private_all18_independent_joint_measurement_acceptance.json').read_text())
    source = json.loads((units / 'private_independent_source_units_and_stage_summary.json').read_text())
    assert len(rows) == len({r['archive_sample_ordinal'] for r in rows}) == 18
    assert all(r['status'] == 'ACCEPTED_EXACT_MATRIX_COORDINATE_CELL_KEY_MULTISET' and r['exact_cell_key_alignment'] and r['finite_coordinates'] for r in rows)
    panel = json.loads((joint / 'private_whitelist_gene_panel_coverage.json').read_text())
    for r in panel:
        assert r['measured_archive_files'] == sum(e['whitelisted_gene_copies'][r['gene']] == 1 for e in rows)
    target_unmeasured = all(r['CXCL13_known_alias_matches'] == r['CXCL13_Ensembl_id_matches'] == r['CXCL13_prefix_name_candidates'] == 0 for r in rows)
    marker_set = ['CD3D','CD3E','TRAC','CD4','CD8A','FOXP3','IL2RA','CTLA4','TNFRSF9','LAG3','PDCD1','TOX','HAVCR2','ENTPD1','MKI67','TOP2A','UBE2C','PCNA','TCF7','IL7R','EPCAM','SOX2']
    measured = {g: sum(r['whitelisted_gene_copies'][g] == 1 for r in rows) for g in marker_set}
    assert source['source_case_diagnoses']['HNSCC_or_squamous'] == 10 and source['source_sample_diagnoses']['HNSCC_or_squamous'] == 17
    assert source['source_case_diagnoses']['ameloblastoma'] == source['source_sample_diagnoses']['ameloblastoma'] == 1
    review = dict(status='ACCEPTED_SOURCE_JOINT_MEASUREMENT_SCOPE_NOT_CELL_STATE_OR_SPATIAL_RESULT', completed_UTC=datetime.now(timezone.utc).isoformat(), source='GSE300147', source_SHA256=sha(__file__), joint_acceptance_SHA256=sha(joint / 'private_terminal.json'), units_acceptance_SHA256=sha(units / 'private_terminal.json'), all18_matrix_coordinate_keys_and_finite_coordinates_accepted=True, original_short_transport_failure_recovered_without_dropping_any_slot=True, source_cohort=source, per_marker_measured_archive_files=measured, CXCL13_unmeasured_by_symbol_official_alias_ID_and_prefix_checks=target_unmeasured, field_presence_not_expression_detection_or_coexpression=True, no_per_cell_annotation_or_mechanism_adoption=True, scientific_decision='Do not claim direct measurement of CXCL13-positive Tex. Measured regulatory/cycling/dysfunction marker combinations are separate operational measurement questions, requiring actual coexpression and quality checks before any spatial inference.', next_action='Freeze observed marker definitions and capture-quality contrasts. Quantify all source states and eligibility before estimating state-specific spatial arrangement, retain technical replication within source cases, distinguish recurrence from anatomical tumor class and ICI timepoint.', source_case_units_not_cross_study_patient_review=True, p16_unknown_category_includes_unrecognized_literals_not_imputed=True, no_pathology_tumor_mask_TLS_or_ICI_endpoint_accepted=True, source_site_and_therapy_phase_still_need_explicit_scope=True, individual_source_files_cells_or_patient_rows_not_exported=True)
    out.mkdir(mode=0o700)
    (out / 'private_source_measurement_scope_review.json').write_text(json.dumps(review, indent=2)+'\n')
    terminal = dict(status=review['status'], completed_UTC=review['completed_UTC'], source_SHA256=review['source_SHA256'], scientific_review_SHA256=sha(out / 'private_source_measurement_scope_review.json'), no_new_spatial_or_clinical_inference=True)
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'], target_gene_unmeasured=target_unmeasured, all18_input_opportunities_accepted=True)))


if __name__ == '__main__':
    main()
