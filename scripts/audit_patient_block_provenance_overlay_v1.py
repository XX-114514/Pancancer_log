"""Independent row-wise evidence replay; no producer import or identity changes."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    start = datetime.now(timezone.utc).isoformat()
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    family = root / 'immune_same_name_family_source_review_v1'
    out = family / 'patient_block_provenance_overlay_v1'
    assert not (out / 'independent_audit.json').exists()
    inventory = json.loads((out / 'inventory.json').read_text())
    for name, value in inventory['outputs'].items():
        assert sha(out / name) == value
    current = family / 'family_clinical_context_review_v1/review_v2'
    parent_inputs = root / 'CM_interface_cancer_association_v1/source_unit_review_v1'
    plan = json.loads((root / 'fullgene_receiver_target_code_v1/patient_block_provenance_overlay_plan_v1.json').read_text())
    paths = dict(current_scope=current / 'all1290_family_context_scope.tsv', current_audit=current / 'independent_audit.json',
                 author_mapping=current / 'exact87_author_sample_context_mapping.tsv', parent_sidecar=parent_inputs / 'public_library_patient_tissue_processing_sidecar.tsv',
                 parent_audit=parent_inputs / 'independent_audit.json', frozen_scope=family / 'frozen1290_scope.tsv')
    for key, path in paths.items():
        assert sha(path) == inventory['input_sha256'][key] == plan['input_sha256'][key], key
    assert sha(__file__) == plan['code_sha256'][Path(__file__).name]
    original = pd.read_csv(current / 'all1290_family_context_scope.tsv', sep='\t', keep_default_na=False)
    overlay = pd.read_csv(out / 'private_all1290_patient_provenance_overlay.tsv', sep='\t', keep_default_na=False)
    pd.testing.assert_frame_equal(overlay[original.columns], original)
    parent = root / 'CM_interface_cancer_association_v1/source_unit_review_v1'
    side = pd.read_csv(parent / 'public_library_patient_tissue_processing_sidecar.tsv', sep='\t', keep_default_na=False).set_index('sampleID')
    author = pd.read_csv(current / 'exact87_author_sample_context_mapping.tsv', sep='\t', keep_default_na=False).set_index('sampleID')
    frozen = pd.read_csv(family / 'frozen1290_scope.tsv', sep='\t', keep_default_na=False).set_index('sampleID')
    counts = {}
    patient_blocks = set()
    source_blocks = {}
    mismatch = 0
    for row in overlay.itertuples():
        assert row.effective_resampling_block == row.analysis_block
        if row.patient_review_status == 'reviewed_identity':
            assert row.analysis_block == frozen.loc[row.sampleID, 'analysis_block']
            status = 'reviewed_frozen_identity'
        elif row.patient_review_status == 'reviewed_author_patient':
            assert row.analysis_block == author.loc[row.sampleID, 'analysis_block']
            status = 'reviewed_author_patient'
        elif row.unit_kind == 'reviewed_patient':
            evidence = side.loc[row.sampleID]
            assert evidence.mapping_status == 'exact_public_library_match'
            assert row.gse_id == evidence.gse_id and row.analysis_block == evidence.reviewed_analysis_block
            status = 'reviewed_parent_official_metadata'
            mismatch += 1
        else:
            status = 'unreviewed_sample_fallback'
        assert row.effective_patient_evidence_status == status
        counts[status] = counts.get(status, 0) + 1
        if row.unit_kind == 'reviewed_patient':
            assert row.effective_patient_block == row.analysis_block
            patient_blocks.add(row.analysis_block)
            source_blocks.setdefault(row.gse_id, set()).add(row.analysis_block)
        else:
            assert row.effective_patient_block == ''
    assert counts == inventory['summary']['effective_evidence_status_counts']
    assert mismatch == 194 and len(patient_blocks) == 175
    sources = pd.read_csv(out / 'source_patient_block_summary.tsv', sep='\t')
    assert set(sources.gse_id) == set(source_blocks)
    for row in sources.itertuples():
        selected = original[original.gse_id.eq(row.gse_id) & original.unit_kind.eq('reviewed_patient')]
        assert row.source_patient_blocks == len(source_blocks[row.gse_id])
        assert row.reviewed_libraries == len(selected) and row.legacy_analysis_unit_keys == selected.analysis_unit.nunique()
    audit = dict(status='PASS_ALL1290_ROW_IDENTITY_PRESERVATION_OFFICIAL_PATIENT_BLOCK_PROVENANCE_AND175_SOURCE_PATIENT_BLOCKS',
                 started_UTC=start, finished_UTC=datetime.now(timezone.utc).isoformat(),
                 auditor_sha256=sha(__file__), producer_inventory_sha256=sha(out / 'inventory.json'),
                 rows=1290, inherited_official_metadata_rows=194, reviewed_patient_libraries=402,
                 source_scoped_patient_blocks=175, original_patient_block_assignments_changed=0,
                 original_analysis_unit_keys_do_not_count_patient_blocks=True,
                 cross_GSE_independent_people_not_verified=True,
                 new_statistical_models=False, formal_states=False, canonical_Zarr_changed=False,
                 scientific_release=False, axis_validation_complete=False)
    (out / 'independent_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit), flush=True)


if __name__ == '__main__':
    main()
