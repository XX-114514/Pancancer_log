"""New metadata overlay; original units, identities and frozen tables stay intact.

Correct effective evidence labels only after matching audited primary sidecars.
Resampling uses original audited analysis_block, not stale per-library unit keys.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    started = datetime.now(timezone.utc).isoformat()
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    family = root / 'immune_same_name_family_source_review_v1'
    out = family / 'patient_block_provenance_overlay_v1'
    assert not out.exists()
    current = family / 'family_clinical_context_review_v1/review_v2'
    parent = root / 'CM_interface_cancer_association_v1/source_unit_review_v1'
    paths = dict(current_scope=current / 'all1290_family_context_scope.tsv',
                 current_audit=current / 'independent_audit.json',
                 author_mapping=current / 'exact87_author_sample_context_mapping.tsv',
                 parent_sidecar=parent / 'public_library_patient_tissue_processing_sidecar.tsv',
                 parent_audit=parent / 'independent_audit.json',
                 frozen_scope=family / 'frozen1290_scope.tsv')
    plan = json.loads((root / 'fullgene_receiver_target_code_v1/patient_block_provenance_overlay_plan_v1.json').read_text())
    for key, path in paths.items():
        assert sha(path) == plan['input_sha256'][key], key
    assert sha(__file__) == plan['code_sha256'][Path(__file__).name]
    parent_audit = json.loads(paths['parent_audit'].read_text())
    assert parent_audit['status'] == 'PASS_ALL265_PUBLIC_METADATA_SIDE_KEYS_EXACT_GSM_LIBRARY_PATIENT_SPECIMEN_TISSUE_PROCESSING_AND_MMR_VALUES'
    assert parent_audit['n_exact_matches'] == 194 and parent_audit['n_unresolved_US_aliases'] == 71
    current_audit = json.loads(paths['current_audit'].read_text())
    assert current_audit['scope']['reviewed_samples'] == 87
    scope = pd.read_csv(paths['current_scope'], sep='\t', keep_default_na=False)
    sidecar = pd.read_csv(paths['parent_sidecar'], sep='\t', keep_default_na=False)
    author = pd.read_csv(paths['author_mapping'], sep='\t', keep_default_na=False)
    frozen = pd.read_csv(paths['frozen_scope'], sep='\t', keep_default_na=False)
    assert len(scope) == len(frozen) == 1290 and scope.sampleID.is_unique
    assert scope.sampleID.tolist() == frozen.sampleID.tolist()
    assert sidecar.sampleID.is_unique and author.sampleID.is_unique
    exact = sidecar[sidecar.mapping_status.eq('exact_public_library_match')].set_index('sampleID')
    lookup = author.set_index('sampleID')
    inherited = scope.unit_kind.eq('reviewed_patient') & scope.patient_review_status.eq('legacy_unreviewed')
    assert int(inherited.sum()) == 194
    effective = []
    for row in scope.itertuples():
        if row.patient_review_status == 'reviewed_identity':
            original = frozen[frozen.sampleID.eq(row.sampleID)].iloc[0]
            assert row.analysis_block == original.analysis_block and row.unit_kind == 'reviewed_patient'
            value = 'reviewed_frozen_identity'
        elif row.patient_review_status == 'reviewed_author_patient':
            evidence = lookup.loc[row.sampleID]
            assert row.gse_id == evidence.gse_id and row.analysis_block == evidence.analysis_block
            assert row.analysis_unit == evidence.analysis_unit and row.unit_kind == 'reviewed_patient'
            value = 'reviewed_author_patient'
        elif row.unit_kind == 'reviewed_patient':
            evidence = exact.loc[row.sampleID]
            assert row.gse_id == evidence.gse_id and row.analysis_block == evidence.reviewed_analysis_block
            assert row.gse_id in ['GSE178341', 'GSE253013']
            value = 'reviewed_parent_official_metadata'
        else:
            assert row.unit_kind == 'sample_fallback'
            value = 'unreviewed_sample_fallback'
        effective.append(value)
    overlay = scope.copy()
    overlay['effective_patient_evidence_status'] = effective
    overlay['effective_resampling_block'] = overlay.analysis_block
    overlay['effective_patient_block'] = overlay.analysis_block.where(overlay.unit_kind.eq('reviewed_patient'), '')
    overlay['original_analysis_unit_is_not_always_patient_key'] = overlay.analysis_unit.ne(overlay.analysis_block.str.removeprefix('patient::')) & inherited
    reviewed = overlay[overlay.unit_kind.eq('reviewed_patient')]
    summary = dict(libraries=1290, reviewed_patient_libraries=len(reviewed),
                   sample_fallback_libraries=int(overlay.unit_kind.eq('sample_fallback').sum()),
                   source_scoped_patient_blocks=int(reviewed.effective_resampling_block.nunique()),
                   effective_evidence_status_counts=overlay.effective_patient_evidence_status.value_counts().to_dict(),
                   inherited_review_status_labels_corrected=int(inherited.sum()),
                   original_patient_block_assignments_changed=0,
                   independent_cross_GSE_people_not_verified=True)
    assert summary['reviewed_patient_libraries'] == 402 and summary['source_scoped_patient_blocks'] == 175
    sources = reviewed.groupby('gse_id').agg(reviewed_libraries=('sampleID', 'size'),
                                           source_patient_blocks=('effective_resampling_block', 'nunique'),
                                           legacy_analysis_unit_keys=('analysis_unit', 'nunique')).reset_index()
    out.mkdir()
    overlay.to_csv(out / 'private_all1290_patient_provenance_overlay.tsv', sep='\t', index=False)
    sources.to_csv(out / 'source_patient_block_summary.tsv', sep='\t', index=False)
    manifest = dict(status='METADATA_OVERLAY_COMPLETE_PENDING_INDEPENDENT_AUDIT', started_UTC=started,
                    finished_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__),
                    input_sha256={name: sha(path) for name, path in paths.items()}, summary=summary,
                    outputs={p.name: sha(p) for p in out.iterdir() if p.is_file()},
                    new_statistical_models=False, formal_states=False, canonical_Zarr_changed=False,
                    scientific_release=False, axis_validation_complete=False)
    (out / 'inventory.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(dict(status=manifest['status'], summary=summary)), flush=True)


if __name__ == '__main__':
    main()
