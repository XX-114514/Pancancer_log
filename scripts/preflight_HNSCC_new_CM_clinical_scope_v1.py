"""Read-only clinical scope audit for all accepted new atlas CM discoveries.

No clinical effect or NMF refit. Keep original October2026 v20261010 9CM
separate from the new atlas scopes and forbid post-treatment H as baseline.
"""
import hashlib, json, os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd


def main():
    project = Path(os.environ['PANCANCER_PROJECT_ROOT'])
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    atlas = project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    out = root/'HNSCC_new_CM_clinical_scope_20261016_v1'
    assert not out.exists()
    gate = json.loads((atlas/'manifests/immune_final_interface_audit_v1_public.json').read_text())
    assert gate['status'] == 'PASS_FINAL_ARCHIVE_INTERFACES'
    index_path = atlas/'manifests/CM_recovery_v3_public.json'
    index = json.loads(index_path.read_text())
    assert index['all_ready_scopes_completed'] and len(index['results']) == 15
    author_path = root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    author = pd.read_csv(author_path, sep='\t', keep_default_na=False)
    assert len(author) == 57 and author.library_id.is_unique
    # Canonical source reader constructs sampleID as dataset::library_id.
    author['sampleID'] = 'GSE286827::'+author.library_id.astype(str)
    assert author.sampleID.is_unique
    author = author.set_index('sampleID')
    summaries = []; matches = []
    for discovery in index['results']:
        folder = atlas/'CM'/discovery['scope']
        accepted = folder/discovery['discovery_attempt']
        meta = pd.read_csv(folder/'inputs/sample_metadata_private.tsv', sep='\t', keep_default_na=False)
        freq = pd.read_csv(folder/'inputs/frequency_private.tsv', sep='\t', index_col=0)
        W = pd.read_csv(accepted/'W_private.tsv', sep='\t', index_col=0)
        H = pd.read_csv(accepted/'H_patient_block_private.tsv', sep='\t', index_col=0)
        selected = meta[meta.dataset.eq('GSE286827')].copy()
        assert selected.sampleID.is_unique
        missing = selected[~selected.sampleID.isin(author.index)]
        if len(missing):
            raise ValueError('Exact source reader sampleID to original author library key is incomplete; no guessed aliases')
        baseline = selected[selected.timepoint.eq('pre')]
        if len(baseline):
            joined = baseline.merge(author[['library_id', 'patient_registry', 'time', 'Tumor.regression', 'neo_type', 'HPV/p16']],
                left_on='sampleID', right_index=True, validate='one_to_one')
            assert joined.time.eq('Pre').all()
            for patient, group in joined.groupby('patient_registry'):
                assert group['Tumor.regression'].nunique() == 1
                assert group['neo_type'].nunique() == 1 and group['HPV/p16'].nunique() == 1
                matches.append(dict(scope=discovery['scope'], patient_registry=patient,
                    baseline_libraries=len(group), independent_units=group.independent_unit.nunique()))
        else:
            joined = pd.DataFrame()
        assert W.index.equals(freq.index), 'NMF W node order differs from source frequencies'
        assert list(W.columns) == list(H.index), 'NMF component keys differ'
        assert set(freq.columns[1:]) == set(meta.sampleID)
        repeats = meta.groupby('independent_unit').timepoint.nunique()
        mixed = set(repeats[repeats > 1].index)
        summaries.append(dict(scope=discovery['scope'], rank=W.shape[1],
            supported_network_modules=discovery['n_modules_with_network_members'],
            HNSCC_libraries=len(selected), HNSCC_baseline_libraries=len(baseline),
            baseline_reviewed_author_patients=int(joined.patient_registry.nunique()) if len(joined) else 0,
            baseline_patients_with_mixed_time_H=int(joined[joined.independent_unit.isin(mixed)].patient_registry.nunique()) if len(joined) else 0,
            full_W_frequency_node_order_equal=True,
            W_sha256=hashlib.sha256((accepted/'W_private.tsv').read_bytes()).hexdigest(),
            input_frequency_sha256=hashlib.sha256((folder/'inputs/frequency_private.tsv').read_bytes()).hexdigest(),
            H_not_accepted_as_baseline_predictor=True))
    out.mkdir()
    pd.DataFrame(matches).to_csv(out/'private_baseline_scope_patient_alignment.tsv', sep='\t', index=False)
    pd.DataFrame(summaries).to_csv(out/'private_all15_CM_clinical_scope.tsv', sep='\t', index=False)
    report = dict(status='ALL15_NEW_ATLAS_CM_CLINICAL_SCOPE_AUDITED_NO_EFFECT_TEST',
        finished_UTC=datetime.now(timezone.utc).isoformat(),
        input_index_sha256=hashlib.sha256(index_path.read_bytes()).hexdigest(),
        original_author_table_sha256=hashlib.sha256(author_path.read_bytes()).hexdigest(),
        scopes=summaries, original_20261010_9CM_changed=False,
        next_contract='All eligible new scopes/components frozen before clinical effects. Baseline frequencies require frozen-discovery normalization and W projection; mixed-time H cannot be used as baseline clinical predictor.',
        new_model_fit_or_clinical_P=False, new_expression_download=False)
    (out/'private_terminal.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
