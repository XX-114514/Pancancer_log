"""Outcome-free qualification of clinical keys for the current frozen RNA MCP.

Read clinical definitions and identity metadata only, never programme scores.
Patient/library rows remain private. This is not a new clinical-effect family.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 4 * 1024**2)
    scope = root/'immune_RNA_MCP_discovery_scope_20261018_v1'
    table = root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    clinical = root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1/private_patient_clinical_definition.tsv'
    partition = scope/'private_metadata_only_case_source_partition.tsv'
    assert sha(scope/'private_terminal.json') == '0f4833a07bd5cdf5e8f6806d06805777dddf78adad60a6151eca991dfdda6196'
    assert sha(table) == '9942706870518c20c87a45fed44caff6bec22a1a512b6f194dc17d406df5cd22'
    assert sha(clinical) == 'ce3c924ead635e51d322669e9b3504b5abee38a30044394eb9c4f7778f949db6'
    terminal = json.loads((scope/'private_terminal.json').read_text())
    assert sha(partition) == terminal['files'][partition.name]['sha256']
    out = root/'HNSCC_immune_RNA_MCP_clinical_metadata_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        resource_guard=guard, input_sha256={str(p.relative_to(root)):sha(p) for p in [partition, table, clinical]},
        question='Exact original author clinical keys and pretreatment coverage of current frozen RNA discovery cases.',
        no_programme_score_or_clinical_effect_read=True,
        authority='Previously reviewed author mmc2 clinical table; preserve author binary conflicts; continuous pathological regression is not RECIST.',
        selection='All current selected HNSCC cases, then explicit pre only. Never substitute a post case or add a previously unselected baseline silently.',
        limits='HNSCC contributed RNA discovery. Outcomes were not used in RNA fitting, but the project has inspected them previously. Future associations are exploratory, not independent clinical validation.')
    (out/'plan_before_clinical_metadata_join.json').write_text(json.dumps(plan, indent=2)+'\n')
    author = pd.read_csv(table, sep='\t', keep_default_na=False)
    fields = ['library_id', 'patient_registry', 'timepoint_registry', 'Tumor.regression', 'neo_type', 'HPV/p16', 'Partial.pathologic.response']
    author = author[fields].copy()
    assert len(author) == 57 and author.library_id.is_unique
    author['sampleID'] = 'GSE286827::'+author.library_id.astype(str)
    chosen = pd.read_csv(partition, sep='\t', keep_default_na=False)
    chosen = chosen[chosen.dataset.eq('GSE286827')].copy()
    assert len(chosen) == 28 and chosen.case_block.is_unique and chosen.fit_sample.is_unique
    assert chosen.reviewed_patient.astype(str).str.lower().eq('true').all()
    joined = chosen.merge(author.drop(columns='library_id'), on='sampleID', how='left', validate='one_to_one', indicator='clinical_join')
    assert joined.clinical_join.eq('both').all()
    # A checked one-to-one patient bridge; do not invent string aliases.
    assert joined.groupby('patient_block').patient_registry.nunique().eq(1).all()
    assert joined.groupby('patient_registry').patient_block.nunique().eq(1).all()
    assert joined.patient_registry.is_unique
    for _, group in author.groupby('patient_registry'):
        for col in ['Tumor.regression', 'neo_type', 'HPV/p16']:
            assert group[col].nunique(dropna=False) == 1
    old = pd.read_csv(clinical, sep='\t', keep_default_na=False).set_index('patient_id')
    assert old.index.is_unique
    c = old.loc[joined.patient_registry]
    regression = pd.to_numeric(joined['Tumor.regression'], errors='raise').to_numpy(float)
    assert np.isfinite(regression).all() and ((regression>=0)&(regression<=100)).all()
    assert np.array_equal(regression, c.regression.to_numpy(float))
    assert np.array_equal(joined.neo_type.to_numpy(str), c.treatment.to_numpy(str))
    assert np.array_equal(joined['HPV/p16'].to_numpy(str), c.HPV.to_numpy(str))
    joined['regression'] = regression
    joined['derived_response'] = (regression>=50).astype(int)
    joined['author_binary_consistent'] = c.author_binary_consistent.to_numpy()
    joined['author_binary_concordant'] = c.author_binary_concordant.to_numpy()
    joined['baseline_eligible'] = joined.phase_class.eq('pre')
    baseline = joined[joined.baseline_eligible].copy()
    assert baseline.patient_registry.is_unique
    selected_fields = ['sampleID','fit_sample','patient_registry','patient_block','phase_class','timepoint_registry','partition','regression','neo_type','HPV/p16','derived_response','author_binary_consistent','author_binary_concordant','baseline_eligible']
    joined[selected_fields].to_csv(out/'private_all_selected_case_clinical_keys.tsv', sep='\t', index=False)
    baseline[selected_fields].to_csv(out/'private_baseline_case_clinical_keys.tsv', sep='\t', index=False)
    def counts(series):
        return {str(k):int(v) for k,v in series.value_counts(dropna=False).items()}
    summary = dict(dataset='GSE286827', selected_reviewed_cases=len(joined), baseline_reviewed_cases=len(baseline),
        selected_phase_counts=counts(joined.phase_class), baseline_treatment_counts=counts(baseline.neo_type),
        baseline_HPV_counts=counts(baseline['HPV/p16']), baseline_derived_pathology_groups=counts(baseline.derived_response),
        baseline_regression_nonmissing=int(np.isfinite(baseline.regression).sum()),
        baseline_regression_unique_values=int(baseline.regression.nunique()),
        baseline_author_binary_consistent=int(baseline.author_binary_consistent.astype(str).str.lower().eq('true').sum()),
        baseline_author_binary_concordant=int(baseline.author_binary_concordant.astype(str).str.lower().eq('true').sum()),
        missing_author_library_keys=0, duplicate_selected_patients=0, no_programme_scores_or_effects_computed=True,
        clinical_independent_validation=False, outcome='continuous pathological regression; author binary retained separately')
    (out/'private_source_metadata_qualification.json').write_text(json.dumps(summary, indent=2)+'\n')
    files = {p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    assert sum(v['bytes'] for v in files.values()) < 4*1024**2
    t = dict(status='COMPLETE_CURRENT_RNA_CASE_CLINICAL_METADATA_ONLY_INDEPENDENT_KEY_AUDIT_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n')
    print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__ == '__main__':
    main()
