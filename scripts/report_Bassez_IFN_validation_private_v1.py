"""Verify all control descriptors independently and combine private evidence."""
import hashlib, json, os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd
from audit_Bassez_IFN_validation_v1 import fit, close


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs = root/'Bassez_IFN_validation_preparation_20261016_v1'
    result = root/'Bassez_IFN_validation_results_20261016_v1'
    diagnostic = root/'Bassez_IFN_control_diagnostics_20261016_v1'
    out = root/'Bassez_IFN_private_evidence_20261016_v1'
    assert not out.exists()
    audit = json.loads((result/'private_independent_audit.json').read_text())
    assert audit['status'] == 'PASS_ALL8_BASSEZ_IFN_EFFECTS_PERMUTATIONS_BOOTSTRAPS_AND_BH'
    primary_file = result/'all8_private_IFN_validation_statistics.tsv'
    assert hashlib.sha256(primary_file.read_bytes()).hexdigest() == audit['result_sha256']
    frame = pd.read_csv(inputs/'private_patient_scores_and_covariates.tsv', sep='\t')
    plan = json.loads((result/'execution_contract_before_effects.json').read_text())
    questions = {q['name']: q for q in plan['questions']}
    descriptors = pd.read_csv(diagnostic/'all24_private_control_and_LOO_descriptors.tsv', sep='\t', keep_default_na=False)
    assert len(descriptors) == 24
    for row in descriptors.itertuples():
        q = dict(questions[row.question], numeric=row.numeric_controls.split('|') if row.numeric_controls else [])
        levels = plan['categorical_levels'] if str(row.categorical_controls).lower() == 'true' else {}
        try:
            point = fit(frame, q, levels, True)[0]
            close(point, float(row.beta))
        except (ValueError, np.linalg.LinAlgError):
            assert row.beta in ['NA', '']
        leave = []
        for i in range(len(frame)):
            try:
                leave.append(fit(frame.drop(index=frame.index[i]), q, levels, True)[0])
            except (ValueError, np.linalg.LinAlgError):
                pass
        assert len(leave) == row.LOO_valid
        if leave:
            close(min(leave), float(row.LOO_min)); close(max(leave), float(row.LOO_max))
    out.mkdir()
    text = '''# Independent Bassez IFN validation (private)

## Prespecified scope and interpretation

All eligible native baseline author T/myeloid cells and original E/NE labels
were retained. E denotes T-clonal expansion, not pathology/RECIST or survival.
These are fixed public partial functional readouts with source-positive
measurement evidence; this does not rescue incomplete old MCP transfer or
establish an accepted old172/C1QC identity mapping.

Three primary questions were frozen before scoring/effects. Adjusted models
include treatment, cancer subtype and required RNA/count controls. The sender
model additionally controls receiver own IFNG. Interface columns are excluded
from both response signatures and rank backgrounds. All primary slots remain
in BH3; unadjusted BH3 and treatment-stratified Wilcoxon BH2 are sensitivities.

The rank-biserial Wilcoxon contrast is weighted within treatment strata and
permuted within those strata. Both positive and negative results are retained;
it is not a replacement for the fully adjusted model.

## Original effect families and independent numerical audit

'''
    statistics = pd.read_csv(primary_file, sep='\t')
    text += '```text\n'+statistics.to_string(index=False)+'\n```\n'
    text += '\n```json\n'+json.dumps(audit, indent=2)+'\n```\n'
    text += '''

## Post-result diagnostics: no new P or BH

All prespecified questions are decomposed into full control blocks and all
leave-one-patient fits. Removing a rare categorical level can make the fixed
design unidentifiable; failed fits stay unavailable. Original bootstrap has
the same limitation. Missing CI is not evidence of a narrow/zero effect; no
replacement bootstrap family is created after observing results.

'''
    text += '```text\n'+descriptors.to_string(index=False)+'\n```\n'
    for filename in ['private_terminal.json', 'private_readout_control_correlations.tsv']:
        text += '\n### '+filename+'\n\n```text\n'+(diagnostic/filename).read_text()+'\n```\n'
    text += '\n## Measurement and compatible receptor RNA\n\n```json\n'+(inputs/'private_terminal.json').read_text()+'\n```\n'
    text += '\n```json\n'+(result/'private_descriptive_interface_and_readout.json').read_text()+'\n```\n'
    next_scope = root/'HNSCC_new_CM_clinical_scope_20261016_v1/private_terminal.json'
    if next_scope.exists():
        text += '\n## Next module clinical input audit; no effect test\n\n```json\n'+next_scope.read_text()+'\n```\n'
    text += '''

## Decision rules

No unadjusted or single sensitivity P can substitute for primary evidence.
Control decomposition is descriptive and cannot prove whether shared RNA
complexity is purely technical or biological. Receptor coexpression provides
RNA compatibility only. Neither score correlation nor an E/NE difference
demonstrates direct binding, a causal sender or clinical prediction.

Preserve the earlier HNSCC single-end lead and all current counterevidence.
Do not keep varying IFNG thresholds to obtain significance. Next assess
capturable multicellular module clinical/cohort value using frozen weights
and baseline-only composition; patient-block H that mixes pre/post cannot
be used as a baseline predictor. New atlas modules remain distinct from
October2026 v20261010 9CM. Any new clinical family must be frozen before effects.

The full biological/clinical cell-interaction objective remains active.
'''
    (out/'BASSEZ_IFN_PRIVATE_EVIDENCE.md').write_text(text)
    terminal = dict(status='BASSEZ_ALL8_TESTS_AND_ALL24_CONTROL_LOO_REPLAY_EVIDENCE_COMPLETE_PRIVATE',
        finished_UTC=datetime.now(timezone.utc).isoformat(),
        all24_control_descriptors_independent_manual_rank_QR=True,
        primary_result_sha256=audit['result_sha256'],
        report_sha256=hashlib.sha256((out/'BASSEZ_IFN_PRIVATE_EVIDENCE.md').read_bytes()).hexdigest(),
        new_significance_tests_in_report=False, new_expression_download=False,
        direct_interaction_verified=False, goal_completed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')
    print(json.dumps(terminal))


if __name__ == '__main__':
    main()
