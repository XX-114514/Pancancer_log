"""Receipt of actual review of accepted source-level F01 background summaries.

No refit, new P, gene weights/lists, case rows or new association results.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    resource(root, 1024**2)
    src = root / 'NSCLC_native_MFA_interpretation_20261016_v1'
    proof = src / 'private_independent_interpretation_audit.json'
    assert sha(proof) == '1ae54e6dcc3034d8fd40b6a12cff24784c987c563a8ef5b132f974bd53e73b8c'
    assert json.loads(proof.read_text())['producer_terminal_sha256'] == sha(src / 'private_terminal.json')
    manifest = json.loads((src / 'private_terminal.json').read_text())
    tables = {}
    for name in ['private_all40_subtype_mixture_predictions.tsv',
                 'private_all_covariate_diagnostics.tsv',
                 'private_all_author_subtype_diagnostics.tsv',
                 'private_all60_composition_and_sex_sensitivity.tsv']:
        path = src / name
        assert sha(path) == manifest['files'][name]['sha256']
        data = pd.read_csv(path, sep='\t')
        selected = data[data.factor.eq('NSCLC_MFA_20261016_F01')]
        if name == 'private_all_covariate_diagnostics.tsv':
            selected = selected[selected.covariate.isin([
                'source_age', 'source_female', 'global_logN',
                'global_mean_logRNA', 'global_mean_loggenes'])]
        tables[name] = dict(source_sha256=sha(path), aggregate_rows=len(selected),
            reviewed_aggregates=json.loads(selected.to_json(orient='records', double_precision=15)))
    report = dict(status='ACTUALLY_REVIEWED_ACCEPTED_F01_NATIVE_COMPOSITION_SEX_AND_SUBTYPE_SUMMARIES',
        recorded_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__),
        original_interpretation_audit_sha256=sha(proof), dataset='GSE243013',
        scope='Four original views, all six original pair comparisons and all original subtype rows for F01; five specified background diagnostics per view. No learned gene weights or case data.',
        tables=tables,
        interpretation_limits=['Descriptive prior-result review, not new inference or clinical validation.',
            'A shared native factor and an externally supported pair do not establish pair-specific signalling.',
            'Author subtype labels and programmes are different evidence levels; no formal state identity assigned.',
            'Original post-neoadjuvant surgical RNA cannot establish pretreatment predictive performance.',
            'Published spatial source qualification remains separate; no expression download or new spatial model in this review.'],
        new_P_models_or_clinical_outcomes_loaded=False,
        new_native_interface_association_effects_loaded=False)
    path = src / 'private_F01_background_actual_review_20261018_v1.json'
    assert not path.exists()
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(status=report['status'], recorded_UTC=report['recorded_UTC'],
                          receipt_sha256=sha(path))))


if __name__ == '__main__':
    main()
