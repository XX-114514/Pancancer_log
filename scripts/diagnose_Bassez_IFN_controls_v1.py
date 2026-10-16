"""Post-result full control decomposition; no new P/BH or threshold changes."""
import json, os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from test_Bassez_IFN_validation_v1 import fitted, CLINICAL, COUPLING


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs = root/'Bassez_IFN_validation_preparation_20261016_v1'
    result = root/'Bassez_IFN_validation_results_20261016_v1'
    out = root/'Bassez_IFN_control_diagnostics_20261016_v1'
    assert not out.exists()
    frame = pd.read_csv(inputs/'private_patient_scores_and_covariates.tsv', sep='\t')
    plan = json.loads((result/'execution_contract_before_effects.json').read_text())
    levels = plan['categorical_levels']
    rows = []
    quality = ['T_logRNA', 'T_loggenes', 'My_logRNA', 'My_loggenes']
    counts = ['T_logcells', 'My_logcells']
    for question in plan['questions']:
        variants = [('unadjusted', [], False), ('treatment_cancer_subtype', [], True)]
        if question['binary']:
            variants += [('quality_only', CLINICAL[:2], False), ('count_only', CLINICAL[2:], False),
                ('treatment_cancer_subtype_quality', CLINICAL[:2], True), ('treatment_cancer_subtype_count', CLINICAL[2:], True),
                ('full_primary', CLINICAL, True)]
        else:
            variants += [('quality_only', quality, False), ('count_only', counts, False), ('receiver_IFNG_only', ['My_IFNG'], False),
                ('treatment_cancer_subtype_quality', quality, True), ('treatment_cancer_subtype_count', counts, True),
                ('treatment_cancer_subtype_receiver_IFNG', ['My_IFNG'], True),
                ('full_without_receiver_IFNG', quality+counts, True), ('full_primary', COUPLING, True)]
        for name, numeric, categorical in variants:
            q = dict(question, numeric=numeric)
            cats = levels if categorical else {}
            point = np.nan; loo = []; status = 'rank_or_df_insufficient'
            try:
                point = fitted(frame, q, cats, True)[0]; status = 'descriptive_estimated'
            except (ValueError, np.linalg.LinAlgError):
                pass
            for i in range(len(frame)):
                try:
                    loo.append(fitted(frame.drop(index=frame.index[i]), q, cats, True)[0])
                except (ValueError, np.linalg.LinAlgError):
                    pass
            rows.append(dict(question=q['name'], variant=name, numeric_controls='|'.join(numeric),
                categorical_controls=categorical, n_patients=len(frame), beta=point,
                LOO_valid=len(loo), LOO_min=min(loo) if loo else np.nan, LOO_max=max(loo) if loo else np.nan, status=status))
    assert len(rows) == 24
    correlations = []
    for endpoint in ['T_IFNG', 'My_HALLMARK_INTERFERON_GAMMA_RESPONSE']:
        for control in COUPLING:
            correlations.append(dict(readout=endpoint, control=control, rho=float(spearmanr(frame[endpoint], frame[control]).statistic)))
    source_units = pd.read_csv(root/'Bassez_fixed_function_readiness_20261014_v1/private_patient_phase_capture.tsv', sep='\t', keep_default_na=False)
    concordance = {}
    for field in ['outcome', 'Treatment.or.Mode.of.Action', 'Cancer_type_update']:
        count = 0
        for patient in frame.patient_token:
            group = source_units[source_units.patient_token.eq(patient)]
            count += group[field].nunique() > 1
        concordance[field] = int(count)
    assert not any(concordance.values()), 'Frozen outcome/treatment/subtype differs across visits'
    out.mkdir()
    pd.DataFrame(rows).to_csv(out/'all24_private_control_and_LOO_descriptors.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(correlations).to_csv(out/'private_readout_control_correlations.tsv', sep='\t', index=False)
    report = dict(status='ALL24_BASSEZ_CONTROL_DESCRIPTORS_AND_CROSSVISIT_KEYS_COMPLETE_NO_NEW_P',
        finished_UTC=datetime.now(timezone.utc).isoformat(),
        cancer_subtype_counts={str(k): int(v) for k, v in frame.cancer_type.value_counts().items()},
        treatment_counts={str(k): int(v) for k, v in frame.treatment.value_counts().items()},
        crossvisit_metadata_conflicts=concordance, all_LOO_rank_failures_retained=True,
        bootstrap_interpretation='Original E/NE-stratified bootstrap may omit a rare fixed categorical level. Original invalid counts and unavailable CI are retained; no replacement family.',
        new_P_or_BH=False)
    (out/'private_terminal.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report)); print(pd.DataFrame(rows).to_json(orient='records'))


if __name__ == '__main__':
    main()
