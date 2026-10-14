"""Finalize completed v1 readiness artifacts; repair NumPy JSON scalar only.

The raw-count scan already finished. Reuse its observed column set and all
coverage/capture tables; do not rerun scans, alter gates or fit effects.
"""
import hashlib, json, os
from datetime import datetime, timezone
from pathlib import Path
import h5py, numpy as np, pandas as pd
from anndata._io.specs import read_elem


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    project = Path(os.environ['PANCANCER_PROJECT_ROOT'])
    prior = root / 'Bassez_fixed_function_readiness_20261014_v1'
    out = root / 'Bassez_fixed_function_readiness_20261014_v2'
    assert not out.exists() and not (prior / 'private_terminal.json').exists()
    files = ['plan_before_expression_review.json', 'private_all50_measurement_coverage.tsv',
             'private_reference_columns.tsv', 'private_patient_phase_capture.tsv',
             'private_source_positive_confirmed_columns.npy']
    manifest = []
    for name in files:
        p = prior / name
        assert p.exists() and p.stat().st_size > 0
        manifest.append(dict(path=str(p.relative_to(root)), sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    plan = json.loads((prior / files[0]).read_text())
    assert plan['minimum_reference_observed_fraction'] == .8 and plan['minimum_observed_genes'] == 10
    cover = pd.read_csv(prior / files[1], sep='\t')
    assert len(cover) == 50 and cover.reference.is_unique
    expected = cover.source_positive_confirmed_fraction.ge(.8) & cover.source_positive_confirmed_genes.ge(10)
    assert np.array_equal(expected.to_numpy(), cover.fixed_partial_readout_ready.to_numpy())
    columns = np.load(prior / files[4])
    assert np.array_equal(columns, np.unique(columns))
    units = pd.read_csv(prior / files[3], sep='\t', keep_default_na=False)
    assert not units.duplicated(['patient_token', 'phase']).any()
    baseline = units[units.phase.eq('Pre')]
    eligible = baseline[(baseline.author_T_cells >= 20) & (baseline.author_myeloid_cells >= 20)]
    conflicts = eligible[[c for c in eligible if c.endswith('_conflict')]]
    # Explicit parsing survives CSV Boolean representation changes.
    conflict_mask = conflicts.apply(lambda c: c.astype(str).str.lower().eq('true')).any(axis=1)
    labelled = eligible[~conflict_mask & eligible.outcome.isin(['E', 'NE'])]
    contract_path = project / 'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    assert hashlib.sha256(contract_path.read_bytes()).hexdigest() == plan['source_contract_sha256']
    rec = next(r for r in json.loads(contract_path.read_text())['records'] if r['dataset'] == 'ICB_Bassez')
    with h5py.File(rec['path']) as handle:
        var = read_elem(handle['raw/var'])
        masks = [k for k in ['measured_in_source', 'feature_presence', 'source_gene_mask'] if k in handle.get('varm', {})]
    assert (columns >= 0).all() and (columns < len(var)).all()
    treatment = pd.crosstab(labelled['Treatment.or.Mode.of.Action'], labelled.outcome)
    report = dict(status='BASSEZ_FIXED_PUBLIC_FUNCTION_READINESS_AUDITED_NO_EFFECT_TEST',
        finished_UTC=datetime.now(timezone.utc).isoformat(), source_cells=int(len(np.load(rec['row_positions_ref']))),
        native_global_features=int(len(var)), source_positive_confirmed_features=int(len(columns)),
        per_source_measurement_mask_fields=masks, author_baseline_tokens=int(len(baseline)),
        both_lineages_atleast20=int(len(eligible)), conflict_free_E_NE_tokens=int(len(labelled)),
        treatment_by_expansion={str(c): {str(k): int(v) for k, v in treatment[c].items()} for c in treatment},
        all50_fixed_partial_readouts_ready=int(expected.sum()),
        prospective_IFN_readouts=cover[cover.reference.str.contains('INTERFERON')].to_dict('records'),
        accepted_old_fine_node_mapping=False, clinical_response_validation=False,
        old_MCP_weights_changed=False, new_expression_download=False,
        repair='v1 completed all count/coverage/capture writes then failed serializing NumPy int64; v2 reuses exact artifacts and casts JSON counts.',
        prior_artifact_manifest=manifest)
    out.mkdir()
    (out / 'private_terminal.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'prior_artifact_manifest'}))


if __name__ == '__main__':
    main()
