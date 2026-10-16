"""Independent all-slot same-library, patient and numeric aggregation audit."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_external_lineage_inputs_v2 import load


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    resource(root, 1048576)
    src = root / 'NSCLC_external_programme_pair_inputs_20261016_v1'
    tpath = src / 'private_terminal.json'
    terminal = json.loads(tpath.read_text())
    assert terminal['scientific_effects_or_P_computed'] is False
    for name, item in terminal['files'].items():
        assert sha(src/name) == item['sha256']
    plan = json.loads((src/'plan_before_any_external_pair_association.json').read_text())
    assert plan['min_units'] == 20 and plan['min_residual_df'] == 8
    assert plan['permutations'] == 19999 and plan['bootstrap'] == 2000
    assay = root/'NSCLC_external_fixed_gene_programme_assay_20261016_v1'
    units = pd.read_csv(assay/'private_all6760_programme_units_metadata_quality_TMM.tsv.gz', sep='\t', keep_default_na=False)
    with np.load(assay/'private_all6760_fixed_programme_indices.npz', allow_pickle=False) as z:
        scores = {k: z[k] for k in z.files}
    slots = pd.read_csv(src/'private_all_frozen_pair_slots_and_coverage.tsv', sep='\t', keep_default_na=False)
    cases = pd.read_csv(src/'private_all_pair_analysis_units.tsv.gz', sep='\t', keep_default_na=False)
    members = pd.read_csv(src/'private_same_library_pair_membership.tsv.gz', sep='\t')
    excluded = pd.read_csv(src/'private_excluded_mixed_patient_contexts.tsv.gz', sep='\t')
    lib = load()['lib']
    denominator = lib.annotated_scope_cells.to_numpy(float)
    assert len(slots) == 7800 and slots.slot.is_unique
    assert slots.groupby('context_id').size().eq(60).all()
    assert cases.row_id.tolist() == list(range(len(cases)))
    assert not cases.duplicated(['slot', 'resampling_block']).any()
    assert not members.duplicated(['slot', 'library_index']).any()
    membership = {int(k): g for k, g in members.groupby('row_id')}
    view_order = units.iloc[:4]['view'].tolist()
    maximum = 0.
    for row in cases.itertuples():
        g = membership[row.row_id]
        a, b = g.unit_A.to_numpy(int), g.unit_B.to_numpy(int)
        li = g.library_index.to_numpy(int)
        assert np.array_equal(a//4, li) and np.array_equal(b//4, li)
        assert len(g) == row.paired_libraries and g.slot.eq(row.slot).all()
        assert units.iloc[a].resampling_block.eq(row.resampling_block).all()
        assert units.iloc[b].resampling_block.eq(row.resampling_block).all()
        assert units.iloc[a].reviewed_patient.eq(row.reviewed_patient).all()
        assert units.iloc[b].reviewed_patient.eq(row.reviewed_patient).all()
        fields = ['tissue', 'treatment_arm', 'capture_scope', 'ex_vivo']
        keys = [json.dumps([str(units.iloc[i][c]) for c in fields], ensure_ascii=False, separators=(',', ':')) for i in a]
        assert set(keys) == {row.block_key}
        # Factor is fixed by the slot, never inferred from score similarity.
        factor = row.slot.split('|')[1]
        k = list(scores['factors'].astype(str)).index(factor)
        for ids in [a, b]:
            assert units.iloc[ids].RNA_measurement_available.all()
            assert scores['source_gene_gate'][ids, k].all()
            assert units.iloc[ids].n_cells.ge(20).all()
        expected = np.r_[scores['raw'][a, k].mean(), scores['raw'][b, k].mean(),
                         scores['quality'][a].mean(axis=0), scores['quality'][b].mean(axis=0),
                         np.mean(units.n_cells.to_numpy()[a]/denominator[li]),
                         np.mean(units.n_cells.to_numpy()[b]/denominator[li])]
        observed = np.array([row.raw_A, row.raw_B] + [getattr(row, f'Q_A_{j}') for j in range(5)] +
                            [getattr(row, f'Q_B_{j}') for j in range(5)] + [row.fraction_A, row.fraction_B], float)
        assert np.allclose(expected, observed, rtol=1e-12, atol=1e-12)
        maximum = max(maximum, float(np.max(np.abs(expected-observed))))
        for ids, field in [(a, 'source_QC_A'), (b, 'source_QC_B')]:
            x = scores['source_QC'][ids, k]
            expected_value = x.mean() if np.isfinite(x).all() else np.nan
            observed_value = pd.to_numeric(getattr(row, field), errors='coerce')
            assert np.isclose(expected_value, observed_value, equal_nan=True, atol=1e-12, rtol=1e-12)
    case_groups = {k: g for k, g in cases.groupby('slot')}
    ex_groups = {k: g for k, g in excluded.groupby('slot')}
    # Reconstruct the eligible library denominator independently of membership.
    def phase(text):
        s = str(text).strip().casefold().replace('-', '_').replace(' ', '_')
        if s in {'', 'unknown', 'not_reviewed', 'nan', 'none'}: return 'unknown_recorded_time'
        if s in {'pre', 'pre_treatment', 'pretreatment', 'baseline'}: return 'pre_treatment'
        if s in {'post', 'post_treatment', 'posttreatment'}: return 'post_treatment'
        if s in {'on_trt', 'on_treatment', 'ontreatment'}: return 'on_treatment'
        return s
    base = units.iloc[::4].copy().reset_index(drop=True)
    base['phase_key'] = base.timepoint.map(phase)
    contexts = {tuple(str(v) for v in key): g.index.to_numpy(int) for key, g in base.groupby(plan['context'], sort=True, dropna=False)}
    for row in slots.itertuples():
        li = contexts[tuple(str(getattr(row, c)) for c in plan['context'])]
        a, b = li*4+view_order.index(row.view_A), li*4+view_order.index(row.view_B)
        k = list(scores['factors'].astype(str)).index(row.factor)
        rna = units.iloc[a].RNA_measurement_available.to_numpy() & units.iloc[b].RNA_measurement_available.to_numpy()
        good = rna & scores['source_gene_gate'][a, k] & scores['source_gene_gate'][b, k]
        assert row.context_libraries == len(li) and row.both_RNA_measurable_libraries == rna.sum()
        assert row.both_gene_gate_libraries == good.sum()
        paired = base.iloc[li[good]]
        expected_groups = set(paired.resampling_block.astype(str))
        actual_groups = set(case_groups[row.slot].resampling_block.astype(str)) if row.slot in case_groups else set()
        rejected_groups = set(ex_groups[row.slot].resampling_block.astype(str)) if row.slot in ex_groups else set()
        assert not (actual_groups & rejected_groups)
        assert expected_groups == actual_groups | rejected_groups
        for block, g in paired.groupby('resampling_block', sort=True):
            keys = {tuple(str(r[c]) for c in plan['nuisance_exchangeability_block']) for _, r in g.iterrows()}
            assert (block in rejected_groups) == (len(keys)>1)
        observed_cases = case_groups.get(row.slot)
        assert row.all_units == (len(observed_cases) if observed_cases is not None else 0)
        assert row.reviewed_units == (int(observed_cases.reviewed_patient.sum()) if observed_cases is not None else 0)
        assert row.mixed_block_units_excluded == len(rejected_groups)
        assert row.in_model_source == (row.dataset == 'GSE243013')
    proof = src/'private_independent_pair_input_audit.json'
    assert not proof.exists()
    report = dict(status='PASS_ALL7800_SLOTS_SAME_LIBRARY_MEMBERSHIP_PATIENT_AGGREGATES_AND_MISSING_DENOMINATORS',
                  finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                  producer_terminal_sha256=sha(tpath), all_slots_checked=len(slots),
                  all_analysis_rows_checked=len(cases), maximum_aggregation_error=maximum,
                  scientific_effects_or_P_computed=False)
    proof.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status=report['status'], finished_UTC=report['finished_UTC'], audit_sha256=sha(proof))))


if __name__ == '__main__':
    main()
