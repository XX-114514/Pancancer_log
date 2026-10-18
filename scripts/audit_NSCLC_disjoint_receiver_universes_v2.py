"""Independent QR/direct-kept-feature full replay of every calibration score."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.stats import spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import sha, VIEWS
from calibrate_NSCLC_disjoint_receiver_universes_v2 import inputs
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from review_NSCLC_native_interface_associations_v1 import bound_wait


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'NSCLC_disjoint_receiver_universe_calibration_20261018_v2'
    terminal = out / 'private_terminal.json'
    bound_wait(root, 'calibrate_NSCLC_disjoint_receiver_universes_v2', terminal)
    resource(root, 1024**2)
    manifest = json.loads(terminal.read_text())
    for name, item in manifest['files'].items(): assert sha(out / name) == item['sha256']
    state = inputs()
    data, native, universe, baseline, native_ok, slots, masks, whole, source, gates, flat, reference, hashes = state
    assert json.loads((out / 'plan_before_disjoint_receiver_universe_calibration.json').read_text())['inputs'] == hashes
    table = pd.read_csv(out / 'private_all76800_disjoint_universe_calibrations.tsv.gz', sep='\t', keep_default_na=False, na_values=['NA']).set_index(['slot_index', 'mask_index'])
    assert len(table) == 76800 and table.index.is_unique
    with np.load(out / 'private_all1600x48x208_disjoint_receiver_scores.npz', allow_pickle=False) as z: saved = z['scores']
    assert saved.shape == (1600, 48, 208) and saved.dtype == np.float64
    train = data['train']; cases = data['cases'].sampleID.astype(str).tolist()
    flat_sets = flat.groupby('from')['to'].agg(set).to_dict()
    ref_sets = reference.groupby('ligand').receptor.agg(set).to_dict()
    max_error = 0.; max_metric = 0.; checked = 0
    checkpoint = root / 'HNSCC_receiver_followup_job_control_20261016_v2' / Path(__file__).stem / 'private_audit_checkpoint.json'
    for vi, view in enumerate(VIEWS):
        p = data['prepared'][vi]; columns = p['selected_native_columns']; genes = data['genes'][columns]
        w = data['model']['weights_'+str(vi)]; sw = np.sqrt(data['model']['precision_'+str(vi)])
        q, r = np.linalg.qr(np.vstack([sw[:, None]*w, np.eye(10)]), mode='reduced')
        projection = solve_triangular(r, q[:len(w)].T*sw)
        fixed = p['design']@p['regression_coefficients']; rms = float(p['training_view_RMS'])
        path = root / 'NSCLC_native_MFA_inputs_20261016_v1' / ('private_'+view+'_pseudobulk.tsv.gz')
        declared = json.loads((path.parent / 'private_terminal.json').read_text())['files'][path.name]
        assert isinstance(declared, str) and sha(path) == declared
        counts = pd.read_csv(path, sep='\t', index_col=0)
        assert counts.columns.tolist() == cases and counts.index.astype(str).tolist() == data['genes'].tolist()
        raw = counts.iloc[columns].to_numpy().T
        norm = pd.read_csv(universe / ('private_'+view+'_all_universe_TMM.tsv.gz'), sep='\t', dtype={'sampleID': str})
        local = slots[slots.view.eq(view) & slots.ligand.ne('')]
        ids = local.slot_index.to_numpy(int)
        factors = np.array([int(v.rsplit('F', 1)[1])-1 for v in local.factor])
        kept = projection[factors].T.copy(); keep_masks = []
        for j, row in enumerate(local.itertuples()):
            omitted = {row.ligand} | flat_sets.get(row.ligand, set()) | ref_sets.get(row.ligand, set())
            omitted = omitted | {g for text in omitted for g in str(text).split('_')}
            take = np.array([name not in omitted for name in genes])
            kept[~take, j] = 0.; keep_masks.append(take)
        reference_values = baseline[ids].T; sd = np.std(reference_values[train], axis=0, ddof=0)
        valid_reference = np.isfinite(reference_values).all(axis=0) & (sd > 1e-12)
        for mi, mask in enumerate(masks):
            take = mask[columns]
            scale = norm[norm.mask_id.eq('mask_%03d' % mi)].set_index('sampleID').loc[cases].effective_library_size.to_numpy(float)
            assert (scale > 0).all() and np.isfinite(scale).all()
            value = (np.log2(1+raw[:, take].astype(float)*1e6/scale[:, None])-fixed[:, take])/rms
            prediction = value@kept[take]
            actual = saved[ids, mi, :].T
            delta = float(np.max(abs(prediction-actual))); max_error = max(max_error, delta)
            assert np.allclose(prediction, actual, atol=1e-8, rtol=1e-9), delta
            if mi == 0: assert np.allclose(prediction, reference_values, atol=1e-8, rtol=1e-9)
            for j, row in enumerate(local.itertuples()):
                observed = table.loc[(row.slot_index, mi)]
                assert observed['view'] == view and observed.factor == row.factor and observed.direction == row.direction and observed.ligand == row.ligand
                expected = {}; passed = bool(valid_reference[j])
                for label, partition in [('discovery', train), ('internal_holdout', ~train)]:
                    a, b = reference_values[partition, j], prediction[partition, j]
                    if valid_reference[j]:
                        a0, b0 = a-a.mean(), b-b.mean()
                        variable = np.std(b, ddof=0) > 1e-12
                        rho = float(a0@b0/np.sqrt((a0@a0)*(b0@b0))) if variable and a0@a0 > 0 else np.nan
                        rank = float(spearmanr(a, b).statistic) if variable and np.ptp(a) > 0 else np.nan
                        d = b-a; rmse = float(np.sqrt(d@d/len(d))/sd[j]); shift = float(d.mean()/sd[j])
                    else: rho = rank = rmse = shift = np.nan
                    expected.update({label+'_Pearson': rho, label+'_Spearman': rank, label+'_RMSE_SD': rmse, label+'_mean_shift_SD': shift})
                    passed = bool(passed and np.isfinite(rho) and np.isfinite(rmse) and rho >= .95 and rmse <= .25)
                for name, v in expected.items():
                    ov = float(observed[name])
                    assert np.allclose(v, ov, atol=1e-8, rtol=1e-8, equal_nan=True), (name, v, ov)
                    if np.isfinite(v): max_metric = max(max_metric, abs(v-ov))
                assert np.allclose(observed.native_disjoint_discovery_populationSD, sd[j], atol=1e-12, rtol=1e-12)
                assert observed.retained_selected_features == int(np.count_nonzero(take & keep_masks[j]))
                assert observed.native_response_ok == native_ok[row.slot_index]
                assert observed.preserves_native_disjoint_response == passed
                reason = 'CALIBRATION_PASSED'
                if not native_ok[row.slot_index]: reason = 'ORIGINAL_NATIVE_RESPONSE_GATE_FAILED'
                elif not valid_reference[j]: reason = 'CONSTANT_OR_INVALID_NATIVE_DISJOINT_RESPONSE'
                elif not passed: reason = 'DISJOINT_RESPONSE_UNIVERSE_FIDELITY_FAILED'
                assert observed.status == reason
                checked += 1
            checkpoint.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), view=view, completed_masks=mi+1, total_masks=48))+'\n')
        for row in slots[slots.view.eq(view) & slots.ligand.eq('')].itertuples():
            assert np.isnan(saved[row.slot_index]).all()
            for mi in range(48):
                observed = table.loc[(row.slot_index, mi)]
                assert observed.status == 'EMPTY_ORIGINAL_FIXED_SLOT' and not observed.native_response_ok and not observed.preserves_native_disjoint_response
                assert observed.retained_selected_features == 0
                assert observed['view'] == row.view and observed.factor == row.factor and observed.direction == row.direction and observed.ligand == ''
                assert np.isnan(observed.native_disjoint_discovery_populationSD)
                assert all(np.isnan(observed[key]) for key in expected)
                checked += 1
    assert checked == 76800
    observed = pd.read_csv(out / 'private_all97600_source_disjoint_receiver_calibrations.tsv.gz', sep='\t', keep_default_na=False, na_values=['NA'])
    assert len(observed) == 97600 and not observed.duplicated(['slot_index', 'arm', 'dataset']).any()
    by_source = {(r.arm, r.dataset): int(r.mask_index) for r in source.itertuples()}
    old = {(r.arm, r.dataset, r.view, r.factor): bool(r.native_RNA_universe_gate) for r in gates.itertuples()}
    for row in observed.itertuples():
        mi = by_source[(row.arm, row.dataset)]; assert row.mask_index == mi
        base = table.loc[(row.slot_index, mi)]
        for name in ['view', 'factor', 'direction', 'ligand', 'native_response_ok', 'preserves_native_disjoint_response']:
            assert getattr(row, name) == base[name]
        full_gate = old[(row.arm, row.dataset, row.view, row.factor)]
        assert row.original_full_source_universe_gate == full_gate
        assert row.source_disjoint_receiver_universe_gate == (base.native_response_ok and base.preserves_native_disjoint_response and full_gate)
        reason = base.status
        if reason == 'CALIBRATION_PASSED' and not full_gate: reason = 'ORIGINAL_FULL_SOURCE_UNIVERSE_GATE_FAILED'
        assert row.status == reason
    report = dict(status='PASS_ALL1600x48x208_SCORES76800_CALIBRATIONS97600_SOURCE_GATES_QR_DIRECT_REPLAY',
        finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        producer_terminal_sha256=sha(terminal), all_score_positions_checked=int(saved.size),
        calibration_rows_checked=checked, source_rows_checked=len(observed),
        maximum_score_error=max_error, maximum_metric_error=max_metric,
        new_associations_clinical_or_external_biological_validation=False)
    proof = out / 'private_independent_disjoint_universe_audit.json'
    assert not proof.exists(); proof.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status=report['status'], finished_UTC=report['finished_UTC'], audit_sha256=sha(proof))))


if __name__ == '__main__':
    main()
