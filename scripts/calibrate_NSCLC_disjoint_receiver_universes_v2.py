"""All fixed LR-disjoint receiver slots under every accepted gene universe.

Only measurement fidelity on original frozen cases; no new association effects,
clinical outcome, external expression, feature selection or P is loaded.
"""
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import rankdata

from NSCLC_native_MFA_interpretation_inputs_v1 import load, sha, VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_fixed_programme_continuous_priors_v1 import read_rds, PRIORS
from preflight_NSCLC_native_interfaces_external_opportunities_v1 import accepted


def inputs():
    data = load(); root = data['root']
    native, native_proof = accepted(root,
        'NSCLC_continuous_prior_native_interfaces_20261016_v2',
        'private_independent_native_interface_audit.json',
        '1cf72d62c2df5f8a2afea70b46a5842c2da14737396b19faf2b9f8e4a5c39565')
    universe, universe_proof = accepted(root,
        'NSCLC_native_MFA_external_RNA_universe_calibration_20261016_v1',
        'private_independent_native_RNA_universe_audit.json',
        'f9822656c7be74fa7953b04c8807b5da718bc2084bf3a946e7943504ed31d5b3')
    with np.load(native / 'private_all_native_interface_measurements.npz', allow_pickle=False) as z:
        baseline = z['response']; native_ok = z['response_ok']
        assert np.array_equal(z['cases'], data['cases'].sampleID.to_numpy(str))
        assert np.array_equal(z['train'], data['train'])
    slots = pd.read_csv(native / 'private_all1600_original_slots.tsv.gz', sep='\t', keep_default_na=False)
    assert slots.slot_index.tolist() == list(range(1600))
    with np.load(universe / 'private_all_universe_native_scores.npz', allow_pickle=False) as z:
        masks = z['masks']; whole = z['scores']
    assert masks.shape == (48, 31831) and masks.dtype == bool and masks[0].all()
    assert whole.shape == (4, 48, 208, 10) and baseline.shape == (1600, 208)
    source = pd.read_csv(universe / 'private_all61_source_native_universe_map.tsv', sep='\t')
    gates = pd.read_csv(universe / 'private_all2440_native_RNA_universe_calibrations.tsv', sep='\t')
    assert len(source) == 61 and len(gates) == 2440 and gates.native_RNA_universe_gate.dtype == bool
    prior = root / 'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    path = prior / 'lr_network_human_21122021.rds'
    assert sha(path) == PRIORS[path.name]
    flat = read_rds(path); flat['from'] = flat['from'].astype(str); flat['to'] = flat['to'].astype(str)
    path = root / 'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv'
    assert sha(path) == 'df574533d866438e646423f29bd555cb283173dbceadfe9f8264e1622a86fd06'
    reference = pd.read_csv(path, sep='\t', keep_default_na=False)
    hashes = dict(native_audit=sha(native_proof), universe_audit=sha(universe_proof),
        slots=sha(native / 'private_all1600_original_slots.tsv.gz'), model=sha(data['model_path']))
    return data, native, universe, baseline, native_ok, slots, masks, whole, source, gates, flat, reference, hashes


def metrics(prediction, baseline, train):
    # Population SD matches the original native LR-removal calibrate() scale.
    sd = baseline[train].std(axis=0, ddof=0)
    valid_base = np.isfinite(baseline).all(axis=0) & np.isfinite(sd) & (sd > 1e-12)
    result = {}; passed = valid_base.copy()
    for label, take in [('discovery', train), ('internal_holdout', ~train)]:
        a, b = baseline[take], prediction[take]
        aa, bb = a-a.mean(axis=0), b-b.mean(axis=0)
        ar, br = rankdata(a, axis=0), rankdata(b, axis=0)
        ar -= ar.mean(axis=0); br -= br.mean(axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            pearson = (aa*bb).sum(axis=0)/(np.linalg.norm(aa, axis=0)*np.linalg.norm(bb, axis=0))
            spearman = (ar*br).sum(axis=0)/(np.linalg.norm(ar, axis=0)*np.linalg.norm(br, axis=0))
            rmse = np.sqrt(np.mean((b-a)**2, axis=0))/sd
            shift = (b-a).mean(axis=0)/sd
        prediction_variable = b.std(axis=0, ddof=0) > 1e-12
        pearson[~valid_base | ~prediction_variable] = np.nan
        spearman[~valid_base | ~prediction_variable] = np.nan
        rmse[~valid_base] = np.nan; shift[~valid_base] = np.nan
        for name, values in [('Pearson', pearson), ('Spearman', spearman), ('RMSE_SD', rmse), ('mean_shift_SD', shift)]:
            result[label+'_'+name] = values
        passed &= np.isfinite(pearson) & np.isfinite(rmse) & (pearson >= .95) & (rmse <= .25)
    return result, passed, valid_base, sd


def source_table(table, source, gates):
    columns = ['slot_index', 'view', 'factor', 'direction', 'ligand', 'mask_index',
               'native_response_ok', 'preserves_native_disjoint_response', 'status']
    result = table[columns].merge(source[['arm', 'dataset', 'mask_index']], on='mask_index', how='inner', validate='many_to_many')
    old = gates[['arm', 'dataset', 'view', 'factor', 'native_RNA_universe_gate']].rename(
        columns={'native_RNA_universe_gate': 'original_full_source_universe_gate'})
    result = result.merge(old, on=['arm', 'dataset', 'view', 'factor'], how='left', validate='many_to_one')
    assert result.original_full_source_universe_gate.notna().all()
    result['source_disjoint_receiver_universe_gate'] = (
        result.native_response_ok & result.preserves_native_disjoint_response & result.original_full_source_universe_gate)
    result['status'] = [
        status if status != 'CALIBRATION_PASSED' else
        ('CALIBRATION_PASSED' if old else 'ORIGINAL_FULL_SOURCE_UNIVERSE_GATE_FAILED')
        for status, old in zip(result.status, result.original_full_source_universe_gate)]
    assert len(result) == 97600 and not result.duplicated(['slot_index', 'arm', 'dataset']).any()
    return result.sort_values(['slot_index', 'arm', 'dataset']).reset_index(drop=True)


def main():
    state = inputs()
    data, native, universe, baseline, native_ok, slots, masks, whole, source, gates, flat, reference, hashes = state
    root = data['root']; budget = resource(root, 160*1024**2)
    available = int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')))*1024
    assert available > 2*1024**3
    out = root / 'NSCLC_disjoint_receiver_universe_calibration_20261018_v2'
    assert not out.exists(); out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), inputs=hashes,
        budget=budget, slots=1600, masks=48, source_arms=61, cases=208, discovery=145, internal_holdout=63,
        measurement='Original frozen gene-universe TMM, quality coefficients/design and RMS. Whole-universe own score minus measured literal LR-excluded contributions. No retained-feature renormalization or missing-expression zero claim.',
        gate='Both original partitions Pearson>=.95 and RMSE<=.25 native LR-disjoint discovery populationSD(ddof0); native response_ok and original full source universe gate also required. No old-gate rescue.',
        scientific_limits='Only native measurement-universe perturbation, not external pool identity, ligand assay, same-cell receptor, association or clinical validation.',
        no_new_association_results_P_outcomes_external_expression_or_downloads=True)
    (out / 'plan_before_disjoint_receiver_universe_calibration.json').write_text(json.dumps(plan, indent=2)+'\n')
    all_scores = np.full((1600, 48, 208), np.nan)
    rows = []; baseline_error = 0.; inp = root / 'NSCLC_native_MFA_inputs_20261016_v1'
    original_inputs = json.loads((inp / 'private_terminal.json').read_text())
    cases = data['cases'].sampleID.astype(str).tolist(); train = data['train']
    for vi, view in enumerate(VIEWS):
        p = data['prepared'][vi]; selected = p['selected_native_columns']; genes = data['genes'][selected]
        w = data['model']['weights_'+str(vi)]; tau = data['model']['precision_'+str(vi)]
        projection = np.linalg.solve(np.eye(10)+w.T@(tau[:, None]*w), w.T*tau)
        fixed_quality = p['design']@p['regression_coefficients']; rms = float(p['training_view_RMS'])
        path = inp / ('private_'+view+'_pseudobulk.tsv.gz')
        declared = original_inputs['files'][path.name]
        assert isinstance(declared, str) and sha(path) == declared
        counts = pd.read_csv(path, sep='\t', index_col=0)
        assert counts.columns.tolist() == cases and counts.index.astype(str).tolist() == data['genes'].tolist()
        raw = counts.iloc[selected].to_numpy().T
        norm = pd.read_csv(universe / ('private_'+view+'_all_universe_TMM.tsv.gz'), sep='\t', dtype={'sampleID': str})
        local = slots[slots.view.eq(view) & slots.ligand.ne('')].copy()
        ids = local.slot_index.to_numpy(int); factors = local.factor.str.rsplit('F', n=1).str[-1].astype(int).to_numpy()-1
        assert len(ids) > 0
        omitted = {}
        for ligand in local.ligand.unique():
            names = {ligand, *flat.loc[flat['from'].eq(ligand), 'to'], *reference.loc[reference.ligand.eq(ligand), 'receptor']}
            omit = set(names)
            for name in names: omit.update(str(name).split('_'))
            omitted[ligand] = np.flatnonzero(np.isin(genes, list(omit)))
        rr, cc, vv = [], [], []
        for j, row in enumerate(local.itertuples()):
            columns = omitted[row.ligand]
            rr.extend(columns); cc.extend([j]*len(columns)); vv.extend(projection[factors[j], columns])
        removed = csr_matrix((vv, (rr, cc)), shape=(len(selected), len(ids)))
        local_baseline = baseline[ids].T
        for mi, mask in enumerate(masks):
            take = mask[selected]
            normalization = norm[norm.mask_id.eq('mask_%03d' % mi)].set_index('sampleID').loc[cases]
            scale = normalization.effective_library_size.to_numpy(float)
            assert np.isfinite(scale).all() and (scale > 0).all()
            value = (np.log2(1+raw[:, take].astype(float)*1e6/scale[:, None])-fixed_quality[:, take])/rms
            prediction = whole[vi, mi][:, factors] - (removed[take].T@value.T).T
            assert np.isfinite(prediction).all()
            all_scores[ids, mi, :] = prediction.T
            if mi == 0:
                delta = float(np.max(abs(prediction-local_baseline))); baseline_error = max(baseline_error, delta)
                assert delta < 1e-8
            values, passed, variable, sd = metrics(prediction, local_baseline, train)
            for j, row in enumerate(local.itertuples()):
                reason = 'CALIBRATION_PASSED'
                if not native_ok[row.slot_index]: reason = 'ORIGINAL_NATIVE_RESPONSE_GATE_FAILED'
                elif not variable[j]: reason = 'CONSTANT_OR_INVALID_NATIVE_DISJOINT_RESPONSE'
                elif not passed[j]: reason = 'DISJOINT_RESPONSE_UNIVERSE_FIDELITY_FAILED'
                rows.append(dict(slot_index=int(row.slot_index), view=view, factor=row.factor,
                    direction=int(row.direction), ligand=row.ligand, mask_index=mi,
                    native_response_ok=bool(native_ok[row.slot_index]),
                    native_disjoint_discovery_populationSD=float(sd[j]),
                    retained_selected_features=int(take.sum()-np.count_nonzero(take[omitted[row.ligand]])),
                    preserves_native_disjoint_response=bool(passed[j]), status=reason,
                    **{key: float(val[j]) for key, val in values.items()}))
            for row in slots[slots.view.eq(view) & slots.ligand.eq('')].itertuples():
                rows.append(dict(slot_index=int(row.slot_index), view=view, factor=row.factor,
                    direction=int(row.direction), ligand='', mask_index=mi, native_response_ok=False,
                    native_disjoint_discovery_populationSD=np.nan, retained_selected_features=0,
                    preserves_native_disjoint_response=False, status='EMPTY_ORIGINAL_FIXED_SLOT',
                    **{key: np.nan for key in values}))
        (out / 'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), completed_view=view, completed_views=vi+1, total_views=4))+'\n')
    table = pd.DataFrame(rows).sort_values(['slot_index', 'mask_index']).reset_index(drop=True)
    assert len(table) == 76800 and not table.duplicated(['slot_index', 'mask_index']).any()
    assert np.isnan(all_scores[slots.ligand.eq('').to_numpy()]).all()
    resource(root, 160*1024**2)
    table.to_csv(out / 'private_all76800_disjoint_universe_calibrations.tsv.gz', sep='\t', index=False, compression='gzip', na_rep='NA')
    source_table(table, source, gates).to_csv(out / 'private_all97600_source_disjoint_receiver_calibrations.tsv.gz', sep='\t', index=False, compression='gzip', na_rep='NA')
    np.savez_compressed(out / 'private_all1600x48x208_disjoint_receiver_scores.npz', scores=all_scores)
    report = dict(status='COMPLETE_ALL_DISJOINT_RECEIVER_UNIVERSE_CALIBRATIONS_INDEPENDENT_AUDIT_PENDING',
        finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        full_universe_maximum_native_response_error=baseline_error,
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        new_associations_P_clinical_or_external_biological_validation=False)
    assert sum(x['bytes'] for x in report['files'].values()) < budget['output_cap_bytes']
    (out / 'private_terminal.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status=report['status'], finished_UTC=report['finished_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
