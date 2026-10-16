"""Complete prespecified all-source coordination family, no outcome selection."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist, spearmanr
from statsmodels.stats.multitest import multipletests

import NSCLC_external_programme_coordination_math_v2 as kernel
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

QCOLS = [f'Q_{side}_{j}' for side in ['A', 'B'] for j in range(5)]


def read_inputs(root):
    src = root/'NSCLC_external_programme_pair_inputs_20261016_v1'
    proof = src/'private_independent_pair_input_audit.json'
    audit = json.loads(proof.read_text())
    assert audit['status'].startswith('PASS_ALL7800_SLOTS_SAME_LIBRARY_')
    assert audit['producer_terminal_sha256'] == sha(src/'private_terminal.json')
    terminal = json.loads((src/'private_terminal.json').read_text())
    for name, item in terminal['files'].items():
        assert sha(src/name) == item['sha256']
    contract = json.loads((src/'plan_before_any_external_pair_association.json').read_text())
    assert (contract['min_units'], contract['min_residual_df'], contract['permutations'], contract['bootstrap']) == (20, 8, 19999, 2000)
    slots = pd.read_csv(src/'private_all_frozen_pair_slots_and_coverage.tsv', sep='\t', keep_default_na=False)
    cases = pd.read_csv(src/'private_all_pair_analysis_units.tsv.gz', sep='\t', keep_default_na=False)
    for c in QCOLS + ['raw_A', 'raw_B', 'fraction_A', 'fraction_B']:
        cases[c] = pd.to_numeric(cases[c], errors='raise')
    for c in ['source_QC_A', 'source_QC_B']:
        cases[c] = pd.to_numeric(cases[c], errors='coerce')
    assert pd.api.types.is_bool_dtype(cases.reviewed_patient)
    return src, proof, contract, slots, cases


def qualify(slots, cases):
    records, prepared, frames = [], [], []
    for row in slots.itertuples():
        g0 = cases.iloc[row.row_start:row.row_stop]
        assert len(g0) == row.all_units and (not len(g0) or g0.slot.eq(row.slot).all())
        for family in ['all_units', 'reviewed_patients']:
            g = g0 if family == 'all_units' else g0[g0.reviewed_patient]
            record = row.slot+'|'+family
            item = dict(record=record, family=family, slot=row.slot, candidate_units=len(g),
                        candidate_reviewed_units=int(g.reviewed_patient.sum()),
                        status='ELIGIBLE_FOR_FIXED_PRIMARY_TEST', active_units=0, singleton_units=0,
                        active_reviewed_units=0, quality_rank=0, residual_df=0, test_index=-1)
            if row.in_model_source:
                item['status'] = 'ORIGINAL_MODEL_SOURCE_NOT_EXTERNAL_FAMILY'
            else:
                try:
                    p = kernel.prepare(g.raw_A.to_numpy(), g.raw_B.to_numpy(), g[QCOLS].to_numpy(), g.block_key.to_numpy(str))
                    active = g.iloc[np.flatnonzero(p['keep'])]
                    item.update(active_units=len(active), singleton_units=int((~p['keep']).sum()),
                                active_reviewed_units=int(active.reviewed_patient.sum()), quality_rank=p['design_rank'],
                                residual_df=p['residual_df'], test_index=len(prepared))
                    prepared.append(p)
                    frames.append(active)
                except ValueError as e:
                    item['status'] = str(e)
                    counts = g.block_key.value_counts()
                    item['singleton_units'] = int(g.block_key.map(counts).eq(1).sum())
            records.append(item)
    return pd.DataFrame(records), prepared, frames


def secondary(a, b, q, blocks):
    try:
        return kernel.effect(kernel.prepare(a, b, q, blocks)), 'ESTIMABLE'
    except ValueError as e:
        return np.nan, str(e)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 64*1024**2)
    src, proof, contract, slots, cases = read_inputs(root)
    kernel_proof = root/'NSCLC_external_coordination_kernel_check_20261016_v2/private_independent_kernel_check.json'
    checked = json.loads(kernel_proof.read_text())
    assert checked['status'].startswith('PASS_SYNTHETIC_QR_VS_OLS_')
    assert checked['math_sha256'] == sha(kernel.__file__)
    out = root/'NSCLC_external_programme_coordination_20261016_v2'
    assert not out.exists()
    out.mkdir(mode=0o700)
    scope, prepared, frames = qualify(slots, cases)
    # Write eligibility and complete multiplicity denominators before any rho/P.
    scope.to_csv(out/'private_all15600_primary_test_eligibility.tsv', sep='\t', index=False)
    sizes = {f: int(((scope.family == f) & (scope.test_index >= 0)).sum()) for f in scope.family.unique()}
    assert sum(sizes.values()) == len(prepared)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__),
                math_sha256=sha(kernel.__file__), input_audit_sha256=sha(proof),
                kernel_check_sha256=sha(kernel_proof),
                original_frozen_contract_sha256=sha(src/'plan_before_any_external_pair_association.json'),
                eligibility_sha256=sha(out/'private_all15600_primary_test_eligibility.tsv'),
                family_sizes_before_any_effect=sizes, contract=contract, budget=budget,
                clinical_outcomes_loaded=False, direct_communication_test=False)
    (out/'plan_before_all_external_coordination_effects.json').write_text(json.dumps(plan, indent=2)+'\n')
    boots = np.full((len(prepared), 2000), np.nan)
    loo_arrays = []
    rows = []
    metadata = slots.set_index('slot')
    for item in scope[scope.test_index.ge(0)].itertuples():
        i = item.test_index
        p, frame = prepared[i], frames[i]
        rho = kernel.effect(p)
        count = kernel.permutation_count(p, 19999, kernel.seed(item.record, 'Freedman_Lane'))
        pv = (1+count)/20000
        upper = 1. if count == 19999 else float(beta_dist.isf(.01/sizes[item.family], count+1, 19999-count))
        upper = min(1., (1+19999*upper)/20000)
        boots[i] = kernel.bootstrap(p, 2000, kernel.seed(item.record, 'bootstrap'))
        good = np.isfinite(boots[i])
        lo, hi = np.quantile(boots[i, good], [.025, .975]) if good.sum() >= 1900 else [np.nan, np.nan]
        loo = kernel.leave_one(p)
        loo_arrays.append(loo)
        empty = np.empty((len(frame), 0))
        block_effect, block_status = secondary(p['a'], p['b'], empty, p['blocks'])
        qc_effect, qc_status = secondary(frame.source_QC_A.to_numpy(), frame.source_QC_B.to_numpy(), empty, p['blocks'])
        abund = np.column_stack([p['quality'], frame[['fraction_A', 'fraction_B']].to_numpy()])
        abundance_effect, abundance_status = secondary(p['a'], p['b'], abund, p['blocks'])
        details = metadata.loc[item.slot].to_dict()
        details['context_reviewed_units'] = details.pop('reviewed_units')
        rows.append(dict(record=item.record, family=item.family, slot=item.slot, test_index=i,
                         **details, active_units=len(frame), reviewed_units=int(frame.reviewed_patient.sum()),
                         quality_rank=p['design_rank'], residual_df=p['residual_df'],
                         status='ESTIMATED_EXPLORATORY_COHORT_PROGRAMME_COORDINATION',
                         partial_rank_rho=rho, P=pv, MC_upper_P=upper, permutation_exceedances=count,
                         CI_low=float(lo), CI_high=float(hi), valid_bootstraps=int(good.sum()),
                         raw_Spearman=float(spearmanr(p['a'], p['b']).statistic),
                         block_only_rho=block_effect, block_only_status=block_status,
                         source_QC_block_rho=qc_effect, source_QC_status=qc_status,
                         abundance_adjusted_rho=abundance_effect, abundance_status=abundance_status,
                         LOO_valid=int(np.isfinite(loo).sum()), LOO_all_valid=bool(np.isfinite(loo).all()),
                         LOO_same_sign=bool(np.isfinite(loo).all() and (np.sign(loo) == np.sign(rho)).all()),
                         LOO_min=float(np.nanmin(loo)) if np.isfinite(loo).any() else np.nan,
                         LOO_max=float(np.nanmax(loo)) if np.isfinite(loo).any() else np.nan))
        if (i+1) % 10 == 0 or i+1 == len(prepared):
            (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), completed_tests=i+1, total_tests=len(prepared), stage='fixed_full_family_resampling'))+'\n')
    result = pd.DataFrame(rows)
    for family, g in result.groupby('family', sort=True):
        assert len(g) == sizes[family]
        for target, method, source in [('BH_q', 'fdr_bh', 'P'), ('BY_q', 'fdr_by', 'P'), ('MC_upper_BH_q', 'fdr_bh', 'MC_upper_P')]:
            result.loc[g.index, target] = multipletests(g[source].to_numpy(), method=method)[1]
    result.to_csv(out/'private_all_estimated_programme_coordination.tsv', sep='\t', index=False, na_rep='NA')
    offsets = np.r_[0, np.cumsum([len(a) for a in loo_arrays])]
    np.savez_compressed(out/'private_all_bootstrap_and_LOO.npz', bootstrap=boots,
                        LOO=np.concatenate(loo_arrays) if loo_arrays else np.empty(0), LOO_offsets=offsets,
                        records=result.record.to_numpy(str))
    report = dict(status='COMPLETE_ALL_EXTERNAL_COHORT_PROGRAMME_COORDINATION_FAMILIES_AUDIT_PENDING',
                  finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                  files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
                  new_direct_LR_clinical_or_independent_patient_validation=False)
    assert sum(x['bytes'] for x in report['files'].values()) < budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status=report['status'], finished_UTC=report['finished_UTC'], terminal_sha256=sha(out/'private_terminal.json'))))


if __name__ == '__main__':
    main()
