"""Complete frozen external family inference, with no outcome/endpoint search."""
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from external_core_rank_model_v2 import point, wild_count, monte_carlo_upper, bootstrap_rhos, leave_one_out


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def wait_proof(root, relative, job):
    path = root / relative
    receipt = json.loads((root / 'HNSCC_receiver_followup_job_control_20261016_v2' / job / 'launch_receipt.json').read_text())
    deadline = time.monotonic() + 3600
    while not path.exists():
        p = Path('/proc') / str(receipt['pid']) / 'stat'
        assert p.exists()
        s = p.read_text().rsplit(')', 1)[1].split()
        assert s[19] == str(receipt['proc_start_ticks']) and s[0] != 'Z' and time.monotonic() < deadline
        time.sleep(5)
    assert json.loads(path.read_text())['status'].startswith('PASS_')
    return path


def bridge_discovery(root, cfg):
    parents = pd.read_csv(root / 'NSCLC_automatic_ligand_response_review_20261016_v1/private_all_primary_BH_supported_rows.tsv', sep='\t')
    original = root / 'NSCLC_automatic_ligand_response_preflight_20261016_v1'
    original_plan = json.loads((original / 'plan_before_any_ligand_response_effect.json').read_text())
    assert len(original_plan['nuisance']) == 18
    names = original_plan['nuisance'][8:]
    assert len(names) == 10 and names[0] == 'sender_fraction' and names[-1] == 'max_qualified_receptor_joint_fraction'
    original_hashes = pd.read_csv(original / 'private_all_ready_design_SHA256.tsv', sep='\t').set_index('record').sha256
    tf_root = root / 'NSCLC_gene_disjoint_receiver_TF_measurement_20261016_v1'
    tf_scope = pd.read_csv(tf_root / 'private_all_TF_measurement_scope.tsv', sep='\t')
    tf_hashes = pd.read_csv(tf_root / 'private_ready_design_SHA256.tsv', sep='\t').set_index('TF_record').sha256
    result = []
    for parent in parents.itertuples():
        path = original / 'private_frozen_designs' / (parent.record+'.npz')
        assert sha(path) == original_hashes.loc[parent.record]
        a = np.load(path, allow_pickle=False)
        assert a['covariates'].ndim == 2 and a['covariates'].shape[1] == 18
        readouts = [('masked_programme', a['Y'])]
        if parent.ligand == cfg['secondary_ligand']:
            tf = tf_scope[tf_scope.record.eq(parent.record) & tf_scope.TF.eq(cfg['secondary_regulator'])].iloc[0]
            p = tf_root / 'private_frozen_designs' / (tf.TF_record+'.npz')
            assert sha(p) == tf_hashes.loc[tf.TF_record]
            data = np.load(p, allow_pickle=False)
            assert np.array_equal(data['X'], a['X']) and np.array_equal(data['cases'], a['cases'])
            readouts.append(('fixed_functional_followup', data['Y']))
        for kind, y in readouts:
            full = point(a['X'], y, a['covariates'])
            core = point(a['X'], y, a['covariates'][:, 8:])
            result.append(dict(original_record=parent.record, ligand=parent.ligand, readout=kind,
                               original_sender_author_subtype=parent.sender_author_subtype, original_receiver_author_subtype=parent.receiver_author_subtype,
                               cases=len(y), full18_partial_rank_rho=full['partial_rank_rho'], core10_partial_rank_rho=core['partial_rank_rho'],
                               no_new_discovery_P_or_BH=True, different_external_nodes_and_sampling_phase=True))
    return pd.DataFrame(result)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    pre = root / 'GSE305078_external_ligand_programme_core_20261016_v1'
    input_audit = wait_proof(root, pre.name + '/private_independent_input_audit.json', 'audit_GSE305078_external_core_inputs_v2')
    calibration = wait_proof(root, 'external_core_rank_model_calibration_20261016_v2/private_terminal.json', 'calibrate_external_core_rank_model_v2')
    kernel = Path(__file__).with_name('external_core_rank_model_v2.py')
    old_plan = json.loads((pre / 'plan_before_any_external_association.json').read_text())
    assert sha(kernel) == old_plan['kernel_sha256'] == json.loads(calibration.read_text())['kernel_sha256']
    cfg_file = root / 'external_ligand_programme_core_followup_config_20261016_v1.json'
    assert sha(cfg_file) == old_plan['config_sha256']
    cfg = json.loads(cfg_file.read_text())
    inp = json.loads((pre / 'private_input_terminal.json').read_text())
    for name, digest in inp['files'].items():
        assert sha(pre / name) == digest
    assert json.loads(input_audit.read_text())['producer_input_terminal_sha256'] == sha(pre / 'private_input_terminal.json')
    table = pd.read_csv(pre / 'private_all676_frozen_model_scope.tsv', sep='\t')
    hashes = pd.read_csv(pre / 'private_ready_design_SHA256.tsv', sep='\t').set_index('model_id').sha256
    out = root / 'GSE305078_external_ligand_programme_core_results_20261016_v1'
    assert not out.exists(); out.mkdir(mode=0o700); (out / 'private_resamples').mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), kernel_sha256=sha(kernel),
                input_audit_sha256=sha(input_audit), calibration_sha256=sha(calibration), frozen_config_sha256=sha(cfg_file),
                full_input_terminal_sha256=sha(pre / 'private_input_terminal.json'), complete_family=cfg['complete_family'],
                wild_draws=cfg['wild_draws'], bootstrap_draws=cfg['bootstrap_draws'], seed=cfg['seed'],
                inferential_units='Original source-stage matched physical libraries; no verified patient/outcome or exact-author-subtype promotion.',
                multiplicity='ALL676 predefined slots; unavailable internalP1, main BH, BY and family99% MonteCarlo-upper BH sensitivity. No later sample/TF/node/threshold replacements.',
                inference='HC3 studentized fixed restricted-HC2 wild approximation; no exact randomization or small-sample size guarantee. Full-control-space rank identified, original clinical covariates unavailable.',
                resampling='Pairs bootstrap2000, rank both outcomes afresh, complete control span; >=1800 valid required. All individual library deletions are descriptions.',
                discovery_bridge='All original primary-screen parents, same10 core controls alongside old18, descriptive only. Endpoint and pre/post-treatment differences explicitly retained.')
    (out / 'plan_before_true_external_effects.json').write_text(json.dumps(plan, indent=2) + '\n')
    bridge_discovery(root, cfg).to_csv(out / 'private_all_discovery_core10_bridge_points.tsv', sep='\t', index=False)
    rows = []; resample_hashes = []
    for i, row in enumerate(table.to_dict('records')):
        row.update(association_status='INPUT_UNAVAILABLE', partial_rank_rho=np.nan, HC3_t=np.nan,
                   wild_exceedances=np.nan, wild_draws=0, wild_P=np.nan, MC_upper_P=np.nan,
                   bootstrap_valid=0, bootstrap_low=np.nan, bootstrap_high=np.nan)
        if row['status'] == 'READY_FOR_INDEPENDENT_INPUT_AUDIT':
            p = pre / 'private_designs' / (row['model_id']+'.npz'); assert sha(p) == hashes.loc[row['model_id']]
            a = np.load(p, allow_pickle=False); x = a['X']; y = a['Y']; c = a['C_raw']
            seed = cfg['seed'] + i*10
            try:
                g = point(x, y, c, cfg['min_units'], cfg['min_residual_df'])
                count = wild_count(g, cfg['wild_draws'], seed)
                boot = bootstrap_rhos(x, y, c, cfg['bootstrap_draws'], seed+1)
                loo = leave_one_out(x, y, c)
                valid = int(np.isfinite(boot).sum()); lo, hi = np.quantile(boot[np.isfinite(boot)], [.025, .975]) if valid >= 1800 else [np.nan, np.nan]
                row.update(association_status='ESTIMATED_FIXED_FAMILY_WILD_APPROXIMATION', partial_rank_rho=g['partial_rank_rho'], coefficient=g['coefficient'],
                           HC3_t=g['HC3_t'], HC3_variance=g['HC3_variance'], wild_exceedances=count, wild_draws=cfg['wild_draws'],
                           wild_P=(1+count)/(cfg['wild_draws']+1), MC_upper_P=monte_carlo_upper(count, cfg['wild_draws'], cfg['complete_family']),
                           bootstrap_valid=valid, bootstrap_low=lo, bootstrap_high=hi, unadjusted_spearman=float(spearmanr(x,y).statistic),
                           LOO_valid=int(np.isfinite(loo).sum()), LOO_min=float(np.nanmin(loo)), LOO_max=float(np.nanmax(loo)),
                           LOO_direction_stable=bool(np.isfinite(loo).all() and np.all(np.sign(loo) == np.sign(g['partial_rank_rho']))))
                for name, extra in [('sender_programme_adjusted_rho', a['sender_programme']), ('sender_same_regulator_adjusted_rho', a['sender_same_regulator'])]:
                    row[name] = np.nan
                    if np.isfinite(extra).all():
                        try:
                            row[name] = point(x, y, np.column_stack([c, extra]), min_n=1, min_df=1, require_leverage=False)['partial_rank_rho']
                        except ValueError:
                            pass
                p = out / 'private_resamples' / (row['model_id']+'.npz')
                np.savez_compressed(p, bootstrap_rho=boot, LOO_rho=loo)
                resample_hashes.append(dict(model_id=row['model_id'], sha256=sha(p)))
            except ValueError as error:
                allowed = ('INSUFFICIENT_LIBRARY_UNITS', 'MISSING_INPUT_NOT_ZERO', 'CONSTANT_LIGAND',
                           'LIGAND_NOT_IDENTIFIABLE_GIVEN_COMPLETE_CONTROLS', 'LOW_EFFECTIVE_DF_OR_ILL_CONDITIONED',
                           'EXTREME_LEVERAGE', 'MISSING_OR_CONSTANT_RECEIVER_READOUT', 'NO_RECEIVER_VARIATION_AFTER_COMPLETE_CONTROLS',
                           'DEGENERATE_HC3_VARIANCE', 'NUMERICALLY_AMBIGUOUS_CONTROL_RANK', 'NUMERICALLY_AMBIGUOUS_AUGMENTED_RANK')
                if str(error) not in allowed:
                    raise
                row['association_status'] = 'UNESTIMABLE_' + str(error)
        rows.append(row)
        if (i+1) % 50 == 0 or i+1 == len(table):
            (out / 'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), processed=i+1, complete_family=len(table), estimated_so_far=len(resample_hashes))) + '\n')
    result = pd.DataFrame(rows)
    internal_p = result.wild_P.fillna(1).to_numpy(); upper_p = result.MC_upper_P.fillna(1).to_numpy()
    assert len(result) == cfg['complete_family']
    result['BH_q_complete676'] = multipletests(internal_p, method='fdr_bh')[1]
    result['BY_q_complete676'] = multipletests(internal_p, method='fdr_by')[1]
    result['MC_upper_BH_q_complete676'] = multipletests(upper_p, method='fdr_bh')[1]
    unavailable = result.wild_P.isna()
    result.loc[unavailable, ['BH_q_complete676', 'BY_q_complete676', 'MC_upper_BH_q_complete676']] = np.nan
    result.to_csv(out / 'private_all676_external_core_results.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(resample_hashes).to_csv(out / 'private_resample_SHA256.tsv', sep='\t', index=False)
    terminal = dict(status='COMPLETE_FULL676_EXTERNAL_CORE_INFERENCE_INDEPENDENT_AUDIT_PENDING_NOT_MECHANISM_OR_CLINICAL_PREDICTION',
                    finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), kernel_sha256=sha(kernel),
                    complete_family=len(result), estimated_models=len(resample_hashes), wild_draws_each=cfg['wild_draws'], bootstrap_draws_each=cfg['bootstrap_draws'],
                    files={p.name: sha(p) for p in out.iterdir() if p.is_file()}, no_endpoint_or_outcome_or_original_q_changes=True)
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n'); print(json.dumps(terminal))


if __name__ == '__main__':
    main()
