"""Independent all-slot input reconstruction; never calculate an X-Y effect."""
import hashlib
import itertools
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from scipy.linalg import qr


NAMES = ['sender_fraction', 'sender_logN', 'sender_mean_logRNA', 'sender_mean_loggenes',
         'receiver_fraction', 'receiver_logN', 'receiver_mean_logRNA', 'receiver_mean_loggenes',
         'receiver_own_ligand', 'max_qualified_receptor_joint_fraction']


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'GSE305078_external_ligand_programme_core_20261016_v1'
    receipt = json.loads((root / 'HNSCC_receiver_followup_job_control_20261016_v2/preflight_GSE305078_external_core_v1/launch_receipt.json').read_text())
    deadline = time.monotonic() + 3600
    while not (out / 'private_input_terminal.json').exists():
        p = Path('/proc') / str(receipt['pid']) / 'stat'
        assert p.exists()
        s = p.read_text().rsplit(')', 1)[1].split()
        assert s[19] == str(receipt['proc_start_ticks']) and s[0] != 'Z' and time.monotonic() < deadline
        time.sleep(5)
    terminal = json.loads((out / 'private_input_terminal.json').read_text())
    for name, digest in terminal['files'].items():
        assert sha(out / name) == digest
    cfg_file = root / 'external_ligand_programme_core_followup_config_20261016_v1.json'
    cfg = json.loads(cfg_file.read_text())
    plan = json.loads((out / 'plan_before_any_external_association.json').read_text())
    assert sha(cfg_file) == plan['config_sha256']
    scope = pd.read_csv(out / 'private_all676_frozen_model_scope.tsv', sep='\t', keep_default_na=False)
    member = pd.read_csv(out / 'private_all_model_library_membership.tsv.gz', sep='\t', dtype={'library_id': str})
    hashes = pd.read_csv(out / 'private_ready_design_SHA256.tsv', sep='\t').set_index('model_id').sha256
    native = root / 'existing_immune_fine_native_signatures_20261016_v4'
    unit = pd.read_csv(native / 'private_all_original_fine_units_and_quality.tsv', sep='\t', dtype={'library_id': str})
    unit = unit[unit.dataset.eq(cfg['source'])].copy()
    denominator = unit.groupby(['library_id', 'major']).n_cells.sum()
    assert np.allclose(unit.within_original_major_fraction, [r.n_cells / denominator.loc[(r.library_id, r.major)] for r in unit.itertuples()])
    side = pd.read_csv(root / 'GSE305078_official_sampling_sidecar_20261016_v2/private_exact_official_sampling_sidecar.tsv', sep='\t', dtype={'library_id': str})
    opportunity = root / 'external_fine_pair_case_opportunities_20261016_v1'
    families = pd.read_csv(opportunity / 'private_all_deduplicated_extension_families.tsv', sep='\t').set_index('extension_family')
    pairs = pd.read_csv(opportunity / 'private_all_source_fine_pair_coverage.tsv', sep='\t')
    pairs = pairs[pairs.dataset.eq(cfg['source'])].sort_values(['extension_family', 'sender_node', 'receiver_node'])
    pm = pd.read_csv(native / 'private_all_source_parent_programme_scope.tsv.gz', sep='\t')
    pm = pm[pm.dataset.eq(cfg['source'])]
    tm = pd.read_csv(native / 'private_all_source_parent_TF_scope.tsv.gz', sep='\t')
    tm = tm[tm.dataset.eq(cfg['source'])]
    axis = pd.read_csv(native / 'private_all_source_panel_gene_axis.tsv.gz', sep='\t')
    axis = axis[axis.dataset.eq(cfg['source'])].set_index('gene')
    a = np.load(native / 'private_frozen_programme_and_TF_scores.npz'); P = a['programme']; T = a['TF']
    a = np.load(native / 'private_limited_native_panel.npz'); L = a['mean_logCP10K']; C = a['positive_cell_counts']
    expected = []
    for stage, pair in itertools.product(cfg['stages'], pairs.itertuples()):
        kinds = ['masked_programme'] + (['fixed_functional_followup'] if families.loc[pair.extension_family, 'ligand'] == cfg['secondary_ligand'] else [])
        for kind in kinds:
            expected.append((stage, pair.extension_family, pair.sender_node, pair.receiver_node, kind))
    assert len(expected) == len(scope) == 676
    assert expected == list(scope[['stage', 'extension_family', 'sender_node', 'receiver_node', 'readout']].itertuples(index=False, name=None))
    ready_count = 0
    for row in scope.itertuples():
        f = families.loc[row.extension_family]; parents = f.original_parent_records.split('|')
        ids = set(side.loc[side.official_stage.eq(row.stage), 'library_id'])
        sub = unit[unit.library_id.isin(ids) & unit.n_cells.ge(20)]
        ss = sub[sub.node_id.eq(row.sender_node)].set_index('library_id')
        rr = sub[sub.node_id.eq(row.receiver_node)].set_index('library_id')
        lib = sorted(set(ss.index) & set(rr.index))
        assert row.library_units == len(lib) and row.original_parent_records == f.original_parent_records
        membership = member[member.model_id.eq(row.model_id)]
        assert membership.library_id.tolist() == lib
        if len(lib):
            assert np.array_equal(membership.sender_unit, ss.loc[lib, 'unit_index'])
            assert np.array_equal(membership.receiver_unit, rr.loc[lib, 'unit_index'])
        if len(lib) < cfg['min_units']:
            assert row.status == 'INSUFFICIENT_LIBRARY_UNITS'; continue
        si = ss.loc[lib, 'unit_index'].to_numpy(int); ri = rr.loc[lib, 'unit_index'].to_numpy(int)
        lc = int(axis.loc[f.ligand, 'panel_column'])
        x = L[si, lc]; own = L[ri, lc]
        assert int(np.sum(x > 0)) == row.ligand_positive_libraries
        if np.sum(x > 0) < 10:
            assert row.status == 'INSUFFICIENT_LIGAND_POSITIVE_LIBRARIES'; continue
        receptor_columns = []
        for receptor in f.qualified_receptors.split('|'):
            assert '_' not in receptor and axis.loc[receptor, 'source_native_unambiguous']
            receptor_columns.append(int(axis.loc[receptor, 'panel_column']))
        fractions = C[np.ix_(ri, receptor_columns)] / rr.loc[lib, 'n_cells'].to_numpy()[:, None]
        pool = C[np.ix_(ri, receptor_columns)].sum(axis=0) / rr.loc[lib, 'n_cells'].sum()
        assert np.isclose(row.maximum_original_receptor_pooled_RNA_fraction, pool.max())
        if pool.max() < .01:
            assert row.status == 'NO_ORIGINAL_RECEPTOR_ONE_PERCENT_POOLED_RNA'; continue
        roles = {role: int(pm.loc[pm.record.eq(parents[0]) & pm.role.eq(role), 'score_column'].iloc[0]) for role in ['sender', 'receiver_masked']}
        sender_programme = P[si, roles['sender']]
        same = np.full(len(lib), np.nan)
        if row.readout == 'masked_programme':
            y = P[ri, roles['receiver_masked']]
        else:
            frow = tm[tm.record.eq(parents[0]) & tm.TF.eq(cfg['secondary_regulator'])].iloc[0]
            if not frow.status.startswith('MEASURED_') or frow.target_retained_fraction < .95 or frow.background_retained_fraction < .95:
                assert row.status == 'INSUFFICIENT_MATCHED_FUNCTIONAL_TARGET_OR_BACKGROUND'; continue
            y = T[ri, int(frow.score_column)]; same = T[si, int(frow.score_column)]
        cv = np.array([[getattr(g, name) if name != 'logN' else np.log1p(g.n_cells)
                        for g in [ss.loc[b], rr.loc[b]] for name in ['within_original_major_fraction', 'logN', 'mean_logRNA', 'mean_loggenes']]
                       + [own[j], fractions[j].max()] for j, b in enumerate(lib)], dtype=float)
        sd = cv.std(axis=0); z = np.column_stack([np.ones(len(lib)), (cv-cv.mean(axis=0))/np.where(sd > 1e-12, sd, 1)])
        u, s, _ = np.linalg.svd(z, full_matrices=False); rank_z = int((s > s[0]*1e-12).sum())
        xr = rankdata(x).astype(float); xr = (xr-xr.mean())/xr.std() if xr.std() else xr*0
        d = np.column_stack([z, xr]); singular = np.linalg.svd(d, compute_uv=False); rank_d = int((singular > singular[0]*1e-12).sum())
        condition = singular[0]/singular[rank_d-1]
        if len(lib)-rank_d < 10 or condition >= 1e8:
            assert row.status == 'LOW_EFFECTIVE_DF_OR_ILL_CONDITIONED'; continue
        assert rank_d == rank_z+1
        h = np.sum(d * np.linalg.pinv(d, rcond=1e-12).T, axis=1)
        if h.max() >= .999999:
            assert row.status == 'EXTREME_LEVERAGE'; continue
        assert row.status == 'READY_FOR_INDEPENDENT_INPUT_AUDIT'
        q, _, _ = qr(z, mode='economic', pivoting=True)
        residual_y = rankdata(y) - q[:, :rank_z] @ (q[:, :rank_z].T @ rankdata(y))
        assert np.linalg.norm(residual_y) > 1e-8
        assert int(row.effective_rank) == rank_d and int(row.nuisance_rank) == rank_z and int(row.residual_df) == len(lib)-rank_d
        assert bool(row.strict_full_column_rank) == (rank_d == 12)
        path = out / 'private_designs' / (row.model_id+'.npz')
        assert sha(path) == hashes.loc[row.model_id]
        saved = np.load(path, allow_pickle=False)
        assert saved['covariate_names'].astype(str).tolist() == NAMES
        for name, value in [('X', x), ('Y', y), ('C_raw', cv), ('sender_programme', sender_programme), ('sender_same_regulator', same)]:
            assert np.allclose(saved[name], value, atol=1e-12, rtol=1e-12, equal_nan=True), name
        assert saved['library_ids'].astype(str).tolist() == lib
        ready_count += 1
    assert ready_count == len(hashes)
    proof = dict(status='PASS_ALL676_SCOPE_LIBRARY_MEMBERSHIP_NATIVE_COLUMNS_FROZEN_READOUTS_COMPLETE_CONTROLS_AND_READY_DESIGNS_NO_EFFECTS',
                 finished_UTC=datetime.now(timezone.utc).isoformat(), audit_source_sha256=sha(__file__), producer_input_terminal_sha256=sha(out / 'private_input_terminal.json'),
                 complete_family=len(scope), ready_designs=ready_count, no_external_associations_P_BH_or_outcomes=True)
    p = out / 'private_independent_input_audit.json'; assert not p.exists()
    p.write_text(json.dumps(proof, indent=2) + '\n'); print(json.dumps(proof))


if __name__ == '__main__':
    main()
