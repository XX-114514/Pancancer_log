"""Independently reconstruct every new pooled expression score and design."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON
from NSCLC92_receiver_function_math_v1 import slots, prepare


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 2 * 1024**2)
    src = root / 'NSCLC92_receiver_function_inputs_20261019_v1'
    out = root / 'NSCLC92_receiver_function_input_audit_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    t = json.loads((src / 'private_terminal.json').read_text())
    for name, d in t['files'].items(): assert sha(src / name) == d['sha256']
    a = dict(np.load(src / 'private_complete_function_inputs.npz'))
    p = dict(np.load(src / 'private_discovery_only_gene_and_function_scaling.npz'))
    target = pd.read_csv(src / 'private_complete_old_BH_receiver_scope.tsv', sep='\t')
    defs = json.loads((root / 'NSCLC92_receiver_function_reference_scope_20261019_v1/private_all50_LR_disjoint_reference_definitions.json').read_text())
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    u = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t')
    index = pd.MultiIndex.from_frame(u[['sampleID', 'sub_cell_type']])
    mean = np.load(native / 'mean_logCP10K.npy', mmap_mode='r')
    original = root / 'NSCLC_communication_tensor_family_20261019_v1'
    tensor = dict(np.load(original / 'private_tensor_axes.npz'))
    assert a['record_ids'].tolist() == tensor['discovery_ids'].tolist() + tensor['heldout_ids'].tolist()
    blind = pd.read_csv(root / 'NSCLC_communication_tensor_inputs_20261019_v1/private_all243_blind_record_partitions.tsv', sep='\t', keep_default_na=False).set_index('sampleID').loc[a['record_ids']]
    assert np.array_equal(blind.partition, a['partition']) and np.array_equal(blind.profile_block, a['profile_blocks'])
    clinical_path = root / 'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv'
    old_terminal = json.loads((root / 'NSCLC_frozen92_axis_source_inputs_20261019_v2/private_terminal.json').read_text())
    assert sha(clinical_path) == old_terminal['clinical_source_sha256']
    clinical = pd.read_csv(clinical_path, sep='\t', keep_default_na=False).set_index('sampleID').loc[a['record_ids']]
    response = pd.to_numeric(clinical.response, errors='coerce').to_numpy(float)
    np.testing.assert_array_equal(a['response'], response)
    levels = np.load(root / 'NSCLC_frozen92_axis_source_inputs_20261019_v2/private_all129_all92_native_axis_and_design_inputs.npz')['all_agent_levels']
    base = np.column_stack([clinical[COMMON].apply(pd.to_numeric, errors='coerce'), [[float(value == level) for level in levels[1:]] for value in blind.agent_identity]])
    total = u.groupby('sampleID').n_cells.sum()
    globalQ = np.column_stack([np.log1p(total.reindex(a['record_ids']))] + [(u.assign(v=u.n_cells * u[f]).groupby('sampleID').v.sum() / total).reindex(a['record_ids']) for f in ['mean_logRNA', 'mean_loggenes']])
    parent = u.groupby(['sampleID', 'major_cell_type']).n_cells.sum()
    vocab = json.loads((original / 'private_selected_LR_vocabulary.json').read_text())
    vocab = {v['canonical_resource_LR_id']: v for v in vocab}
    old = dict(np.load(root / 'NSCLC_frozen92_axis_source_inputs_20261019_v2/private_all129_all92_native_axis_and_design_inputs.npz'))
    maximum = 0.
    for k, row in enumerate(target.itertuples()):
        states = [row.sender_original_author_type, row.receiver_original_author_type]
        ids = np.stack([index.get_indexer(pd.MultiIndex.from_tuples([(case, state) for case in a['record_ids']])) for state in states], axis=1)
        np.testing.assert_array_equal(ids, a['unit_indices'][:, k])
        n = np.zeros_like(ids); Q = np.full((216, 2, 2), np.nan); fraction = np.full((216, 2), np.nan)
        for i, case in enumerate(a['record_ids']):
            for e in range(2):
                if ids[i, e] >= 0:
                    rr = u.iloc[ids[i, e]]; n[i, e] = int(rr.n_cells)
                    Q[i, e] = rr[['mean_logRNA', 'mean_loggenes']]
                    fraction[i, e] = n[i, e] / parent.loc[(case, rr.major_cell_type)]
        np.testing.assert_array_equal(n, a['unit_cells'][:, k])
        np.testing.assert_array_equal((n >= 20).all(1), a['endpoint_measurable'][:, k])
        control = np.column_stack([base, globalQ, np.log1p(n[:, 0]), np.log1p(n[:, 1]), Q[:, 0], Q[:, 1], fraction[:, 0], fraction[:, 1]])
        np.testing.assert_allclose(control, a['controls'][:, k], atol=1e-12, rtol=1e-12)
        usable = n[:, 1] >= 20
        train = usable & (a['partition'] == 'DISCOVERY')
        local = {int(c): z for z, c in enumerate(p['native_columns'])}
        X = np.full((216, len(local)), np.nan)
        X[usable] = mean[np.ix_(ids[usable, 1], p['native_columns'])]
        center = np.mean(X[train], axis=0); scale = np.std(X[train], axis=0)
        np.testing.assert_allclose(center, p['gene_mean'][k], atol=1e-12, rtol=1e-12)
        np.testing.assert_allclose(scale, p['gene_SD'][k], atol=1e-12, rtol=1e-12)
        for f, d in enumerate(defs):
            cols = [local[c] for c in d['unique_native_columns'].values()]
            variable = [c for c in cols if scale[c] > 1e-8]
            eligible = bool(d['complement_measurement_eligible'] and len(variable) >= 20 and len(variable) >= .8 * len(cols))
            expected = np.full(216, np.nan)
            if eligible:
                # Gene-by-gene scalar summation, separately centred at reference-only records.
                value = np.zeros(216)
                for c in variable: value += (X[:, c] - center[c]) / scale[c]
                value /= len(variable)
                mu, sd = np.mean(value[train]), np.std(value[train])
                eligible = bool(np.isfinite(sd) and sd > 1e-8)
                np.testing.assert_allclose([mu, sd], p['function_mean_SD'][k, f], atol=1e-12, rtol=1e-12)
                if eligible: expected = (value - mu) / sd
            assert a['function_eligible'][k, f] == eligible
            np.testing.assert_allclose(expected, a['function_score'][k, :, f], atol=1e-10, rtol=1e-10)
            if eligible: maximum = max(maximum, float(np.nanmax(np.abs(expected - a['function_score'][k, :, f]))))
        lr = vocab[row.canonical_resource_LR_id]
        j = old['axis_ids'].tolist().index(row.axis_id)
        axis = np.full(216, np.nan)
        for i in np.flatnonzero(a['endpoint_measurable'][:, k]):
            left = np.min(mean[ids[i, 0], lr['native_ligand_columns']]).astype(np.float32)
            right = np.min(mean[ids[i, 1], lr['native_receptor_columns']]).astype(np.float32)
            raw = np.float32(np.log1p(np.float64(left) * np.float64(right)))
            axis[i] = (float(raw) - old['blind_score_mean_SD'][j, 0]) / old['blind_score_mean_SD'][j, 1]
        np.testing.assert_allclose(axis, a['axis_score'][:, k], atol=0, rtol=0)
    coverage = pd.read_csv(src / 'private_complete_family_design_coverage.tsv', sep='\t').set_index('slot')
    rows = []
    for slot, k, f, partition, question in slots(a):
        m, g, status, valid = prepare(a, k, f, partition, question, independent=True)
        assert coverage.loc[slot, 'status'] == status
        assert coverage.loc[slot, 'author_records'] == valid.sum()
        rows.append(dict(slot=slot, status=status, author_records=int(valid.sum())))
    pd.DataFrame(rows).to_csv(out / 'private_all_independent_design_checks.tsv', sep='\t', index=False)
    proof = dict(status='ACCEPTED_COMPLETE_FUNCTION_INPUTS_ALL_NATIVE_VALUES_SCALES_CONTROLS_AND_DELETED_QR', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), maximum_score_reconstruction_difference=maximum, full_source_raw_replayed=False, inherited_raw_audit_units=64, clinical_function_effects_read=False, resource_guard=guard, scientific_outputs_public=False)
    (out / 'private_terminal.json').write_text(json.dumps(proof, indent=2) + '\n')
    print(json.dumps(dict(status=proof['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
