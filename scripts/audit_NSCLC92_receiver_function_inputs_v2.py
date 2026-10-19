"""Bind previous complete assay proof and independently verify added components."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC92_receiver_function_math_v2 import prepare, slots


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 1024**2)
    out = root / 'NSCLC92_receiver_function_input_audit_20261019_v2'
    assert not out.exists(); out.mkdir(mode=0o700)
    src = root / 'NSCLC92_receiver_function_inputs_20261019_v2'
    t = json.loads((src / 'private_terminal.json').read_text())
    for name, d in t['files'].items(): assert sha(src / name) == d['sha256']
    old = root / 'NSCLC92_receiver_function_inputs_20261019_v1'
    oldproof = root / 'NSCLC92_receiver_function_input_audit_20261019_v1/private_terminal.json'
    assert sha(oldproof) == t['old_input_audit_sha256'] == '9f756ac9ee50959fc1c840e7925fa763ea4d4094bddaa95d46b0d392ba959e7a'
    assert json.loads(oldproof.read_text())['producer_terminal_sha256'] == sha(old / 'private_terminal.json')
    a = dict(np.load(src / 'private_complete_function_inputs.npz'))
    baseline = dict(np.load(old / 'private_complete_function_inputs.npz'))
    for name, value in baseline.items(): np.testing.assert_array_equal(value, a[name])
    params = np.load(src / 'private_discovery_component_scaling.npz')['component_mean_SD']
    target = pd.read_csv(old / 'private_complete_old_BH_receiver_scope.tsv', sep='\t')
    vocab = json.loads((root / 'NSCLC_communication_tensor_family_20261019_v1/private_selected_LR_vocabulary.json').read_text())
    vocab = {v['canonical_resource_LR_id']: v for v in vocab}
    mean = np.load(root / 'NSCLC_full_native_gene_units_20261016_v1/mean_logCP10K.npy', mmap_mode='r')
    for k, row in enumerate(target.itertuples()):
        lr = vocab[row.canonical_resource_LR_id]
        for e, name in enumerate(['sender_expression_score', 'receiver_expression_score']):
            columns = lr['native_ligand_columns'] if e == 0 else lr['native_receptor_columns']
            valid = a['unit_cells'][:, k, e] >= 20
            values = np.full(216, np.nan)
            for i in np.flatnonzero(valid):
                components = [float(mean[a['unit_indices'][i, k, e], c]) for c in columns]
                values[i] = np.float32(sorted(components)[0])
            train = valid & (a['partition'] == 'DISCOVERY')
            mu = np.sum(values[train]) / train.sum()
            sd = np.sqrt(np.mean((values[train] - mu)**2))
            np.testing.assert_allclose([mu, sd], params[k, e], atol=1e-12, rtol=1e-12)
            expected = (values - mu) / sd if sd > 1e-8 else np.full(216, np.nan)
            np.testing.assert_allclose(expected, a[name][:, k], atol=1e-12, rtol=1e-12)
    coverage = pd.read_csv(src / 'private_complete_family_design_coverage.tsv', sep='\t').set_index('slot')
    for slot, k, f, partition, question in slots(a):
        m, g, status, valid = prepare(a, k, f, partition, question, independent=True)
        assert coverage.loc[slot, 'status'] == status
        assert coverage.loc[slot, 'author_records'] == valid.sum()
    proof = dict(status='ACCEPTED_AMENDED_FUNCTION_INPUTS_COMPLETE_COMPONENT_SCALES_AND_THREE_QUESTION_DESIGNS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), original_full_input_audit_sha256=sha(oldproof), original_assay_and_covariates_bitwise_unchanged=True, added_components_completely_reconstructed=True, all_new_deleted_QR_designs_checked=True, clinical_function_effects_read=False, resource_guard=guard, scientific_outputs_public=False)
    (out / 'private_terminal.json').write_text(json.dumps(proof, indent=2) + '\n')
    print(json.dumps(dict(status=proof['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
