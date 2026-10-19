"""Pre-effect addition: sender RNA conditional on receiver RNA and response."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC92_receiver_function_math_v2 import slots, prepare


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 1024**2)
    out = root / 'NSCLC92_receiver_function_inputs_20261019_v2'
    assert not out.exists(); out.mkdir(mode=0o700)
    old = root / 'NSCLC92_receiver_function_inputs_20261019_v1'
    proof = root / 'NSCLC92_receiver_function_input_audit_20261019_v1/private_terminal.json'
    assert sha(proof) == '9f756ac9ee50959fc1c840e7925fa763ea4d4094bddaa95d46b0d392ba959e7a'
    assert json.loads(proof.read_text())['producer_terminal_sha256'] == sha(old / 'private_terminal.json')
    assert not (root / 'NSCLC92_receiver_function_family_20261019_v1').exists()
    save = lambda n, d: (out / n).write_text(json.dumps(d, indent=2) + '\n')
    save('plan_before_component_expression_and_any_true_effects.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), old_full_input_audit_sha256=sha(proof), resource_guard=guard, reason='Before any functional effect: axis-function association may be receptor cell-state proxy. Add sender expression target controlling receiver expression and original response, plus all original quality/abundance/clinical controls. This is not causal identification.', amended_family='All50 x2original partitions x3questions xall oldBH receiver scopes; jointBH/BY with failures, no new reference/endpoints/seed search.', component_scales='Percomponent original129 all own endpoint>=20 records without opposite endpoint/response filtering; populationSD, original mincomponent float32 meanlog expression. Apply unchanged to original87.', all_true_function_effects_unread=True))
    t = json.loads((old / 'private_terminal.json').read_text())
    for name, d in t['files'].items(): assert sha(old / name) == d['sha256']
    a = dict(np.load(old / 'private_complete_function_inputs.npz'))
    targets = pd.read_csv(old / 'private_complete_old_BH_receiver_scope.tsv', sep='\t')
    vocab = json.loads((root / 'NSCLC_communication_tensor_family_20261019_v1/private_selected_LR_vocabulary.json').read_text())
    vocab = {v['canonical_resource_LR_id']: v for v in vocab}
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    mean = np.load(native / 'mean_logCP10K.npy', mmap_mode='r')
    standardized = np.full((216, len(targets), 2), np.nan)
    parameters = np.full((len(targets), 2, 2), np.nan)
    for k, target in enumerate(targets.itertuples()):
        lr = vocab[target.canonical_resource_LR_id]
        for e, columns in enumerate([lr['native_ligand_columns'], lr['native_receptor_columns']]):
            measured = a['unit_cells'][:, k, e] >= 20
            raw = np.full(216, np.nan)
            for i in np.flatnonzero(measured): raw[i] = float(np.float32(np.min(mean[a['unit_indices'][i, k, e], columns])))
            train = measured & (a['partition'] == 'DISCOVERY')
            mu, sd = float(raw[train].mean()), float(raw[train].std(ddof=0))
            parameters[k, e] = [mu, sd]
            if sd > 1e-8: standardized[measured, k, e] = (raw[measured] - mu) / sd
    a['sender_expression_score'] = standardized[:, :, 0]
    a['receiver_expression_score'] = standardized[:, :, 1]
    np.savez_compressed(out / 'private_complete_function_inputs.npz', **a)
    np.savez_compressed(out / 'private_discovery_component_scaling.npz', component_mean_SD=parameters)
    coverage = []
    for slot, k, f, partition, question in slots(a):
        m, g, status, valid = prepare(a, k, f, partition, question)
        coverage.append(dict(slot=slot, scope_id=a['scope_ids'][k], reference=a['references'][f], partition=partition, question=question, status=status, author_records=int(valid.sum()), nonMPR=int((a['response'][valid] == 0).sum()), pCR_or_MPR=int((a['response'][valid] == 1).sum()), profile_blocks=len(set(a['profile_blocks'][valid])), residual_df=m['df'] if m else np.nan, max_block_leverage=g['max_block_leverage'] if g else np.nan))
    pd.DataFrame(coverage).to_csv(out / 'private_complete_family_design_coverage.tsv', sep='\t', index=False)
    files = {p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}
    assert sum(d['bytes'] for d in files.values()) < 1024**2
    save('private_terminal.json', dict(status='COMPLETE_AMENDED_THREE_QUESTION_INPUTS_NO_TRUE_FUNCTION_EFFECTS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), old_input_audit_sha256=sha(proof), files=files, true_functional_effects_read=False, scientific_outputs_public=False))
    print(json.dumps(dict(status='COMPLETE_AMENDED_THREE_QUESTION_INPUTS_NO_TRUE_FUNCTION_EFFECTS', terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
