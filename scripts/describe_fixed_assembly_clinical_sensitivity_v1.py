"""Fixed assembly follow-up; two descriptive effects, no renewed P/FDR."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_frozen92_axis_io_v2 import controls, eligibility
from run_gamma_receiver_identity_sensitivity_v1 import fit_pair


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 256 * 1024)
    qualification = root / 'native_fixed_ligand_assembly_qualification_20261019_v1/private_terminal.json'
    assert sha(qualification) == '6eb4d53a6d61a75253d3c3eb418fe2da9ecbfe0f5bcec9539a69c2326690f545'
    q = json.loads(qualification.read_text())
    request = root / 'private_complementary_spatial_interface_scope_request_20261019_v1.json'
    assert sha(request) == q['private_fixed_request_sha256']
    fixed = json.loads(request.read_text())
    inp = root / 'NSCLC_frozen92_axis_source_inputs_20261019_v2'
    for name, info in json.loads((inp / 'private_terminal.json').read_text())['files'].items():
        assert sha(inp / name) == info['sha256']
    a = dict(np.load(inp / 'private_all129_all92_native_axis_and_design_inputs.npz', allow_pickle=False))
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    native_plan = json.loads((qualification.parent / 'plan_before_cached_assembly_marginals.json').read_text())
    for name in ['private_native_gene_axis.tsv', 'mean_logCP10K.npy']:
        assert sha(native / name) == native_plan['current_inputs_sha256'][name]
    genes = pd.read_csv(native / 'private_native_gene_axis.tsv', sep='\t').gene.to_numpy(str)
    means = np.load(native / 'mean_logCP10K.npy', mmap_mode='r')
    out = root / 'fixed_assembly_clinical_sensitivity_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_assembly_followup_effect.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), assembly_qualification_sha256=sha(qualification), previous_fixed_request_sha256=sha(request), numeric_model_module_sha256=sha(Path(__file__).with_name('run_gamma_receiver_identity_sensitivity_v1.py')), scope='All previously fixed supported axes: original component vs specified assembly extension; post-result descriptive only. Original endpoint cells,129 records,104 original-complete records where applicable, full controls, old mean/SD, same blocks. No new P or FDR, no clinical candidate promotion.', score_arrays_persisted=False, resource_guard=guard))
    rows = []
    for target in fixed['targets']:
        axis = target['axis_id']
        j = a['axis_ids'].tolist().index(axis)
        valid, _, _ = eligibility(a, j)
        C, names = controls(a, j)
        si = a['original_unit_indices'][:, int(a['sender_index'][j])]
        ri = a['original_unit_indices'][:, int(a['receiver_index'][j])]
        strict = list(dict.fromkeys(target['original_ligand_genes'] + target['strict_ligand_assembly_addition']))
        lc = [int(np.where(genes == g)[0][0]) for g in strict]
        rc = [int(np.where(genes == g)[0][0]) for g in target['original_receptor_genes']]
        raw = np.full(129, np.nan)
        for i in np.flatnonzero(a['endpoint_measurable'][:, j]):
            L = np.float32(min(float(means[si[i], c]) for c in lc))
            R = np.float32(min(float(means[ri[i], c]) for c in rc))
            raw[i] = np.float32(np.log1p(float(L) * float(R)))
        assert np.all(raw[valid] <= a['raw_axis_score'][valid, j] + 1e-12)
        mu, sd = a['blind_score_mean_SD'][j]
        for label, y in [('ORIGINAL_COMPONENT', a['standardized_axis_score'][:, j]), ('SPECIFIED_ASSEMBLY_EXTENSION', (raw - mu) / sd)]:
            effect = fit_pair(a['response'], y, C, names, a['profile_blocks'], valid)
            rows.append(dict(source='GSE243013', axis_id=axis, ligand_proxy=label, **effect))
    assert len(rows) == 2 * len(fixed['targets'])
    pd.DataFrame(rows).to_csv(out / 'private_complete2_fixed_assembly_followup_effects.tsv', sep='\t', index=False)
    receipt = dict(status='COMPLETE_FIXED_ASSEMBLY_SENSITIVITY_ORIGINAL_DESIGN_AND_SEPARATE_QR_NO_NEW_P_FDR', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), assembly_qualification_sha256=sha(qualification), descriptive_positions=len(rows), all_estimable_models_separate_nuisance_and_deleted_QR_checked=True, new_raw_expression_read_or_downloaded=False, full_expression_or_score_arrays_persisted=False, original_inference_changed=False, new_P_or_FDR=False, scientific_outputs_public=False, resource_guard=guard, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in receipt['files'].values()) < 256 * 1024
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
