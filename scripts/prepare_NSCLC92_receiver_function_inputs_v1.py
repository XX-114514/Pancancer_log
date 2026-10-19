"""All-reference native expression assay and fixed two-partition inputs."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON
from NSCLC_frozen92_axis_io_v2 import controls as old_controls
from NSCLC92_receiver_function_math_v1 import prepare, slots


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 16 * 1024**2)
    out = root / 'NSCLC92_receiver_function_inputs_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda n, d: (out / n).write_text(json.dumps(d, indent=2) + '\n')
    ref = root / 'NSCLC92_receiver_function_reference_scope_20261019_v1'
    proof = root / 'NSCLC92_receiver_reference_scope_audit_20261019_v1/private_terminal.json'
    audit = json.loads(proof.read_text())
    assert audit['status'] == 'ACCEPTED_ALL50_REFERENCE_NATIVE_MAPPING_EXCLUSIONS_AND_COVERAGE'
    assert audit['producer_terminal_sha256'] == sha(ref / 'private_terminal.json')
    result_proof = root / 'NSCLC_frozen92_axis_source_result_audit_20261019_v1/private_terminal.json'
    assert sha(result_proof) == '79e8aa29aec12db05d437f6e4ebc24dec7e97772bf62ca594434efe36522b9b7'
    save('plan_before_all_reference_expression_and_effects.json', dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        reference_audit_sha256=sha(proof), total_resource_guard=guard,
        scoring='All50 LRdisjoint complements. Gene mean/populationSD from all original129 receiver>=20 units without outcome/sender filtering. SD>1e-8. >=20variable genes and>=.8 of uniquely mapped complement. Mean genez then discovery-only programme SD. Apply unchanged to original87.',
        family='All original BH92 receivers x50references x2original partitions x2questions, joint BH/BY including failures. WCR39999/WCU19999, frozen methods. All new functional effects unread.',
        raw_proof='Inherited wholecache conservation/overlap and64unit allnativegene independent raw audit, not all native units. No new download or full raw scan.',
        all_scientific_outputs_private=True))
    registry = pd.read_csv(root / 'NSCLC_frozen92_axis_source_inputs_20261019_v2/private_all92_frozen_axis_registry.tsv', sep='\t')
    result_dir = root / 'NSCLC_frozen92_axis_source_consistency_20261019_v1'
    result_file = result_dir / 'private_all92_source_consistency_results.tsv'
    assert sha(result_file) == json.loads((result_dir / 'private_terminal.json').read_text())['files'][result_file.name]['sha256']
    results = pd.read_csv(result_file, sep='\t')
    assert sha(root / 'NSCLC_frozen92_axis_source_consistency_20261019_v1/private_terminal.json') == json.loads(result_proof.read_text())['producer_terminal_sha256']
    selected = registry[registry.axis_id.isin(results.loc[results.BH92 < .05, 'axis_id'])].copy()
    assert len(selected) > 0
    selected.to_csv(out / 'private_complete_old_BH_receiver_scope.tsv', sep='\t', index=False)
    original = root / 'NSCLC_communication_tensor_family_20261019_v1'
    t = json.loads((original / 'private_terminal.json').read_text())
    for n in ['private_tensor_axes.npz', 'private_selected_LR_vocabulary.json']:
        assert sha(original / n) == t['files'][n]['sha256']
    axes = dict(np.load(original / 'private_tensor_axes.npz'))
    cases = axes['discovery_ids'].tolist() + axes['heldout_ids'].tolist()
    assert len(cases) == len(set(cases)) == 216
    vocab = json.loads((original / 'private_selected_LR_vocabulary.json').read_text())
    vi = {v['canonical_resource_LR_id']: i for i, v in enumerate(vocab)}
    old = dict(np.load(root / 'NSCLC_frozen92_axis_source_inputs_20261019_v2/private_all129_all92_native_axis_and_design_inputs.npz'))
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    units = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t')
    index = {(v.sampleID, v.sub_cell_type): int(v.unit_index) for v in units.itertuples()}
    mean = np.load(native / 'mean_logCP10K.npy', mmap_mode='r')
    blind = pd.read_csv(root / 'NSCLC_communication_tensor_inputs_20261019_v1/private_all243_blind_record_partitions.tsv', sep='\t', keep_default_na=False).set_index('sampleID').loc[cases]
    assert blind.partition.tolist() == ['DISCOVERY'] * 129 + ['HELDOUT'] * 87
    for name, item in json.loads((ref / 'private_terminal.json').read_text())['files'].items():
        assert sha(ref / name) == item['sha256']
    definitions = json.loads((ref / 'private_all50_LR_disjoint_reference_definitions.json').read_text())
    columns = sorted(set(c for d in definitions for c in d['unique_native_columns'].values()))
    reverse = {c: i for i, c in enumerate(columns)}
    K, F, N = len(selected), 50, 216
    function = np.full((K, N, F), np.nan)
    eligible = np.zeros((K, F), bool)
    centers = np.full((K, len(columns)), np.nan)
    scales = np.full_like(centers, np.nan)
    function_parameters = np.full((K, F, 2), np.nan)
    endpoint = np.zeros((N, K), bool)
    unit_ids = np.full((N, K, 2), -1, int)
    ncell = np.zeros_like(unit_ids)
    ownquality = np.full((N, K, 2, 2), np.nan)
    fractions = np.full((N, K, 2), np.nan)
    axis_score = np.full((N, K), np.nan)
    total = units.groupby('sampleID').n_cells.sum()
    global_quality = np.column_stack([np.log1p(total.reindex(cases))] + [
        (units.assign(w=units.n_cells * units[v]).groupby('sampleID').w.sum() / total).reindex(cases).to_numpy()
        for v in ['mean_logRNA', 'mean_loggenes']])
    parent = units.groupby(['sampleID', 'major_cell_type']).n_cells.sum()
    rows = []
    for k, target in enumerate(selected.itertuples()):
        expression = np.full((N, len(columns)), np.nan)
        for i, case in enumerate(cases):
            for e, state in enumerate([target.sender_original_author_type, target.receiver_original_author_type]):
                unit = index.get((case, state), -1)
                if unit >= 0:
                    row = units.iloc[unit]
                    unit_ids[i, k, e] = unit
                    ncell[i, k, e] = int(row.n_cells)
                    ownquality[i, k, e] = row[['mean_logRNA', 'mean_loggenes']].to_numpy(float)
                    fractions[i, k, e] = row.n_cells / parent.loc[(case, row.major_cell_type)]
            if ncell[i, k, 1] >= 20:
                expression[i] = mean[unit_ids[i, k, 1], columns]
        training = np.arange(N) < 129
        training &= ncell[:, k, 1] >= 20
        assert training.sum() >= 50
        centers[k] = expression[training].mean(0)
        scales[k] = expression[training].std(0, ddof=0)
        z = (expression - centers[k]) / np.where(scales[k] > 1e-8, scales[k], 1)
        for f, d in enumerate(definitions):
            mapped = [reverse[c] for c in d['unique_native_columns'].values()]
            varying = [c for c in mapped if scales[k, c] > 1e-8]
            good = bool(d['complement_measurement_eligible'] and len(varying) >= 20 and len(varying) >= .8 * len(mapped))
            if good:
                raw = z[:, varying].mean(1)
                center, scale = float(raw[training].mean()), float(raw[training].std(ddof=0))
                function_parameters[k, f] = [center, scale]
                good = np.isfinite(scale) and scale > 1e-8
                if good:
                    function[k, :, f] = (raw - center) / scale
            eligible[k, f] = good
            rows.append(dict(scope_id=target.axis_id, reference=d['reference'], original_genes=len(d['full_reference_genes']), mapped_complement_genes=len(mapped), varying_complement_genes=len(varying), reference_strong_representation=d['strong_reference_representation'], function_measurement_eligible=good))
        endpoint[:, k] = (ncell[:, k] >= 20).all(1)
        lr = vocab[vi[target.canonical_resource_LR_id]]
        j = old['axis_ids'].tolist().index(target.axis_id)
        mu, sd = old['blind_score_mean_SD'][j]
        for i in np.flatnonzero(endpoint[:, k]):
            ligand = np.float32(min(mean[unit_ids[i, k, 0], c] for c in lr['native_ligand_columns']))
            receptor = np.float32(min(mean[unit_ids[i, k, 1], c] for c in lr['native_receptor_columns']))
            raw = float(np.float32(np.log1p(float(ligand) * float(receptor))))
            axis_score[i, k] = (raw - mu) / sd
        np.testing.assert_allclose(axis_score[:129, k], old['standardized_axis_score'][:, j], atol=0, rtol=0)
    # All new gene scales fixed before joining clinical outcomes.
    save('private_before_new_function_outcome_join.json', dict(observed_UTC=datetime.now(timezone.utc).isoformat(), gene_and_function_scales_fixed=True, functional_effects_read=False))
    clinical = pd.read_csv(root / 'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv', sep='\t', keep_default_na=False).set_index('sampleID').loc[cases]
    levels = old['all_agent_levels'].tolist()
    agent = np.asarray([[float(v == level) for level in levels[1:]] for v in blind.agent_identity])
    base = np.column_stack([clinical[COMMON].apply(pd.to_numeric, errors='coerce').to_numpy(float), agent])
    names = COMMON + ['agent:' + v for v in levels[1:]] + ['logN_all_immune', 'all_immune_mean_logRNA', 'all_immune_mean_loggenes', 'logN_sender', 'logN_receiver', 'sender_meanlogRNA', 'sender_meanloggenes', 'receiver_meanlogRNA', 'receiver_meanloggenes', 'sender_within_major_fraction', 'receiver_within_major_fraction']
    control = np.stack([np.column_stack([base, global_quality, np.log1p(ncell[:, k, 0]), np.log1p(ncell[:, k, 1]), ownquality[:, k, 0], ownquality[:, k, 1], fractions[:, k, 0], fractions[:, k, 1]]) for k in range(K)], axis=1)
    for k, target in enumerate(selected.itertuples()):
        j = old['axis_ids'].tolist().index(target.axis_id)
        C, original_names = old_controls(old, j)
        assert names == original_names
        np.testing.assert_allclose(control[:129, k], C, atol=1e-12, rtol=1e-12)
    a = dict(record_ids=np.asarray(cases), partition=blind.partition.to_numpy(str), profile_blocks=blind.profile_block.to_numpy(str), cancer_type=blind.cancer_type.to_numpy(str), response=pd.to_numeric(clinical.response, errors='coerce').to_numpy(float), controls=control, control_names=np.asarray(names), base_control_count=np.asarray(base.shape[1]), scope_ids=selected.axis_id.to_numpy(str), references=np.asarray([d['reference'] for d in definitions]), axis_score=axis_score, function_score=function, function_eligible=eligible, endpoint_measurable=endpoint, unit_indices=unit_ids, unit_cells=ncell)
    np.savez_compressed(out / 'private_complete_function_inputs.npz', **a)
    np.savez_compressed(out / 'private_discovery_only_gene_and_function_scaling.npz', native_columns=np.asarray(columns), gene_mean=centers, gene_SD=scales, function_mean_SD=function_parameters)
    pd.DataFrame(rows).to_csv(out / 'private_all_reference_expression_coverage.tsv', sep='\t', index=False)
    coverage = []
    exclusions = []
    for slot, k, f, partition, question in slots(a):
        m, g, status, valid = prepare(a, k, f, partition, question)
        coverage.append(dict(slot=slot, scope_id=a['scope_ids'][k], reference=a['references'][f], partition=partition, question=question, status=status, author_records=int(valid.sum()), nonMPR=int((a['response'][valid] == 0).sum()), pCR_or_MPR=int((a['response'][valid] == 1).sum()), profile_blocks=len(set(a['profile_blocks'][valid])), residual_df=m['df'] if m else np.nan, max_block_leverage=g['max_block_leverage'] if g else np.nan))
    for k, scope in enumerate(a['scope_ids']):
        for i, record in enumerate(cases):
            why = []
            if not endpoint[i, k]: why.append('SENDER_OR_RECEIVER_FEWER_THAN20_OR_ABSENT')
            if not np.isfinite(np.r_[a['response'][i], control[i, k]]).all(): why.append('MISSING_FIXED_CLINICAL_OR_RNA_CONTROL')
            exclusions.append(dict(scope_id=scope, record_id=record, partition=a['partition'][i], reasons='|'.join(why) if why else 'ENDPOINT_AND_CONTROLS_AVAILABLE_FUNCTION_GATE_SEPARATE'))
    pd.DataFrame(coverage).to_csv(out / 'private_complete_family_design_coverage.tsv', sep='\t', index=False)
    pd.DataFrame(exclusions).to_csv(out / 'private_all_record_eligibility.tsv.gz', sep='\t', index=False, compression='gzip')
    files = {p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}
    assert sum(d['bytes'] for d in files.values()) < 4 * 1024**2
    receipt = dict(status='COMPLETE_ALL_REFERENCE_EXPRESSION_AND_DESIGN_INPUTS_NO_FUNCTION_EFFECTS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), reference_audit_sha256=sha(proof), files=files, clinical_function_effects_read=False, scientific_outputs_public=False)
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
