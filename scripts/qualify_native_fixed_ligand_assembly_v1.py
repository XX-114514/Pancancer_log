"""Previously fixed assembly extension: cached marginals, no clinical tests."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 256 * 1024)
    request = root / 'private_complementary_spatial_interface_scope_request_20261019_v1.json'
    assert sha(request) == 'f31d71f1efdc8450f346020733b0edea354a962eeec57bcc58328ba36bd890d0'
    fixed = json.loads(request.read_text())
    raw_audit = root / 'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json'
    assert sha(raw_audit) == 'bb20365ee83758a7df8478ceadede2596fe749928dab5765e3bc0866d72dd8c1'
    inp = root / 'NSCLC_frozen92_axis_source_inputs_20261019_v2'
    it = json.loads((inp / 'private_terminal.json').read_text())
    for name, info in it['files'].items():
        assert sha(inp / name) == info['sha256']
    a = dict(np.load(inp / 'private_all129_all92_native_axis_and_design_inputs.npz', allow_pickle=False))
    prior = root / 'NSCLC_frozen92_axis_source_consistency_20261019_v1'
    audit = root / 'NSCLC_frozen92_axis_source_result_audit_20261019_v1/private_terminal.json'
    assert sha(audit) == '79e8aa29aec12db05d437f6e4ebc24dec7e97772bf62ca594434efe36522b9b7'
    assert json.loads(audit.read_text())['producer_terminal_sha256'] == sha(prior / 'private_terminal.json')
    pt = json.loads((prior / 'private_terminal.json').read_text())
    result_path = prior / 'private_all92_source_consistency_results.tsv'
    assert sha(result_path) == pt['files'][result_path.name]['sha256']
    selected = pd.read_csv(result_path, sep='\t').query('BH92 < .05').axis_id.tolist()
    assert set(selected) == {t['axis_id'] for t in fixed['targets']}
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    nt = json.loads((native / 'private_terminal.json').read_text())
    used = ['private_all_native_units_and_quality.tsv', 'private_native_gene_axis.tsv', 'positive_cell_counts.npy', 'mean_logCP10K.npy']
    for name in used:
        assert sha(native / name) == nt['outputs'][name]['sha256']
    u = pd.read_csv(native / used[0], sep='\t', keep_default_na=False)
    genes = pd.read_csv(native / used[1], sep='\t', keep_default_na=False).gene.to_numpy(str)
    assert len(genes) == len(set(genes)) == 31831
    positive = np.load(native / used[2], mmap_mode='r')
    means = np.load(native / used[3], mmap_mode='r')
    cases = sorted(u.sampleID.unique())
    assert len(cases) == 243
    lookup = {(r.sampleID, r.sub_cell_type): int(r.unit_index) for r in u.itertuples()}
    out = root / 'native_fixed_ligand_assembly_qualification_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_cached_assembly_marginals.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), private_fixed_request_sha256=sha(request), raw_cache_audit_sha256=sha(raw_audit), current_inputs_sha256={name: nt['outputs'][name]['sha256'] for name in used}, author_records=243, original_separate_label_records=129, scope='All original supported axes and243 source records; original129 and other114 described separately. Fixed assembly definition predates this cache query. No new response or association.', same_cell_bounds='For k required components, max(0,sum(positive)-(k-1)*n) <= same-cell joint positive <= min(positive). No independence assumption or protein assembly conclusion.', cache_raw_scope='Prior full cache producer and accepted audit; independent raw audit covered prespecified64 units/7957 cells, not every new assembly position.', full_expression_or_score_arrays_persisted=False, resource_guard=guard))
    records, summaries = [], []
    for target in fixed['targets']:
        axis = target['axis_id']
        j = a['axis_ids'].tolist().index(axis)
        sender = str(a['author_states'][int(a['sender_index'][j])])
        receiver = str(a['author_states'][int(a['receiver_index'][j])])
        old_ligand = target['original_ligand_genes']
        strict = list(dict.fromkeys(old_ligand + target['strict_ligand_assembly_addition']))
        receptor = target['original_receptor_genes']
        index = {g: int(np.where(genes == g)[0][0]) for g in strict + receptor}
        assert len(index) == len(set(index.values()))
        local = []
        for case in cases:
            si, ri = lookup.get((case, sender)), lookup.get((case, receiver))
            ns, nr = (int(u.iloc[si].n_cells) if si is not None else 0), (int(u.iloc[ri].n_cells) if ri is not None else 0)
            measured = ns >= 20 and nr >= 20
            lower, upper = 0, 0
            ratio = np.nan
            if ns:
                counts = positive[si, [index[g] for g in strict]].astype(np.int64)
                assert (counts >= 0).all() and (counts <= ns).all()
                lower = max(0, int(counts.sum()) - (len(strict) - 1) * ns)
                upper = int(counts.min())
                assert 0 <= lower <= upper <= ns
                original = min(float(means[si, index[g]]) for g in old_ligand)
                strict_mean = min(float(means[si, index[g]]) for g in strict)
                assert 0 <= strict_mean <= original + 1e-12
                if original > 0:
                    ratio = strict_mean / original
            local.append(dict(source='GSE243013', axis_id=axis, author_record_id=case, original129=case in set(a['author_record_ids']), sender_cells=ns, receiver_cells=nr, both_endpoints_atleast20=measured, guaranteed_sender_joint_positive_lower=lower, possible_sender_joint_positive_upper=upper, lower_fraction=lower / ns if ns else np.nan, upper_fraction=upper / ns if ns else np.nan))
            if measured:
                new_score = np.log1p(float(np.float32(min(float(means[si, index[g]]) for g in old_ligand))) * float(np.float32(min(float(means[ri, index[g]]) for g in receptor)))).astype(np.float32).astype(float)
                if case in set(a['author_record_ids']):
                    ii = a['author_record_ids'].tolist().index(case)
                    np.testing.assert_equal(new_score, a['raw_axis_score'][ii, j])
            local[-1]['_ratio_RAM_only'] = ratio
        data = pd.DataFrame(local)
        for partition, take in [('ALL243', np.ones(len(data), bool)), ('ORIGINAL129', data.original129.to_numpy()), ('OTHER114', ~data.original129.to_numpy())]:
            subset = data.loc[take & data.both_endpoints_atleast20]
            ratios = subset._ratio_RAM_only.dropna()
            summaries.append(dict(source='GSE243013', axis_id=axis, partition=partition, author_records=int(take.sum()), both_endpoints_atleast20=len(subset), guaranteed_sender_joint_positive_atleast1_records=int(subset.guaranteed_sender_joint_positive_lower.ge(1).sum()), lower_bound_atleast10percent_sender_records=int(subset.lower_fraction.ge(.1).sum()), median_joint_positive_lower_fraction=float(subset.lower_fraction.median()), median_joint_positive_upper_fraction=float(subset.upper_fraction.median()), positive_original_ligand_mean_records=len(ratios), median_strict_to_original_ligand_mean_ratio=float(ratios.median()), strict_ligand_mean_lower_than_original_records=int(ratios.lt(1 - 1e-10).sum()), reviewed_independent_patients=0, inferential_test=False, actual_joint_cells_or_protein_observed=False))
        records.extend(data.drop(columns='_ratio_RAM_only').to_dict(orient='records'))
    pd.DataFrame(records).to_csv(out / 'private_all243_sender_assembly_marginal_bounds.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(summaries).to_csv(out / 'private_complete3_source_partition_assembly_descriptions.tsv', sep='\t', index=False)
    receipt = dict(status='COMPLETE_ALL243_FIXED_ASSEMBLY_MARGINAL_BOUNDS_AND3_SOURCE_DESCRIPTIONS_NO_NEW_CLINICAL_EFFECT', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), private_fixed_request_sha256=sha(request), cached_native_audit_sha256=sha(raw_audit), author_record_positions=len(records), source_partition_positions=len(summaries), original_proxy_reconstructed_on_all129_measurable_records=True, actual_joint_positive_cells_or_protein_verified=False, full_expression_or_score_arrays_persisted=False, new_raw_data_read_or_downloaded=False, new_clinical_association=False, scientific_outputs_public=False, resource_guard=guard, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in receipt['files'].values()) < 256 * 1024
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
