"""Freeze every prior-compatible TF slot for all original primary-screen rows.

This follow-up is selected in the same NSCLC data and is exploratory. Gene
disjointness supplies a different RNA readout, never independent samples.
This preflight computes no TF activity, association, P or revised LR q value.
"""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8388608), b''):
            h.update(block)
    return h.hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    results = root / 'NSCLC_automatic_ligand_response_results_20261016_v1'
    original_audit = json.loads((results / 'private_independent_audit.json').read_text())
    assert original_audit['status'].startswith('PASS_COMPLETE_UNIQUE_FAMILY_')
    original_path = results / 'private_all_unique_ligand_response_results.tsv'
    assert sha(original_path) == original_audit['result_sha256']
    all_results = pd.read_csv(original_path, sep='\t')
    selected = all_results[all_results.BH_q_complete_unique_family.lt(.05)].copy()
    assert len(selected) > 0
    reference = root / 'fullgene_receiver_target_code_v1'
    prior = reference / 'receiver_TF_prior_v1'
    prior_audit = json.loads((prior / 'independent_audit.json').read_text())
    assert prior_audit['status'].startswith('PASS')
    inventory = json.loads((prior / 'inventory.json').read_text())
    path_file = prior / 'private_all_receptor_TF_paths.tsv'
    assert sha(path_file) == inventory['outputs'][path_file.name]
    tf_file = reference / 'collectri_reference_v1/human_collectri_signed_complexes.tsv'
    official = reference / 'collectri_reference_v1/official_decoupler192_ULM_source.txt'
    params_file = root / 'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'
    params = pd.read_csv(params_file, sep='\t')
    assert not params.duplicated(['context', 'type', 'program', 'gene']).any()
    net = pd.read_csv(tf_file, sep='\t')
    assert not net.duplicated(['source', 'target']).any() and set(net.weight) <= {-1, 1}
    paths = pd.read_csv(path_file, sep='\t')
    assert not paths.duplicated(['receptor', 'TF']).any()
    paths = paths[~paths.sign_ambiguous & paths.expected_sign.isin([-1, 1]) & paths.min_edges.le(3)]
    tf_targets = {tf: set(f.target) for tf, f in net.groupby('source')}
    native_file = root / 'NSCLC_full_native_gene_units_20261016_v1/private_native_gene_axis.tsv'
    native = set(pd.read_csv(native_file, sep='\t').gene)
    original_input = root / 'NSCLC_automatic_ligand_response_preflight_20261016_v1'
    hashes = pd.read_csv(original_input / 'private_all_ready_design_SHA256.tsv', sep='\t').set_index('record').sha256
    out = root / 'NSCLC_gene_disjoint_receiver_TF_preflight_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__),
        original_results_sha256=sha(original_path), original_audit_sha256=sha(results / 'private_independent_audit.json'),
        signed_TF_network_sha256=sha(tf_file), receptor_TF_paths_sha256=sha(path_file), receptor_TF_audit_sha256=sha(prior / 'independent_audit.json'),
        official_ULM_source_sha256=sha(official), private_frozen_parameters_sha256=sha(params_file), native_gene_axis_sha256=sha(native_file),
        selection='ALL original primary BH records; all unambiguous signed prior receptor-to-TF paths <=3 edges, no top-N or expression/association selection.',
        gene_mask='Remove ALL frozen sender and receiver programme genes, literal ligand and ALL original qualified receptor components from BOTH native scoring background and TF targets. Do not export these private learned gene lists.',
        coverage='At least20 remaining native TF targets and500 remaining native background genes. Every failed slot retained as unavailable, never a negative.',
        prospective_scoring='Official decoupler1.9.2 ULM-equivalent signed TF score on native mean per-cell log1pCP10K across the full remaining native background, not the TF transcript; no TF activity computed in this preflight.',
        cases='EXACT original matched case/sender/receiver units, full original nuisance set unchanged; no deletion to increase effect.',
        multiplicity='Scope freeze only. No new P/BH or association model here. Any later inference must use the complete deduplicated parent-record/TF family with unavailable hypotheses retained, and be labeled conditional same-source exploratory follow-up.',
        receptor_ambiguity='Several qualified receptors can reach one TF with opposing signs. Preserve ALL receptor path signs and mark ambiguous union, never select the sign matching an observed effect.',
        molecular_boundary='Native/cofactor diagnostics and literal binding caveats do not disappear. A TF response cannot repair a biochemically incomplete original ligand/receptor identity.',
        original_LR_family_changed=False, independent_samples_or_protein_validation=False)
    (out / 'plan_before_any_TF_activity_or_association.json').write_text(json.dumps(plan, indent=2) + '\n')
    rows, references, cases = [], [], []
    for row in selected.itertuples():
        sender_type = 'typeB' if row.receiver_type == 'typeA' else 'typeA'
        recv = params[params.context.eq(row.context) & params.type.eq(row.receiver_type) & params.program.eq(row.program)]
        send = params[params.context.eq(row.context) & params.type.eq(sender_type) & params.program.eq(row.program)]
        assert len(recv) and len(send)
        receptors = row.qualified_receptors.split('|')
        forbidden = set(recv.gene) | set(send.gene) | {row.ligand} | {g for r in receptors for g in r.split('_')}
        background = native - forbidden
        design = original_input / 'private_frozen_designs' / (row.record + '.npz')
        assert sha(design) == hashes.loc[row.record]
        z = np.load(design, allow_pickle=False)
        assert z['Z'].shape == (int(row.cases), 19)
        cases.append(dict(record=row.record, cases=len(z['cases']), original_design_sha256=sha(design), disjoint_background_genes=len(background), all_frozen_programme_genes_excluded=True))
        candidates = paths[paths.receptor.isin(receptors) & paths.TF.isin(tf_targets)]
        if candidates.empty:
            cases[-1]['scope_status'] = 'NO_UNAMBIGUOUS_REFERENCE_PATH_NOT_BIOLOGICAL_NEGATIVE'
        else:
            cases[-1]['scope_status'] = 'ALL_REFERENCE_PATHS_ACCOUNTED'
        for tf, f in candidates.groupby('TF', sort=True):
            targets = tf_targets[tf] & background
            status = 'READY_FOR_INPUT_AUDIT_NO_ACTIVITY_COMPUTED' if len(targets) >= 20 and len(background) >= 500 else 'INSUFFICIENT_DISJOINT_REFERENCE_TARGETS_OR_BACKGROUND'
            signs = set(f.expected_sign)
            rows.append(dict(record=row.record, TF=tf, context=row.context, program=row.program, sender_author_subtype=row.sender_author_subtype, receiver_author_subtype=row.receiver_author_subtype,
                native_disjoint_targets=len(targets), native_disjoint_background=len(background), qualified_reference_receptor_paths=len(f), sign_ambiguous_across_receptors=len(signs) > 1,
                union_expected_sign=next(iter(signs)) if len(signs) == 1 else 0, status=status, cases=int(row.cases), original_design_sha256=sha(design), same_source_postselection_exploratory=True))
            for r in f.itertuples():
                references.append(dict(record=row.record, TF=tf, receptor=r.receptor, expected_sign=int(r.expected_sign), shortest_edges=int(r.min_edges)))
    table = pd.DataFrame(rows)
    assert len(table) and not table.duplicated(['record', 'TF']).any()
    pd.DataFrame(cases).to_csv(out / 'private_all_parent_case_scope.tsv', sep='\t', index=False)
    table.to_csv(out / 'private_all_deduplicated_TF_slots.tsv', sep='\t', index=False)
    pd.DataFrame(references).to_csv(out / 'private_all_TF_receptor_sign_provenance.tsv', sep='\t', index=False)
    names = ['private_all_parent_case_scope.tsv', 'private_all_deduplicated_TF_slots.tsv', 'private_all_TF_receptor_sign_provenance.tsv', 'plan_before_any_TF_activity_or_association.json']
    terminal = dict(status='COMPLETE_ALL_PRIMARY_PARENT_PRIOR_TF_SCOPE_NO_ACTIVITY_OR_ASSOCIATION_INPUT_AUDIT_PENDING', finished_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__), parents=len(selected), unique_parent_TF_slots=len(table), coverage_status_counts=table.status.value_counts().to_dict(), ambiguous_receptor_union_slots=int(table.sign_ambiguous_across_receptors.sum()), files={name: sha(out / name) for name in names}, TF_scores_or_new_P_or_BH_computed=False, learned_programme_gene_lists_exported=False, independent_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps(terminal))


if __name__ == '__main__':
    main()
