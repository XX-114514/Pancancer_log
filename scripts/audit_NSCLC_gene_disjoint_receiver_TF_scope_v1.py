"""Independent accounting of all parent/TF input slots; no scientific effects."""
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'NSCLC_gene_disjoint_receiver_TF_preflight_20261016_v1'
    job = root / 'HNSCC_receiver_followup_job_control_20261016_v2/preflight_NSCLC_gene_disjoint_receiver_TF_v1/launch_receipt.json'
    launch = json.loads(job.read_text())
    deadline = time.monotonic() + 3600
    while not (out / 'private_terminal.json').exists():
        proc = Path('/proc') / str(launch['pid']) / 'stat'
        assert proc.exists()
        state = proc.read_text().rsplit(')', 1)[1].split()
        assert state[19] == str(launch['proc_start_ticks']) and state[0] != 'Z'
        assert time.monotonic() < deadline
        time.sleep(10)
    dest = out / 'private_independent_audit.json'
    assert not dest.exists()
    terminal = json.loads((out / 'private_terminal.json').read_text())
    for name, digest in terminal['files'].items():
        assert sha(out / name) == digest
    plan = json.loads((out / 'plan_before_any_TF_activity_or_association.json').read_text())
    assert plan['producer_sha256'] == launch['source_sha256']
    rdir = root / 'NSCLC_automatic_ligand_response_results_20261016_v1'
    assert sha(rdir / 'private_independent_audit.json') == plan['original_audit_sha256']
    assert sha(rdir / 'private_all_unique_ligand_response_results.tsv') == plan['original_results_sha256']
    parents = pd.read_csv(rdir / 'private_all_unique_ligand_response_results.tsv', sep='\t')
    parents = parents[parents.BH_q_complete_unique_family.lt(.05)]
    c = root / 'fullgene_receiver_target_code_v1'
    pairs = [(c / 'collectri_reference_v1/human_collectri_signed_complexes.tsv', 'signed_TF_network_sha256'), (c / 'receiver_TF_prior_v1/private_all_receptor_TF_paths.tsv', 'receptor_TF_paths_sha256'), (c / 'receiver_TF_prior_v1/independent_audit.json', 'receptor_TF_audit_sha256'), (c / 'collectri_reference_v1/official_decoupler192_ULM_source.txt', 'official_ULM_source_sha256'), (root / 'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv', 'private_frozen_parameters_sha256'), (root / 'NSCLC_full_native_gene_units_20261016_v1/private_native_gene_axis.tsv', 'native_gene_axis_sha256')]
    for path, key in pairs:
        assert sha(path) == plan[key]
    net = pd.read_csv(pairs[0][0], sep='\t')
    paths = pd.read_csv(pairs[1][0], sep='\t')
    params = pd.read_csv(pairs[4][0], sep='\t')
    native = set(pd.read_csv(pairs[5][0], sep='\t').gene)
    p = paths.query('sign_ambiguous == False and min_edges <= 3 and expected_sign != 0')
    slots = pd.read_csv(out / 'private_all_deduplicated_TF_slots.tsv', sep='\t')
    refs = pd.read_csv(out / 'private_all_TF_receptor_sign_provenance.tsv', sep='\t')
    scopes = pd.read_csv(out / 'private_all_parent_case_scope.tsv', sep='\t').set_index('record')
    assert scopes.index.tolist() == parents.record.tolist()
    assert not slots.duplicated(['record', 'TF']).any()
    assert not refs.duplicated(['record', 'TF', 'receptor']).any()
    expected_keys, expected_refs = set(), set()
    for row in parents.itertuples():
        parameters = params[params.context.eq(row.context) & params.program.eq(row.program)]
        assert set(parameters.type) == {'typeA', 'typeB'}
        receptors = row.qualified_receptors.split('|')
        forbidden = set(parameters.gene) | {row.ligand}
        forbidden.update(g for receptor in receptors for g in receptor.split('_'))
        background = native.difference(forbidden)
        design = root / 'NSCLC_automatic_ligand_response_preflight_20261016_v1/private_frozen_designs' / (row.record + '.npz')
        z = np.load(design, allow_pickle=False)
        scope = scopes.loc[row.record]
        assert scope.cases == len(z['cases']) == int(row.cases)
        assert scope.disjoint_background_genes == len(background)
        assert scope.all_frozen_programme_genes_excluded and scope.original_design_sha256 == sha(design)
        eligible = p[p.receptor.isin(receptors) & p.TF.isin(net.source)]
        expected_status = 'ALL_REFERENCE_PATHS_ACCOUNTED' if len(eligible) else 'NO_UNAMBIGUOUS_REFERENCE_PATH_NOT_BIOLOGICAL_NEGATIVE'
        assert scope.scope_status == expected_status
        for tf in sorted(eligible.TF.unique()):
            expected_keys.add((row.record, tf))
            f = eligible[eligible.TF.eq(tf)]
            targets = net[net.source.eq(tf)].target
            n_targets = int(targets.isin(background).sum())
            observed = slots[slots.record.eq(row.record) & slots.TF.eq(tf)].iloc[0]
            assert observed.native_disjoint_targets == n_targets and observed.native_disjoint_background == len(background)
            assert observed.qualified_reference_receptor_paths == len(f)
            signs = sorted(f.expected_sign.unique())
            assert observed.sign_ambiguous_across_receptors == (len(signs) > 1)
            assert observed.union_expected_sign == (signs[0] if len(signs) == 1 else 0)
            status = 'READY_FOR_INPUT_AUDIT_NO_ACTIVITY_COMPUTED' if n_targets >= 20 and len(background) >= 500 else 'INSUFFICIENT_DISJOINT_REFERENCE_TARGETS_OR_BACKGROUND'
            assert observed.status == status and observed.cases == len(z['cases'])
            assert observed.context == row.context and observed.program == row.program
            assert observed.sender_author_subtype == row.sender_author_subtype and observed.receiver_author_subtype == row.receiver_author_subtype
            assert observed.original_design_sha256 == sha(design) and observed.same_source_postselection_exploratory
            expected_refs.update((row.record, tf, r.receptor, int(r.expected_sign), int(r.min_edges)) for r in f.itertuples())
    assert expected_keys == set(slots[['record', 'TF']].itertuples(index=False, name=None))
    assert expected_refs == set(refs[['record', 'TF', 'receptor', 'expected_sign', 'shortest_edges']].itertuples(index=False, name=None))
    assert terminal['parents'] == len(parents) and terminal['unique_parent_TF_slots'] == len(slots)
    assert terminal['coverage_status_counts'] == slots.status.value_counts().to_dict()
    result = dict(status='PASS_ALL_ORIGINAL_PRIMARY_PARENTS_ALL_REFERENCE_PATHS_DEDUPLICATED_TF_MASKS_COVERAGE_AND_MATCHED_CASE_DESIGNS', finished_UTC=datetime.now(timezone.utc).isoformat(), audit_code_sha256=sha(__file__), producer_terminal_sha256=sha(out / 'private_terminal.json'), parents=len(parents), all_TF_slots=len(slots), all_receptor_path_records=len(refs), no_TF_scores_associations_P_or_BH_computed=True, same_source_postselection_exploratory=True, learned_programme_gene_lists_exported=False)
    dest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
