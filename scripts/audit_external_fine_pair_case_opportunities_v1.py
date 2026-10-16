"""Set-intersection replay of all fine pair/source/stratum/paired counts."""
import hashlib
import itertools
import json
import os
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'external_fine_pair_case_opportunities_20261016_v1'
    terminal = json.loads((out / 'private_terminal.json').read_text())
    for name, digest in terminal['files'].items():
        assert sha(out / name) == digest
    pre = root / 'existing_immune_cohort_candidate_extension_preflight_20261016_v2'
    lib = pd.read_csv(pre / 'private_all_library_pool_input_coverage.tsv', sep='\t', dtype=str, keep_default_na=False)
    libraries = {(r.dataset, r.library_id): r for r in lib.itertuples()}
    units = pd.read_csv(root / 'existing_immune_fine_native_measurement_scope_20261016_v1/private_all_frozen_T_myeloid_node_units.tsv', sep='\t', dtype={'library_id': str})
    sufficient = defaultdict(set)
    for r in units.itertuples():
        if r.n_cells >= 20:
            sufficient[(r.dataset, r.node_id)].add(r.library_id)
    run = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    z = zarr.open_group(str(run / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    text = read_elem(z['uns']['subcluster_catalog_json'])
    if isinstance(text, np.ndarray):
        text = text.item()
    catalog = json.loads(text)
    pools = {'Myeloid': {r['node_id'] for r in catalog if r['major'] == 'Myeloid'},
             'CD8_candidate': {r['node_id'] for r in catalog if r['subtype_candidate'] == 'CD8_Effector' and r['subtype_evidence_status'] == 'multigene_subtype_candidate'},
             'Treg_candidate': {r['node_id'] for r in catalog if r['subtype_candidate'] == 'Treg' and r['subtype_evidence_status'] == 'multigene_subtype_candidate'}}
    families = pd.read_csv(out / 'private_all_deduplicated_extension_families.tsv', sep='\t')
    expected_keys = {(s, f.extension_family, a, b) for s in set(lib.dataset) for f in families.itertuples()
                     for a, b in itertools.product(pools[f.sender_pool], pools[f.receiver_pool])}
    coverage = pd.read_csv(out / 'private_all_source_fine_pair_coverage.tsv', sep='\t')
    assert not coverage.duplicated(['dataset', 'extension_family', 'sender_node', 'receiver_node']).any()
    assert expected_keys == set(zip(coverage.dataset, coverage.extension_family, coverage.sender_node, coverage.receiver_node))
    fields = ['tissue', 'timepoint', 'treatment_arm', 'capture_scope', 'ex_vivo_stimulation']
    expected_strata = {}; paired = {}
    for r in coverage.itertuples():
        found = sufficient[(r.dataset, r.sender_node)] & sufficient[(r.dataset, r.receiver_node)]
        records = [libraries[(r.dataset, b)] for b in found]
        assert r.both20_libraries == len(found)
        assert r.reviewed_patient_blocks == len({a.patient_block for a in records if a.patient_identity_verified == 'True'})
        assert r.unreviewed_library_fallback == sum(a.patient_identity_verified != 'True' for a in records)
        assert r.non_product_non_exvivo_libraries == sum(a.CAR_product_or_ex_vivo_scope == 'False' for a in records)
        assert not r.exact_NSCLC_endpoints_validated
        by_level = defaultdict(list)
        for a in records:
            by_level[tuple(getattr(a, c) for c in fields)].append(a)
        for levels, group in by_level.items():
            key = (r.dataset, r.extension_family, r.sender_node, r.receiver_node, *levels)
            expected_strata[key] = (len(group), len({a.patient_block for a in group if a.patient_identity_verified == 'True'}), sum(a.patient_identity_verified != 'True' for a in group))
        if r.dataset == 'GSE286827':
            records = [a for a in records if a.patient_identity_verified == 'True' and a.CAR_product_or_ex_vivo_scope == 'False']
            pre_patients = {a.patient_block for a in records if a.timepoint.lower() == 'pre'}
            post_patients = {a.patient_block for a in records if a.timepoint.lower() == 'post'}
            repeated = Counter((a.patient_block, a.timepoint.lower()) for a in records)
            paired[(r.extension_family, r.sender_node, r.receiver_node)] = (
                len(pre_patients & post_patients), sum(a.timepoint.lower() == 'pre' for a in records),
                sum(a.timepoint.lower() == 'post' for a in records), sum(n > 1 for n in repeated.values()),
                sum(a.timepoint.lower() not in ['pre', 'post'] for a in records))
    observed = pd.read_csv(out / 'private_all_source_fine_pair_strata.tsv.gz', sep='\t', dtype={c: str for c in fields}, keep_default_na=False)
    assert len(observed) == len(expected_strata)
    for r in observed.itertuples():
        key = (r.dataset, r.extension_family, r.sender_node, r.receiver_node, *(getattr(r, c) for c in fields))
        assert expected_strata.pop(key) == (r.libraries, r.reviewed_patient_blocks, r.unreviewed_library_fallback)
    assert not expected_strata
    paired_table = pd.read_csv(out / 'private_all_HNSCC_paired_fine_opportunities.tsv', sep='\t')
    assert len(paired_table) == len(paired)
    for r in paired_table.itertuples():
        key = (r.extension_family, r.sender_node, r.receiver_node)
        assert paired.pop(key) == (r.reviewed_pre_post_patient_pairs, r.qualified_pre_libraries, r.qualified_post_libraries,
                                   r.repeated_patient_phase_groups, r.recurrence_or_other_libraries)
    assert not paired
    proof = dict(status='PASS_ALL_SOURCE_FINE_PAIR_STRATA_AND_REVIEWED_PRE_POST_COUNTS_BY_INDEPENDENT_SET_INTERSECTIONS',
                 finished_UTC=datetime.now(timezone.utc).isoformat(), audit_source_sha256=sha(__file__), producer_terminal_sha256=sha(out / 'private_terminal.json'),
                 source_pair_rows=len(coverage), stratum_rows=len(observed), paired_rows=len(paired_table), no_effects_or_expression=True)
    p = out / 'private_independent_audit.json'
    assert not p.exists()
    p.write_text(json.dumps(proof, indent=2) + '\n')
    print(json.dumps(proof))


if __name__ == '__main__':
    main()
