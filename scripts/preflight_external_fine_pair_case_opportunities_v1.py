"""Effect-blind fine endpoint and paired-case opportunities, no expression.

Broader original atlas subtype candidates are extensions, never exact NSCLC
author-state mappings. All candidate node pairs retained, without P/results.
"""
import itertools
import hashlib
import json
import os
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
    pre = root / 'existing_immune_cohort_candidate_extension_preflight_20261016_v2'
    assert json.loads((pre / 'private_independent_audit.json').read_text())['status'].startswith('PASS_')
    units_file = root / 'existing_immune_fine_native_measurement_scope_20261016_v1/private_all_frozen_T_myeloid_node_units.tsv'
    units = pd.read_csv(units_file, sep='\t', dtype={'library_id': str, 'dataset': str, 'node_id': str})
    libraries = pd.read_csv(pre / 'private_all_library_pool_input_coverage.tsv', sep='\t', dtype=str, keep_default_na=False)
    assert not libraries.duplicated(['dataset', 'library_id']).any()
    run = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    z = zarr.open_group(str(run / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    text = read_elem(z['uns']['subcluster_catalog_json'])
    if isinstance(text, np.ndarray):
        text = text.item()
    catalog = pd.DataFrame(json.loads(text))
    pools = {'Myeloid': catalog.loc[catalog.major.eq('Myeloid'), 'node_id'].tolist(),
             'CD8_candidate': catalog.loc[catalog.subtype_candidate.eq('CD8_Effector') & catalog.subtype_evidence_status.eq('multigene_subtype_candidate'), 'node_id'].tolist(),
             'Treg_candidate': catalog.loc[catalog.subtype_candidate.eq('Treg') & catalog.subtype_evidence_status.eq('multigene_subtype_candidate'), 'node_id'].tolist()}
    parents_file = root / 'NSCLC_automatic_ligand_response_review_20261016_v1/private_all_primary_BH_supported_rows.tsv'
    parents = pd.read_csv(parents_file, sep='\t')
    parents['sender_pool'] = np.where(parents.sender_author_subtype.str.startswith('Mφ_'), 'Myeloid', 'Treg_candidate')
    parents['receiver_pool'] = np.where(parents.receiver_author_subtype.str.startswith('Mφ_'), 'Myeloid', 'CD8_candidate')
    parents.loc[parents.receiver_author_subtype.str.startswith('CD4T_Treg_'), 'receiver_pool'] = 'Treg_candidate'
    keys = ['context', 'program', 'ligand', 'qualified_receptors', 'sender_pool', 'receiver_pool']
    families = parents.groupby(keys, sort=True).record.agg(lambda v: '|'.join(sorted(v))).reset_index(name='original_parent_records')
    families['extension_family'] = [f'EF{i:02d}' for i in range(1, len(families) + 1)]
    out = root / 'external_fine_pair_case_opportunities_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__),
                original_units_sha256=sha(units_file), parent_records_sha256=sha(parents_file), input_audit_sha256=sha(pre / 'private_independent_audit.json'),
                scope='All original molecular-family diagnostic lineage extensions and ALL individual original candidate nodes. Deduplicate original parents differing only in author fine endpoint; never count them as independent external evidence.',
                case_gate='>=20 cells at EACH fine node in same physical library; report all zero-coverage source/pair rows, reviewed patient blocks separately from unreviewed library fallback.',
                strata='Dataset x original tissue x timepoint x treatment x capture x ex_vivo, no merging tissue categories or guessed blank times.',
                paired='HNSCC only explicit pre/post original timepoints, same reviewed patient; recurrence/post2 separate. Require both fine nodes >=20 at both times. All libraries retained; report repeat libraries rather than assume one.',
                hierarchy='Diagnostic pool counts are opportunities only. No pooled TF computed by averaging nonlinear fine TF scores; no node relabel or effect-selected merge.',
                limits='Coverage opportunity is not full design identifiability, exact NSCLC endpoint acceptance, independent patients for unknown IDs, or external replication. No expression/effects/clinical endpoint read.')
    (out / 'plan_before_any_external_effects.json').write_text(json.dumps(plan, indent=2) + '\n')
    families.to_csv(out / 'private_all_deduplicated_extension_families.tsv', sep='\t', index=False)
    counts = units.pivot_table(index=['dataset', 'library_id'], columns='node_id', values='n_cells', aggfunc='sum', fill_value=0)
    counts = counts.reindex(pd.MultiIndex.from_frame(libraries[['dataset', 'library_id']]), fill_value=0)
    counts = counts.reindex(columns=catalog.node_id, fill_value=0)
    assert counts.to_numpy().sum() == units.n_cells.sum()
    rows, strata, paired = [], [], []
    split = ['tissue', 'timepoint', 'treatment_arm', 'capture_scope', 'ex_vivo_stimulation']
    for family in families.itertuples():
        for sender, receiver in itertools.product(pools[family.sender_pool], pools[family.receiver_pool]):
            ok = (counts[sender].to_numpy() >= 20) & (counts[receiver].to_numpy() >= 20)
            for source, indices in libraries.groupby('dataset', sort=True).indices.items():
                g = libraries.iloc[indices].loc[ok[indices]]
                verified = g.patient_identity_verified.eq('True')
                rows.append(dict(extension_family=family.extension_family, dataset=source, sender_node=sender, receiver_node=receiver,
                                 both20_libraries=len(g), reviewed_patient_blocks=g.loc[verified, 'patient_block'].nunique(),
                                 unreviewed_library_fallback=int((~verified).sum()),
                                 non_product_non_exvivo_libraries=int(g.CAR_product_or_ex_vivo_scope.eq('False').sum()),
                                 exact_NSCLC_endpoints_validated=False))
                for levels, gg in g.groupby(split, sort=True, dropna=False):
                    vv = gg.patient_identity_verified.eq('True')
                    strata.append(dict(extension_family=family.extension_family, dataset=source, sender_node=sender, receiver_node=receiver,
                                       **dict(zip(split, levels)), libraries=len(gg), reviewed_patient_blocks=gg.loc[vv, 'patient_block'].nunique(),
                                       unreviewed_library_fallback=int((~vv).sum())))
                if source == 'GSE286827':
                    hh = g[g.patient_identity_verified.eq('True') & g.CAR_product_or_ex_vivo_scope.eq('False')].copy()
                    hh['phase'] = hh.timepoint.str.lower()
                    phases = hh[hh.phase.isin(['pre', 'post'])].groupby('patient_block').phase.agg(set)
                    both = phases.map(lambda v: {'pre', 'post'} <= v)
                    repeated = hh.groupby(['patient_block', 'phase']).size()
                    paired.append(dict(extension_family=family.extension_family, sender_node=sender, receiver_node=receiver,
                                       reviewed_pre_post_patient_pairs=int(both.sum()), qualified_pre_libraries=int(hh.phase.eq('pre').sum()),
                                       qualified_post_libraries=int(hh.phase.eq('post').sum()),
                                       repeated_patient_phase_groups=int(repeated.gt(1).sum()), recurrence_or_other_libraries=int((~hh.phase.isin(['pre', 'post'])).sum())))
    pd.DataFrame(rows).to_csv(out / 'private_all_source_fine_pair_coverage.tsv', sep='\t', index=False)
    pd.DataFrame(strata).to_csv(out / 'private_all_source_fine_pair_strata.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(paired).to_csv(out / 'private_all_HNSCC_paired_fine_opportunities.tsv', sep='\t', index=False)
    result = dict(status='COMPLETE_ALL_FROZEN_FINE_PAIR_AND_PAIRED_REVIEWED_CASE_OPPORTUNITIES_NO_EFFECTS',
                  finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), extension_families=len(families),
                  source_fine_pair_rows=len(rows), stratum_rows=len(strata), HNSCC_paired_rows=len(paired),
                  files={p.name: sha(p) for p in out.iterdir() if p.is_file()}, no_expression_or_effects_or_endpoint_acceptance=True)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
