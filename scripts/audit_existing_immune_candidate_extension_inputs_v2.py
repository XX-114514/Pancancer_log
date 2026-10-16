"""V2: preserve source metadata strings when parsing CSV audit inputs.

Only parser typing changes; all original scientific gates and inputs retained.
"""
import hashlib
import json
import os
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'existing_immune_cohort_candidate_extension_preflight_20261016_v1'
    terminal = json.loads((out / 'private_terminal.json').read_text())
    for name, digest in terminal['files'].items():
        assert sha(out / name) == digest
    dest = out / 'private_independent_audit.json'
    assert not dest.exists()
    plan = json.loads((out / 'plan_before_external_input_review.json').read_text())
    run = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    assert sha(run / 'manifests/immune_authority_release_public.json') == plan['atlas_release_sha256']
    z = zarr.open_group(str(run / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    source_order = np.asarray(read_elem(z['uns']['source_order'])).astype(str).tolist()
    assert source_order == plan['sources']
    fields = ['dataset', 'library_id', 'sampleID', 'node_id', 'patient_block', 'patient_identity_verified', 'timepoint', 'tissue', 'treatment_arm', 'capture_scope', 'ex_vivo_stimulation']
    data = {k: np.asarray(read_elem(z['obs'][k])).astype(str) for k in fields}
    catalog_text = read_elem(z['uns']['subcluster_catalog_json'])
    if isinstance(catalog_text, np.ndarray):
        catalog_text = catalog_text.item()
    assert hashlib.sha256(str(catalog_text).encode()).hexdigest() == plan['pool_catalog_sha256']
    catalog = json.loads(catalog_text)
    pools = dict(Myeloid={r['node_id'] for r in catalog if r['major'] == 'Myeloid'}, CD8_candidate={r['node_id'] for r in catalog if r['major'] == 'T' and r['subtype_candidate'] == 'CD8_Effector' and r['subtype_evidence_status'] == 'multigene_subtype_candidate'}, Treg_candidate={r['node_id'] for r in catalog if r['major'] == 'T' and r['subtype_candidate'] == 'Treg' and r['subtype_evidence_status'] == 'multigene_subtype_candidate'})
    assert {k: sorted(v) for k, v in pools.items()} == plan['diagnostic_pool_node_ids']
    libraries = pd.read_csv(out / 'private_all_library_pool_input_coverage.tsv', sep='\t', keep_default_na=False, dtype={name: str for name in fields if name not in ['node_id']}).set_index(['dataset', 'library_id'])
    grouped = pd.DataFrame({'dataset': data['dataset'], 'library_id': data['library_id']}).groupby(['dataset', 'library_id'], sort=True).indices
    assert set(grouped) == set(libraries.index) and len(libraries) == 400
    for key, indices in grouped.items():
        row = libraries.loc[key]
        assert row.all_cells == len(indices)
        for name in [f for f in fields if f not in ['dataset', 'library_id', 'node_id']]:
            values = np.unique(data[name][indices])
            expected = str(values[0]) if len(values) == 1 else 'CONFLICTING_LIBRARY_METADATA'
            assert row[name] == expected
        for name, ids in pools.items():
            assert row[name + '_cells'] == np.isin(data['node_id'][indices], list(ids)).sum()
    assert libraries.all_cells.sum() == 2670057
    gene_symbols = np.asarray(read_elem(z['var']['gene_symbol'])).astype(str)
    aliases = np.asarray(read_elem(z['var']['source_gene_aliases'])).astype(str)
    mask = np.asarray(read_elem(z['varm']['measured_in_source']), bool)
    genes = pd.read_csv(out / 'private_all_source_gene_measurement_status.tsv.gz', sep='\t').set_index(['dataset', 'gene'])
    assert genes.index.is_unique
    lookup = {}
    for source_index, source in enumerate(source_order):
        mapping = defaultdict(set)
        for column in np.flatnonzero(mask[:, source_index]):
            mapping[gene_symbols[column]].add(int(column))
            mapping[aliases[column]].add(int(column))
        lookup[source] = mapping
        for (sid, gene), row in genes.loc[[source]].iterrows():
            count = len(mapping[gene])
            assert row.exact_source_physical_matches == count
            assert row.status == ('UNAMBIGUOUS_SOURCE_NATIVE_COLUMN' if count == 1 else 'SYMBOL_OR_ALIAS_REVIEW_REQUIRED_NOT_BIOLOGICAL_NEGATIVE')
    parent_path = root / 'NSCLC_automatic_ligand_response_review_20261016_v1/private_all_primary_BH_supported_rows.tsv'
    params_path = root / 'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'
    assert sha(parent_path) == plan['original_evidence_sha256'] and sha(params_path) == plan['private_frozen_parameters_sha256']
    parents = pd.read_csv(parent_path, sep='\t')
    params = pd.read_csv(params_path, sep='\t')
    counts = pd.read_csv(out / 'private_all_source_parent_coverage.tsv', sep='\t').set_index(['dataset', 'record'])
    assert counts.index.is_unique and set(counts.index) == {(s, r) for s in source_order for r in parents.record}
    expected_genes = set(params.gene) | set(parents.ligand)
    expected_genes.update(g for receptor in parents.qualified_receptors for r in receptor.split('|') for g in r.split('_'))
    assert set(genes.index) == {(s, g) for s in source_order for g in expected_genes}
    lib = libraries.reset_index()
    strata = []
    for source in source_order:
        mapping = lookup[source]
        for parent in parents.itertuples():
            sender_pool = 'Myeloid' if parent.sender_author_subtype.startswith('M\u03c6_') else 'Treg_candidate'
            receiver_pool = 'Myeloid' if parent.receiver_author_subtype.startswith('M\u03c6_') else 'CD8_candidate'
            if parent.receiver_author_subtype.startswith('CD4T_Treg_'):
                receiver_pool = 'Treg_candidate'
            eligible = lib[lib.dataset.eq(source) & lib[sender_pool + '_cells'].ge(20) & lib[receiver_pool + '_cells'].ge(20)]
            observed = counts.loc[(source, parent.record)]
            assert observed.both_pool20_libraries == len(eligible)
            assert observed.reviewed_patient_blocks == eligible.loc[eligible.patient_identity_verified.eq('True'), 'patient_block'].nunique()
            assert observed.unreviewed_library_units == (~eligible.patient_identity_verified.eq('True')).sum()
            assert observed.CAR_or_ex_vivo_libraries == eligible.CAR_product_or_ex_vivo_scope.sum()
            assert observed.ligand_unambiguous_native == (len(mapping[parent.ligand]) == 1)
            for role, typ in [('receiver', parent.receiver_type), ('sender', 'typeB' if parent.receiver_type == 'typeA' else 'typeA')]:
                p = params[params.context.eq(parent.context) & params.type.eq(typ) & params.program.eq(parent.program)]
                indices = [next(iter(mapping[g])) if len(mapping[g]) == 1 else None for g in p.gene]
                cc = Counter(indices)
                total = sum(abs(w) for w in p.weight)
                measured = sum(abs(w) for w, column in zip(p.weight, indices) if column is not None and cc[column] == 1)
                fraction = measured / total
                assert np.isclose(observed[role + '_programme_unambiguous_native_L1_fraction'], fraction, rtol=1e-12, atol=1e-12)
                assert observed[role + '_programme95pct_technical_gate'] == (fraction >= .95)
            assert observed.all_original_receptors_native == all(len(mapping[g]) == 1 for r in parent.qualified_receptors.split('|') for g in r.split('_'))
            keys = ['tissue', 'timepoint', 'treatment_arm', 'capture_scope', 'ex_vivo_stimulation']
            for levels, gg in eligible.groupby(keys, sort=True, dropna=False):
                strata.append(dict(dataset=source, record=parent.record, **dict(zip(keys, levels)), libraries=len(gg), reviewed_patient_blocks=int(gg.loc[gg.patient_identity_verified.eq('True'), 'patient_block'].nunique()), unreviewed_library_units=int((~gg.patient_identity_verified.eq('True')).sum()), unknown_or_conflicting_stratum=any(v in ['', 'unknown', 'Unknown', 'CONFLICTING_LIBRARY_METADATA'] for v in levels), no_effects=True))
    observed = pd.read_csv(out / 'private_all_source_parent_strata.tsv', sep='\t', keep_default_na=False, dtype={name: str for name in ['dataset', 'record', 'tissue', 'timepoint', 'treatment_arm', 'capture_scope', 'ex_vivo_stimulation']})
    pd.testing.assert_frame_equal(observed, pd.DataFrame(strata), check_dtype=False)
    result = dict(status='PASS_ALL_ATLAS_LIBRARY_METADATA_POOL_COUNTS_SOURCE_PHYSICAL_GENE_MATCHES_PROGRAMME_COVERAGE_AND_ALL_PARENT_STRATA', finished_UTC=datetime.now(timezone.utc).isoformat(), audit_source_sha256=sha(__file__), producer_terminal_sha256=sha(out / 'private_terminal.json'), sources=9, libraries=len(libraries), source_parent_rows=len(counts), strata_rows=len(strata), gene_mask_rows=len(genes), no_expression_or_clinical_effects=True, exact_endpoint_mapping_or_independent_validation_accepted=False)
    dest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
