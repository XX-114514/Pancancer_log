"""V2 stable-ID aliases on all nine existing atlas sources; no effects.

Existing subtype candidates define diagnostic pools, not exact NSCLC mapping.
"""
import hashlib
import json
import os
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
    project = root.parents[2]
    atlas_run = project / 'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    assert atlas_run.exists()
    atlas = atlas_run / 'outputs/immune_validation_authority_v1.zarr'
    release_file = atlas_run / 'manifests/immune_authority_release_public.json'
    release = json.loads(release_file.read_text())
    assert release['all_checks_pass'] and release['technical_authority_released']
    z = zarr.open_group(str(atlas), mode='r')
    assert tuple(z['layers']['counts'].attrs['shape']) == (2670057, 49563)
    columns = ['dataset', 'library_id', 'sampleID', 'node_id', 'major', 'patient_block', 'patient_identity_verified', 'timepoint', 'tissue', 'treatment_arm', 'capture_scope', 'ex_vivo_stimulation']
    metadata = pd.DataFrame({name: np.asarray(read_elem(z['obs'][name])).astype(str) for name in columns})
    catalog_text = read_elem(z['uns']['subcluster_catalog_json'])
    if isinstance(catalog_text, np.ndarray):
        catalog_text = catalog_text.item()
    catalog = pd.DataFrame(json.loads(catalog_text))
    assert len(catalog) == 270 and catalog.node_id.is_unique
    pools = dict(Myeloid=set(catalog.loc[catalog.major.eq('Myeloid'), 'node_id']), CD8_candidate=set(catalog.loc[catalog.major.eq('T') & catalog.subtype_candidate.eq('CD8_Effector') & catalog.subtype_evidence_status.eq('multigene_subtype_candidate'), 'node_id']), Treg_candidate=set(catalog.loc[catalog.major.eq('T') & catalog.subtype_candidate.eq('Treg') & catalog.subtype_evidence_status.eq('multigene_subtype_candidate'), 'node_id']))
    assert all(len(values) for values in pools.values())
    source_order = np.asarray(read_elem(z['uns']['source_order'])).astype(str).tolist()
    measured = np.asarray(read_elem(z['varm']['measured_in_source']), bool)
    gene_symbols = np.asarray(read_elem(z['var']['gene_symbol'])).astype(str)
    aliases = np.asarray(read_elem(z['var']['source_gene_aliases'])).astype(str)
    identity_status = np.asarray(read_elem(z['var']['feature_identity_status'])).astype(str)
    canonical = np.asarray(read_elem(z['var']['canonical_gene_identity'])).astype(str)
    assert measured.shape == (len(gene_symbols), len(source_order))
    original_root = root / 'NSCLC_automatic_ligand_response_review_20261016_v1'
    original_file = original_root / 'private_all_primary_BH_supported_rows.tsv'
    original_inventory = json.loads((original_root / 'private_terminal.json').read_text())
    assert sha(original_file) == original_inventory['files'][original_file.name]
    parents = pd.read_csv(original_file, sep='\t')
    params_path = root / 'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'
    params = pd.read_csv(params_path, sep='\t')
    out = root / 'existing_immune_cohort_candidate_extension_preflight_20261016_v2'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(previous_exact_name_only_preflight_sha256=sha(root / 'existing_immune_cohort_candidate_extension_preflight_20261016_v1/private_terminal.json'), alias_policy_reference='Existing check_HNSCC_frozen_gene_axis_v2.py/public_gene_identity_contract.py; same95% gate and collision exclusion.', created_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__), atlas_release_sha256=sha(release_file), original_evidence_sha256=sha(original_file), private_frozen_parameters_sha256=sha(params_path), original_parents=len(parents), sources=source_order, pool_catalog_sha256=hashlib.sha256(str(catalog_text).encode()).hexdigest(), diagnostic_pool_node_ids={k: sorted(v) for k, v in pools.items()},
        scope='ALL nine archived sources, every original parent, all existing physical libraries and original subtype-candidate pools. No RNA expressions or clinical effect read.',
        availability='V2: source-native unambiguous physical columns, canonical ID/exact symbol plus pipe-separated aliases ONLY for stable_Ensembl_ID. Exclude ambiguous source features. No fuzzy or expression-driven matching. Multiple physical matches or two frozen genes to the same column are unresolved, not chosen by weights/effects. Full signed programme95% L1 availability is technical only.',
        pool_gate='>=20 original cells at BOTH diagnostic endpoint pools within a physical library. Pools are lineage-extension opportunities, not exact replication of six NSCLC author CD8 states or VCAN/MARCO labels. No frozen labels changed.',
        strata='Dataset x original tissue x original timepoint x treatment_arm x capture_scope x ex_vivo_stimulation. Unknown/conflicting metadata retained explicitly. Count reviewed patient blocks and unreviewed library-fallback separately; no automatic clinical or independent validation acceptance.',
        controls='Atlas alone lacks some original clinical nuisances. No silent covariate deletion or external model fitted; any different external model requires explicit prospective scope and comparable discovery-source description.',
        clinical='Archived eligibility can be older than later author-table sidecars. This is structural input review, not an assertion that updated source clinical data are absent.',
        no_new_P_BH_scores_expression_or_download=True, original_CM172_identities_changed=False)
    (out / 'plan_before_external_input_review.json').write_text(json.dumps(plan, indent=2) + '\n')
    libraries = []
    fields = ['sampleID', 'patient_block', 'patient_identity_verified', 'timepoint', 'tissue', 'treatment_arm', 'capture_scope', 'ex_vivo_stimulation']
    for (source, library), g in metadata.groupby(['dataset', 'library_id'], sort=True):
        row = dict(dataset=source, library_id=library, all_cells=len(g))
        for name in fields:
            values = sorted(set(g[name]))
            row[name] = values[0] if len(values) == 1 else 'CONFLICTING_LIBRARY_METADATA'
        for name, ids in pools.items():
            row[name + '_cells'] = int(g.node_id.isin(ids).sum())
        row['CAR_product_or_ex_vivo_scope'] = 'product' in row['capture_scope'].lower() or 'ex_vivo' in row['capture_scope'].lower() or row['ex_vivo_stimulation'].lower() in ['true', 'yes', '1']
        libraries.append(row)
    lib = pd.DataFrame(libraries)
    assert len(lib) == 400 and lib.all_cells.sum() == 2670057
    lib.to_csv(out / 'private_all_library_pool_input_coverage.tsv', sep='\t', index=False)
    genes_needed = set(params.gene)
    for row in parents.itertuples():
        genes_needed.add(row.ligand)
        genes_needed.update(g for receptor in row.qualified_receptors.split('|') for g in receptor.split('_'))
    gene_rows, input_rows, stratum_rows = [], [], []
    for source_index, source in enumerate(source_order):
        native_cols = np.flatnonzero(measured[:, source_index])
        names = {}
        for col in native_cols:
            if canonical[col].startswith('AMBIGUOUS::') or identity_status[col] == 'ambiguous_source_feature':
                continue
            known = {gene_symbols[col], canonical[col]}
            if identity_status[col] == 'stable_Ensembl_ID':
                known.update(aliases[col].split('|'))
            for name in known - {'', 'nan', 'None', '__SOURCE_VALUE_MISSING__'}:
                names.setdefault(name, set()).add(int(col))
        lookup = {}
        for gene in genes_needed:
            indices = sorted(names.get(gene, set()))
            lookup[gene] = int(indices[0]) if len(indices) == 1 else None
            gene_rows.append(dict(dataset=source, gene=gene, exact_source_physical_matches=len(indices), native_column=int(indices[0]) if len(indices) == 1 else -1, status='UNAMBIGUOUS_SOURCE_NATIVE_COLUMN' if len(indices) == 1 else 'SYMBOL_OR_ALIAS_REVIEW_REQUIRED_NOT_BIOLOGICAL_NEGATIVE'))
        source_libs = lib[lib.dataset.eq(source)]
        for parent in parents.itertuples():
            sender_pool = 'Myeloid' if parent.sender_author_subtype.startswith('M\u03c6_') else 'Treg_candidate'
            receiver_pool = 'Myeloid' if parent.receiver_author_subtype.startswith('M\u03c6_') else 'CD8_candidate'
            if parent.receiver_author_subtype.startswith('CD4T_Treg_'):
                receiver_pool = 'Treg_candidate'
            both = source_libs[source_libs[sender_pool + '_cells'].ge(20) & source_libs[receiver_pool + '_cells'].ge(20)]
            row = dict(dataset=source, record=parent.record, ligand=parent.ligand, sender_diagnostic_pool=sender_pool, receiver_diagnostic_pool=receiver_pool, both_pool20_libraries=len(both), reviewed_patient_blocks=int(both.loc[both.patient_identity_verified.eq('True'), 'patient_block'].nunique()), unreviewed_library_units=int((~both.patient_identity_verified.eq('True')).sum()), CAR_or_ex_vivo_libraries=int(both.CAR_product_or_ex_vivo_scope.sum()), ligand_unambiguous_native=lookup[parent.ligand] is not None, exact_NSCLC_endpoints_validated=False, clinical_model_or_independence_accepted=False)
            for role, typ in [('receiver', parent.receiver_type), ('sender', 'typeB' if parent.receiver_type == 'typeA' else 'typeA')]:
                p = params[params.context.eq(parent.context) & params.type.eq(typ) & params.program.eq(parent.program)]
                indices = [lookup[g] for g in p.gene]
                duplicates = {i for i in indices if i is not None and indices.count(i) > 1}
                keep = np.asarray([i is not None and i not in duplicates for i in indices])
                fraction = float(p.weight.abs().to_numpy()[keep].sum() / p.weight.abs().sum())
                row[role + '_programme_unambiguous_native_L1_fraction'] = fraction
                row[role + '_programme95pct_technical_gate'] = fraction >= .95
            receptors = parent.qualified_receptors.split('|')
            row['all_original_receptors_native'] = all(all(lookup[g] is not None for g in r.split('_')) for r in receptors)
            input_rows.append(row)
            strata = ['tissue', 'timepoint', 'treatment_arm', 'capture_scope', 'ex_vivo_stimulation']
            for keys, g in both.groupby(strata, dropna=False, sort=True):
                stratum_rows.append(dict(dataset=source, record=parent.record, **dict(zip(strata, keys)), libraries=len(g), reviewed_patient_blocks=int(g.loc[g.patient_identity_verified.eq('True'), 'patient_block'].nunique()), unreviewed_library_units=int((~g.patient_identity_verified.eq('True')).sum()), unknown_or_conflicting_stratum=any(v in ['', 'unknown', 'Unknown', 'CONFLICTING_LIBRARY_METADATA'] for v in keys), no_effects=True))
    pd.DataFrame(gene_rows).to_csv(out / 'private_all_source_gene_measurement_status.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(input_rows).to_csv(out / 'private_all_source_parent_coverage.tsv', sep='\t', index=False)
    pd.DataFrame(stratum_rows).to_csv(out / 'private_all_source_parent_strata.tsv', sep='\t', index=False)
    result = dict(status='COMPLETE_ALL9_EXISTING_SOURCE_STRUCTURAL_POOL_AND_MEASUREMENT_PREFLIGHT_NOT_ENDPOINT_OR_CLINICAL_ACCEPTANCE', finished_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(__file__), original_parents=len(parents), sources=len(source_order), all_libraries=len(lib), source_parent_rows=len(input_rows), stratum_rows=len(stratum_rows), files={p.name: sha(p) for p in out.iterdir() if p.is_file()}, no_expression_scores_models_or_download=True, independent_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
