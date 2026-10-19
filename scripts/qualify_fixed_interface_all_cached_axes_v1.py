"""All frozen source axes and current original-source masks; no RNA or effects."""
import hashlib
import json
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from NSCLC_external_lineage_inputs_v2 import load

T_MAJORS = {'T', 'T cell', 'CD4 T cell', 'CD8 T cell'}
MY_MAJORS = {'Myeloid', 'Myeloid cell', 'Macrophage', 'Monocyte', 'Dendritic cell'}
CONTEXT = ['arm', 'dataset', 'cancer_label', 'material_class', 'phase_class', 'capture_scope', 'ex_vivo']


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def stable(value):
    s = str(value)
    return s.split('.')[0] if re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?', s) else s


def roles(v, present, symbols, identifiers):
    fields = [k for k in ['gene', 'gene_symbol', 'gene_id_original', 'gene_name_original',
                          'feature_name', 'canonical_gene_identity'] if k in v]
    arrays = [np.asarray(v.index, dtype=str)] + [np.asarray(v[k], dtype=str) for k in fields]
    n = len(v)
    status = np.asarray(v.feature_identity_status, dtype=str) if 'feature_identity_status' in v else np.full(n, '', dtype=str)
    canonical = np.asarray(v.canonical_gene_identity, dtype=str) if 'canonical_gene_identity' in v else np.full(n, '', dtype=str)
    valid = present & (status != 'ambiguous_source_feature') & ~np.char.startswith(canonical, 'AMBIGUOUS::') & ~np.char.startswith(arrays[0], 'AMBIGUOUS::')
    aliases = v.source_gene_aliases.astype(str).str.split('|', regex=False).tolist() if 'source_gene_aliases' in v else [[] for _ in range(n)]
    actual = [set() for _ in symbols]
    for i in np.flatnonzero(valid):
        tokens = {a[i] for a in arrays}
        if status[i] == 'stable_Ensembl_ID':
            tokens.update(aliases[i])
        tokens |= {stable(t) for t in tokens}
        for k, (symbol, identifier) in enumerate(zip(symbols, identifiers)):
            if symbol in tokens or identifier in tokens:
                actual[k].add(int(i))
    # Independent vector membership; source I/O and allowed aliases are shared.
    stripped = [np.array([stable(x) for x in a], dtype=str) for a in arrays]
    for k, (symbol, identifier) in enumerate(zip(symbols, identifiers)):
        matched = np.zeros(n, bool)
        for raw, normalized in zip(arrays, stripped):
            matched |= (raw == symbol) | (normalized == identifier)
        alias_hit = np.array([status[i] == 'stable_Ensembl_ID' and any(t == symbol or stable(t) == identifier for t in aliases[i]) for i in range(n)])
        expected = set(np.flatnonzero(valid & (matched | alias_hit)).tolist())
        assert actual[k] == expected
    cardinalities = [len(x) for x in actual]
    unique = [next(iter(x)) for x in actual if len(x) == 1]
    return cardinalities, all(n == 1 for n in cardinalities) and len(set(unique)) == len(symbols)


def main():
    config = json.load(sys.stdin)
    symbols = config['symbols']; identifiers = config['stable_ids']
    assert len(symbols) == len(set(symbols)) == len(identifiers) == 3
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    def save(name, value):
        (out / name).write_text(json.dumps(value, indent=2) + '\n')
    data = load()
    lib = data['lib']; counts = data['counts']
    entries = lib[['arm', 'dataset']].drop_duplicates().sort_values(['arm', 'dataset'])
    assert len(entries) == 61 and len(lib) == 1690
    contract_path = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    assert sha(contract_path) == '8d7f54a385389d0ecbc825155df727ed065cae2ed1c3f3ed584c1a7896647e71'
    records = json.loads(contract_path.read_text())['records']
    assert len(records) == 21
    save('private_plan_before_all_source_feature_and_parent_counts.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(),
         source_sha256=sha(__file__), source_contract_sha256=sha(contract_path), frozen_library_proof=data['hashes'],
         private_fixed_interface=config, current_atlas_arm_sources=61, original_contract_records=21,
         additional_original_NSCLC_calibration_axes=1, T_parent_major_labels=sorted(T_MAJORS),
         broad_myeloid_parent_labels=sorted(MY_MAJORS), contexts=CONTEXT,
         minimum_parent_cells=20, all_contexts_retained=True, all_gene_aliases_fixed_to_existing_authority=True,
         scope='Metadata readiness/parent capture upper bound only; not RNA detection, endpoint or patient promotion. Full per-material and gene metadata RAM only; source-level summaries private.',
         no_new_clinical_source_or_RNA_download=True, output_cap_bytes=4 * 1024 ** 2))
    legacy = pd.read_csv(root / 'all_CM_state_native_input_preflight_20261014_v1/all_source_native_scope.tsv', sep='\t').set_index('gse_id')
    genedir = root.parents[2] / 'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1/by_gse'
    z = zarr.open_group(str(data['atlas']), mode='r')
    avar = read_elem(z['var'])
    order = np.asarray(read_elem(z['uns/source_order'])).astype(str).tolist()
    mask = np.asarray(read_elem(z['varm/measured_in_source']), dtype=bool)
    assert mask.shape == (len(avar), len(order))
    atlas = []
    axis_hashes = {}
    for n, row in enumerate(entries.itertuples(index=False)):
        if row.arm == 'frozen_original_annotation':
            path = genedir / row.dataset / 'genes.tsv'
            assert sha(path) == legacy.loc[row.dataset, 'native_gene_axis_sha256']
            v = pd.read_csv(path, sep='\t', keep_default_na=False)
            assert v.gene.is_unique and v.core_column.is_unique
            v.index = v.gene.to_numpy(str)
            measured = ~v.gene.str.startswith('AMBIGUOUS::').to_numpy()
            axis_sha = sha(path)
        else:
            assert row.arm == 'archived_immune_atlas'
            v = avar; measured = mask[:, order.index(row.dataset)]
            axis_sha = hashlib.sha256(measured.tobytes()).hexdigest()
        cardinalities, ok = roles(v, measured, symbols, identifiers)
        atlas.append(dict(arm=row.arm, dataset=row.dataset, gene_axis_rows=len(v), source_measured_features=int(measured.sum()),
                          role0_columns=cardinalities[0], role1_columns=cardinalities[1], role2_columns=cardinalities[2],
                          all_fixed_roles_unique=ok, feature_axis_sha256=axis_sha, RNA_detection_or_clinical_admitted=False))
        axis_hashes[row.arm + '|' + row.dataset] = axis_sha
        save('private_progress.json', dict(stage='ALL61_CURRENT_AXIS_METADATA', completed=n + 1, total=61,
             observed_UTC=datetime.now(timezone.utc).isoformat()))
    # Source identity for the original NSCLC calibration is retained separately.
    nv_path = root / 'NSCLC_full_native_gene_units_20261016_v1/private_native_gene_axis.tsv'
    nv = pd.read_csv(nv_path, sep='\t', keep_default_na=False)
    nv.index = nv.gene.to_numpy(str)
    cards, ok = roles(nv, np.ones(len(nv), bool), symbols, identifiers)
    save('private_original_NSCLC_calibration_axis.json', dict(dataset='GSE243013', calibration_only=True,
         source_gene_axis_sha256=sha(nv_path), fixed_role_columns=cards, all_fixed_roles_unique=ok))
    table = pd.DataFrame(atlas)
    table.to_csv(out / 'private_all61_current_source_axis_qualification.tsv', sep='\t', index=False)
    # Inspect current original source axes using declared counts/var keys and masks.
    physical_cache = {}; physical = []
    for record in records:
        key = (record['path'], record['var_key'])
        if key not in physical_cache:
            p = Path(record['path'])
            h = h5py.File(p, 'r') if p.suffix == '.h5ad' else zarr.open_group(str(p), mode='r')
            v = read_elem(h[record['var_key']]); sm = None; so = None
            if 'varm' in h and 'measured_in_source' in h['varm']:
                sm = np.asarray(read_elem(h['varm/measured_in_source']), bool)
                so = np.asarray(read_elem(h['uns/source_order'])).astype(str).tolist()
                assert sm.shape == (len(v), len(so))
            if p.suffix == '.h5ad':
                h.close()
            physical_cache[key] = (v, sm, so)
        v, sm, so = physical_cache[key]
        measured = sm[:, so.index(record.get('new_parent_dataset', record['dataset']))] if sm is not None else np.ones(len(v), bool)
        if 'feature_measured_in_source' in v:
            assert set(v.feature_measured_in_source.unique()) <= {True, False}
            measured = measured & v.feature_measured_in_source.to_numpy(bool)
        cards, ok = roles(v, measured, symbols, identifiers)
        matched = table[table.arm.eq('archived_immune_atlas') & table.dataset.eq(record['dataset'])]
        same = None if not len(matched) else bool(matched.iloc[0][['role0_columns', 'role1_columns', 'role2_columns']].astype(int).tolist() == cards)
        physical.append(dict(dataset=record['dataset'], source_role=record['role'],
             physical_axis_hash=hashlib.sha256((record['path'] + '|' + record['var_key']).encode()).hexdigest(),
             declared_var_key=record['var_key'], declared_count_key=record['counts_key'],
             source_mask_applied=sm is not None or 'feature_measured_in_source' in v,
             role0_columns=cards[0], role1_columns=cards[1], role2_columns=cards[2],
             all_fixed_roles_unique=ok, current_atlas_cardinality_agreement=same,
             prior_role_and_clinical_eligibility_preserved=True))
    direct = pd.DataFrame(physical)
    direct.to_csv(out / 'private_all21_original_contract_axis_qualification.tsv', sep='\t', index=False)
    icb = direct[direct.dataset.str.startswith('ICB_')]
    assert len(icb) == 8 and icb.physical_axis_hash.nunique() == 1
    save('private_ICB_shared_original_axis_inheritance.json', dict(original_contract_source_positions=8,
         unique_physical_feature_axes=1, source_datasets=icb.dataset.tolist(),
         all_source_role_cardinalities_identical=bool(icb[['role0_columns', 'role1_columns', 'role2_columns']].drop_duplicates().shape[0] == 1),
         source_role_cardinalities=icb[['role0_columns', 'role1_columns', 'role2_columns']].iloc[0].astype(int).tolist(),
         no_source_expression_or_biological_absence_claim=True))
    # Context counts are upper bounds using frozen parents, never new cell identities.
    by_parent = counts.assign(T=counts.major.isin(T_MAJORS).astype(int) * counts.n_cells,
                              My=counts.major.isin(MY_MAJORS).astype(int) * counts.n_cells)
    summed = by_parent.groupby(['arm', 'dataset', 'library_id'])[['T', 'My']].sum()
    manual = defaultdict(lambda: [0, 0])
    for r in counts.itertuples(index=False):
        target = manual[(r.arm, r.dataset, r.library_id)]
        target[0] += int(r.n_cells) if r.major in T_MAJORS else 0
        target[1] += int(r.n_cells) if r.major in MY_MAJORS else 0
    parent = []
    for r in lib.itertuples(index=False):
        key = (r.arm, r.dataset, r.library_id)
        a = summed.loc[key].to_numpy(int) if key in summed.index else np.zeros(2, int)
        assert a.tolist() == manual[key]
        parent.append(a)
    lib = lib.copy(); lib[['T_parent_cells', 'My_parent_cells']] = np.asarray(parent)
    lib['both_parents_ge20'] = lib.T_parent_cells.ge(20) & lib.My_parent_cells.ge(20)
    scope = []
    for key, g in lib.groupby(CONTEXT, sort=True, dropna=False):
        both = g[g.both_parents_ge20]
        reviewed = both[both.reviewed_patient]
        scope.append(dict(zip(CONTEXT, key), libraries=len(g), same_material_both_parents_ge20=len(both),
                          reviewed_patient_blocks_with_one_eligible_material=int(reviewed.patient_block.nunique()),
                          unreviewed_eligible_library_units=int((~both.reviewed_patient).sum()),
                          RNA_or_fine_identity_or_endpoint_admission=False))
    contexts = pd.DataFrame(scope)
    assert contexts.libraries.sum() == len(lib)
    contexts.to_csv(out / 'private_all_source_context_parent_capture_upper_bounds.tsv', sep='\t', index=False)
    save('private_metadata_input_and_audit_proof.json', dict(frozen_input_hashes=data['hashes'], axis_hashes=axis_hashes,
         full_atlas_gene_metadata_sha256=hashlib.sha256(avar.to_json(orient='split').encode()).hexdigest(),
         full_source_mask_sha256=hashlib.sha256(mask.tobytes()).hexdigest(),
         all61_plus_calibration_role_sets_double_calculated=True, all21_original_axis_roles_double_calculated=True,
         all1690_material_parent_counts_independent_Counter_equal=True,
         original_contract_path_var_mask_roles_preserved=True, all_context_library_denominators_conserved=True))
    files = {p.name: dict(bytes=p.stat().st_size, sha256=sha(p)) for p in out.iterdir() if p.is_file()}
    assert sum(v['bytes'] for v in files.values()) < 4 * 1024 ** 2
    terminal = dict(status='COMPLETE_ALL61_PLUS_CALIBRATION_AND21_ORIGINAL_AXIS_METADATA_NO_EFFECTS',
         completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), files=files,
         current_source_axes=61, additional_calibration_axes=1, original_contract_records=21,
         original_physical_path_var_axes=len(physical_cache), library_records_checked=1690,
         source_context_positions=len(contexts), new_clinical_data_or_expression_read=False,
         RNA_presence_or_clinical_identity_not_inferred=True)
    save('private_terminal.json', terminal)
    print(json.dumps(dict(status=terminal['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
