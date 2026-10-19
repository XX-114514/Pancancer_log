"""Full cached-statistic and schema audit; not new independent raw/TCR truth."""
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import zarr
from describe_gamma_group_lineage_and_cached_orthogonal_scope_v1 import MARKERS, PANELS, GROUPS
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait


def equal(actual, expected):
    if pd.isna(expected):
        assert pd.isna(actual)
    else:
        np.testing.assert_allclose(float(actual), float(expected), rtol=1e-10, atol=1e-12)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 1024**2)
    src = root / 'gamma_group_lineage_orthogonal_description_20261019_v1'
    out = root / 'gamma_group_lineage_orthogonal_description_audit_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    terminal = wait(root, out, 'describe_gamma_group_lineage_and_cached_orthogonal_scope_v1', src / 'private_terminal.json', {'COMPLETE_ALL51_FIXED25_GROUP_PROFILES_ALL100_PAIRED_DESCRIPTIONS_AND21_ORTHOGONAL_SCHEMAS_AUDIT_PENDING'})
    for name, info in terminal['files'].items():
        assert sha(src / name) == info['sha256']
    plan = json.loads((src / 'plan_before_lineage_statistic_or_TCR_schema.json').read_text())
    proof = json.loads((src / 'private_cache_and_fixed_marker_slice_provenance.json').read_text())
    assert plan['panels'] == PANELS and plan['paired_marker_positions'] == 100
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    nt = json.loads((native / 'private_terminal.json').read_text())
    u = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False).set_index('unit_index').loc[range(11701)].reset_index()
    assert sha(native / 'private_all_native_units_and_quality.tsv') == proof['native_unit_table_sha256']
    genes = pd.read_csv(root / 'GSE243013_genes.csv.gz', keep_default_na=False).geneSymbol.str.upper().str.strip().tolist()
    assert len(genes) == len(set(genes)) == 31831
    columns = [genes.index(g) for g in MARKERS]
    assert columns == proof['fixed_marker_native_columns']
    data = {}
    for filename, expected in proof['selected_cache_slices_sha256'].items():
        assert proof['full_current_cache_sha256_verified'][filename] == nt['outputs'][filename]['sha256']
        a = np.load(native / filename, mmap_mode='r', allow_pickle=False)
        assert (native / filename).stat().st_size == nt['outputs'][filename]['bytes']
        sliced = np.column_stack([a[:, col] for col in columns])
        assert hashlib.sha256(sliced.tobytes()).hexdigest() == expected
        data[filename] = sliced
    profiles = pd.read_csv(src / 'private_all51_author_type_fixed25_marker_profiles.tsv', sep='\t').set_index(['original_author_type', 'marker'])
    qualities = pd.read_csv(src / 'private_all51_author_type_quality_and_coexpression_bounds.tsv', sep='\t').set_index('original_author_type')
    assert len(profiles) == 1275 and len(qualities) == 51
    assert set(profiles.index) == {(x, g) for x in u.sub_cell_type.unique() for g in MARKERS}
    unit = np.load(native / 'private_cell_unit_index.npy', allow_pickle=False)
    np.testing.assert_array_equal(np.bincount(unit, minlength=len(u)), u.n_cells)
    hard = dict(np.load(root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1/private_candidate_identity_sidecar.npz', allow_pickle=False))
    counts = np.zeros((len(u), 4), int)
    for j in range(4):
        selected_units = unit[hard['native_row'][np.bitwise_and(hard['rule_bits'], 2**j) != 0]]
        for key, amount in pd.Series(selected_units).value_counts().items():
            counts[int(key), j] = int(amount)
    corepos = data['positive_cell_counts.npy'][:, [MARKERS.index(g) for g in ['CD3D', 'CD3E', 'CD3G', 'CD247', 'TRDC']]]
    low = np.clip(corepos.sum(axis=1) - 4 * u.n_cells.to_numpy(), 0, None)
    high = np.min(corepos, axis=1)
    assert (counts[:, 0] >= low).all() and (counts[:, 0] <= high).all()
    for label, g in u.groupby('sub_cell_type'):
        q = qualities.loc[label]
        ix = g.unit_index.to_numpy(int)
        ncells = int(g.n_cells.sum())
        assert q.native_cells == ncells and q.author_records == len(g) and q.author_records_atleast20 == g.n_cells.ge(20).sum()
        assert q.major == g.major_cell_type.iloc[0] and q.reviewed_independent_patients == 0
        equal(q.average_RNA, g.native_total_RNA.sum() / ncells)
        equal(q.average_detected_genes, g.native_detected_gene_sum.sum() / ncells)
        assert [q.core5_count1, q.supported_count1, q.core5_count2, q.supported_count2] == counts[ix].sum(axis=0).tolist()
        assert q.core5_Frechet_lower_count == low[ix].sum() and q.core5_Frechet_upper_count == high[ix].sum()
    paired = pd.read_csv(src / 'private_all100_within_author_record_marker_contrasts.tsv', sep='\t').set_index(['reference_gamma_group', 'comparison_group', 'marker'])
    assert len(paired) == 100
    role = pd.Series(np.where(u.major_cell_type.eq('T/NK cell'), 'Other_T_NK', 'Non_T_NK'), index=u.index, dtype=object)
    role[u.sub_cell_type.str.startswith('NK_')] = 'Author_NK'
    role[u.sub_cell_type.eq('T_gdT_TRDV1')] = 'Author_Vdelta1'
    role[u.sub_cell_type.eq('T_gdT_TRDV2')] = 'Author_Vdelta2'
    for j, marker in enumerate(MARKERS):
        d = u.assign(raw=data['sum_raw_counts.npy'][:, j], positive=data['positive_cell_counts.npy'][:, j], cp=data['mean_CP10K.npy'][:, j], logcp=data['mean_logCP10K.npy'][:, j], role=role)
        for label, g in d.groupby('sub_cell_type'):
            v = profiles.loc[(label, marker)]
            n = g.n_cells.sum()
            ok = g[g.n_cells >= 20]
            detect = ok.positive.div(ok.n_cells)
            loo = (g.positive.sum() - g.positive).div(n - g.n_cells)
            loo = loo[n - g.n_cells > 0]
            expected = dict(native_cells=n, author_units=len(g), eligible_units=len(ok), raw_count_sum=g.raw.sum(), detected_cells=g.positive.sum(), detected_fraction=g.positive.sum() / n, cell_mean_CP10K=np.average(g.cp, weights=g.n_cells), cell_mean_logCP10K=np.average(g.logcp, weights=g.n_cells), pseudobulk_CP10K=g.raw.sum() * 10000 / g.native_total_RNA.sum(), eligible_unit_detection_median=detect.median(), eligible_unit_detection_q25=detect.quantile(.25), eligible_unit_detection_q75=detect.quantile(.75), eligible_unit_CP10K_median=ok.cp.median(), leave_one_author_record_detection_min=loo.min(), leave_one_author_record_detection_max=loo.max())
            for field, value in expected.items():
                equal(v[field], value)
            assert v.source == 'GSE243013' and v.fixed_panel == next(k for k, panel in PANELS.items() if marker in panel)
        d['weighted_cp'] = d.cp * d.n_cells
        agg = d.groupby(['sampleID', 'role'])[['n_cells', 'positive', 'weighted_cp']].sum()
        for a in GROUPS[:2]:
            for b in ['Author_NK', 'Other_T_NK']:
                aa = agg.xs(a, level='role')
                bb = agg.xs(b, level='role')
                both = aa.join(bb, lsuffix='_a', rsuffix='_b', how='inner')
                both = both[(both.n_cells_a >= 20) & (both.n_cells_b >= 20)]
                diff = both.positive_a / both.n_cells_a - both.positive_b / both.n_cells_b
                expr = both.weighted_cp_a / both.n_cells_a - both.weighted_cp_b / both.n_cells_b
                v = paired.loc[(a, b, marker)]
                assert v.comparable_author_records == len(both) and not v.inferential_test and v.independent_patients_reviewed == 0
                expected = dict(median_detection_difference=diff.median(), mean_detection_difference=diff.mean(), detection_difference_q25=diff.quantile(.25), detection_difference_q75=diff.quantile(.75), fraction_records_positive_difference=diff.gt(0).mean(), fraction_records_equal_difference=diff.eq(0).mean(), median_CP10K_difference=expr.median())
                for field, value in expected.items():
                    equal(v[field], value)
    contract = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    assert sha(contract) == plan['source_contract_sha256']
    records = json.loads(contract.read_text())['records']
    schemas = {v['dataset']: v for v in json.loads((src / 'private_all21_source_orthogonal_schema_inventory.json').read_text())}
    assert len(records) == len(schemas) == 21
    for r in records:
        p = Path(r['path'])
        h = h5py.File(p, 'r') if p.suffix == '.h5ad' else zarr.open_group(str(p), mode='r')
        try:
            v = schemas[r['dataset']]
            assert v['source_role'] == r['role'] and v['row_selector_sha256'] == sha(r['row_positions_ref'])
            assert v['selected_rows'] == len(np.load(r['row_positions_ref'], allow_pickle=False))
            idx = h['obs'].attrs['_index']
            idx = idx.decode() if isinstance(idx, bytes) else str(idx)
            node = h['obs'][idx]
            assert v['raw_metadata_rows'] == (node.shape[0] if hasattr(node, 'shape') else node['codes'].shape[0])
            assert sorted(h['obs'].keys()) == v['obs_fields']
            assert (sorted(h['obsm'].keys()) if 'obsm' in h else []) == sorted(v['obsm_schema'])
            assert (sorted(h['mod'].keys()) if 'mod' in h else []) == v['modality_fields']
            for key, shape in v['obsm_schema'].items():
                item = h['obsm'][key]
                if hasattr(item, 'shape'):
                    assert shape['shape'] == list(item.shape)
                else:
                    assert shape['schema_keys'] == sorted(item.keys())
            for name in ['uns', 'var']:
                assert (sorted(h[name].keys()) if name in h else []) == v[name + '_fields']
            names = v['obs_fields'] + list(v['obsm_schema']) + v['uns_fields'] + v['var_fields'] + v['modality_fields']
            tcr_pattern = r'airr|tcr|vdj|cdr3|clonotyp|productive|(^|_)ir_|(^|_)(tra|trb|trg|trd)_|(^|_)chain($|_)|[vjdc]_gene'
            protein_pattern = r'protein|(^|_)adt($|_)|antibody|facs|cytometry|tcrgd'
            assert v['potential_TCR_schema'] == sorted({s for s in names if re.search(tcr_pattern, s, re.I)})
            assert v['potential_protein_schema'] == sorted({s for s in names if re.search(protein_pattern, s, re.I)})
        finally:
            if isinstance(h, h5py.File):
                h.close()
    inventory = json.loads((src / 'private_complete_adjacent_orthogonal_file_opportunities.json').read_text())
    expected_paths = set()
    for folder in inventory['one_level_folders']:
        for p in Path(folder).iterdir():
            if p.is_file() and re.search(r'tcr|vdj|airr|contig|clonotyp|(^|_)adt|protein', p.name, re.I):
                expected_paths.add((str(p.parent), p.name, p.stat().st_size))
    assert expected_paths == {(v['folder'], v['filename'], v['bytes']) for v in inventory['files']}
    for f in inventory['files']:
        if f.get('header_status') == 'HEADER_ONLY_READ_NO_VALUES':
            p = Path(f['folder']) / f['filename']
            assert pd.read_csv(p, sep='\t' if '.tsv' in p.name.lower() else ',', nrows=0).columns.tolist() == f['header_fields']
        assert not f['contents_read']
    header = json.loads((src / 'private_original_NSCLC_orthogonal_metadata_header.json').read_text())
    p = root / 'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz'
    assert sha(p) == header['metadata_sha256']
    assert pd.read_csv(p, nrows=0).columns.tolist() == header['fields']
    receipt = dict(status='ACCEPTED_ALL1275_CACHED_MARKER_PROFILES_100_PAIRED_DESCRIPTIONS_51_COEXPRESSION_BOUNDS_AND21_SCHEMA_SCOPES', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), scientific_values_public=False, independently_recomputed_marker_positions=1275, independently_recomputed_paired_positions=100, same_cache_slices_exact=True, parent_cache_raw_audit_boundary_unchanged=True, new_raw_RNA_TCR_or_protein_replayed=False, formal_identity_or_clinical_acceptance=False, source_schemas_checked=21, resource_guard=guard)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
