"""Complete prespecified group lineage descriptions and schema-only inventory.

All scientific outputs stay private. No sequences, protein matrix, clinical
outcomes, new RNA download, threshold tuning or formal label writes.
"""
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import zarr
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

PANELS = {
    'T_lineage': 'CD3D CD3E CD3G CD247 CD2 LCK'.split(),
    'gamma_constant': 'TRDC TRGC1 TRGC2'.split(),
    'gamma_variable': 'TRDV1 TRDV2 TRDV3 TRGV9'.split(),
    'shared_alpha_delta_variable': ['TRAV29DV5'],
    'alpha_beta': 'TRAC TRBC1 TRBC2'.split(),
    'NK_innate': 'FCER1G TYROBP KLRD1 NCR1'.split(),
    'cytotoxic': 'NKG7 GNLY PRF1 GZMB'.split(),
}
MARKERS = sum(PANELS.values(), [])
GROUPS = ['Author_Vdelta1', 'Author_Vdelta2', 'Author_NK', 'Other_T_NK', 'Non_T_NK']
TCR = r'airr|tcr|vdj|cdr3|clonotyp|productive|(^|_)ir_|(^|_)(tra|trb|trg|trd)_|(^|_)chain($|_)|[vjdc]_gene'
PROTEIN = r'protein|(^|_)adt($|_)|antibody|facs|cytometry|tcrgd'


def shape_only(node):
    if hasattr(node, 'shape'):
        return dict(shape=[int(v) for v in node.shape], encoding=str(node.attrs.get('encoding-type', 'array')))
    return dict(encoding=str(node.attrs.get('encoding-type', 'group')), schema_keys=sorted(node.keys()))


def metadata_schema(record):
    path = Path(record['path'])
    h = h5py.File(path, 'r') if path.suffix == '.h5ad' else zarr.open_group(str(path), mode='r')
    try:
        obs = h['obs']
        index_name = obs.attrs['_index']
        index_name = index_name.decode() if isinstance(index_name, bytes) else str(index_name)
        index = obs[index_name]
        n = int(index.shape[0] if hasattr(index, 'shape') else index['codes'].shape[0])
        selector = np.load(record['row_positions_ref'], allow_pickle=False)
        assert selector.ndim == 1 and np.issubdtype(selector.dtype, np.integer)
        assert len(selector) == record['n_cells'] and (not len(selector) or (selector.min() >= 0 and selector.max() < n))
        fields = sorted(obs.keys())
        obsm = {key: shape_only(h['obsm'][key]) for key in h.get('obsm', {})}
        uns = sorted(h['uns'].keys()) if 'uns' in h else []
        var = sorted(h['var'].keys()) if 'var' in h else []
        mod = sorted(h['mod'].keys()) if 'mod' in h else []
        names = fields + list(obsm) + uns + var + mod
        return dict(dataset=record['dataset'], source_role=record['role'], raw_metadata_rows=n, selected_rows=len(selector), source_path=str(path), row_selector_sha256=sha(record['row_positions_ref']), obs_fields=fields, obsm_schema=obsm, uns_fields=uns, var_fields=var, modality_fields=mod, potential_TCR_schema=sorted({v for v in names if re.search(TCR, v, re.I)}), potential_protein_schema=sorted({v for v in names if re.search(PROTEIN, v, re.I)}), scope='Schema only; no chain/productive/sequence/protein content or exact cell bridge accepted')
    finally:
        if isinstance(h, h5py.File):
            h.close()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 2 * 1024**2)
    evidence = root / 'gamma_candidate_evidence_reconciliation_audit_20261019_v1/private_terminal.json'
    raw_audit = root / 'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json'
    assert sha(evidence) == 'a5dabd6d1003f0b4576ebe2c9dcdeaa003b875a10fd55bea344a60852d3ab706'
    assert sha(raw_audit) == 'bb20365ee83758a7df8478ceadede2596fe749928dab5765e3bc0866d72dd8c1'
    prior = root / 'NSCLC_cached_annotation_feature_audit_20261019_v1/private_terminal.json'
    assert sha(prior) == '4b9ae548f959684e4d1653ba0d36fc0a38cf234547601f6784c8ccf88d4b3f80'
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    t = json.loads((native / 'private_terminal.json').read_text())
    assert json.loads(raw_audit.read_text())['producer_terminal_sha256'] == sha(native / 'private_terminal.json')
    contract = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    config = json.loads(contract.read_text())
    previous = json.loads((root / 'NSCLC_cached_annotation_feature_qualification_20261019_v1/private_terminal.json').read_text())
    assert sha(contract) == previous['source_contract_sha256'] and len(config['records']) == 21
    out = root / 'gamma_group_lineage_orthogonal_description_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_lineage_statistic_or_TCR_schema.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), resource_guard=guard, previous_evidence_sha256=sha(evidence), native_raw_audit_sha256=sha(raw_audit), source_contract_sha256=sha(contract), panels=PANELS, author_groups=GROUPS, paired_author_record_minimum_cells_each=20, paired_comparisons=4, paired_marker_positions=100, scope='Every51 original author types/all11701units, complete fixed25markers,21 existing source schemas and one-level adjacent cache folders. No outcome/sequence/protein values, inference or formal classification.', full_cache_raw_replay=False))
    arrays = {}
    cache_hashes = {}
    genes = pd.read_csv(root / 'GSE243013_genes.csv.gz', keep_default_na=False).geneSymbol.str.strip().str.upper()
    assert len(genes) == 31831 and genes.is_unique and len(MARKERS) == len(set(MARKERS)) == 25
    columns = pd.Index(genes).get_indexer(MARKERS)
    assert (columns >= 0).all()
    for name in ['sum_raw_counts.npy', 'positive_cell_counts.npy', 'mean_CP10K.npy', 'mean_logCP10K.npy']:
        path = native / name
        assert sha(path) == t['outputs'][name]['sha256']
        cache_hashes[name] = t['outputs'][name]['sha256']
        a = np.load(path, mmap_mode='r', allow_pickle=False)
        assert a.shape == (11701, 31831)
        arrays[name] = np.array(a[:, columns], copy=True)
        save('private_progress.json', dict(stage='CHECKING_ACCEPTED_CACHE_HASHES_AND_FIXED_MARKER_SLICES', verified_arrays=len(arrays), observed_UTC=datetime.now(timezone.utc).isoformat()))
    upath = native / 'private_all_native_units_and_quality.tsv'
    assert sha(upath) == t['outputs'][upath.name]['sha256']
    u = pd.read_csv(upath, sep='\t', keep_default_na=False).sort_values('unit_index').reset_index(drop=True)
    assert u.unit_index.tolist() == list(range(11701)) and not u.duplicated(['sampleID', 'sub_cell_type']).any()
    n = u.n_cells.to_numpy(np.int64)
    assert n.sum() == 1254749 and (n > 0).all() and u.sub_cell_type.nunique() == 51
    positive = arrays['positive_cell_counts.npy']
    assert (positive >= 0).all() and (positive <= n[:, None]).all()
    raw = arrays['sum_raw_counts.npy']
    cp = arrays['mean_CP10K.npy']
    logcp = arrays['mean_logCP10K.npy']
    assert (raw >= positive).all() and np.isfinite(cp).all() and np.isfinite(logcp).all()
    assert (cp >= 0).all() and (logcp >= 0).all() and (logcp <= np.log1p(cp) + 1e-10).all()
    for a in [raw, cp, logcp]:
        assert np.array_equal(a > 0, positive > 0)
    unit = np.load(native / 'private_cell_unit_index.npy', allow_pickle=False)
    h = dict(np.load(root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1/private_candidate_identity_sidecar.npz', allow_pickle=False))
    rule = np.zeros((len(u), 4), np.int64)
    for j in range(4):
        rule[:, j] = np.bincount(unit[h['native_row'][(h['rule_bits'] & (1 << j)) > 0]], minlength=len(u))
    core = positive[:, pd.Index(MARKERS).get_indexer(['CD3D', 'CD3E', 'CD3G', 'CD247', 'TRDC'])]
    lower = np.maximum(0, core.sum(axis=1) - 4 * n)
    upper = core.min(axis=1)
    assert (rule[:, 0] >= lower).all() and (rule[:, 0] <= upper).all()
    assert (rule[:, 2] <= rule[:, 0]).all() and (rule[:, 3] <= rule[:, 1]).all()
    role = np.full(len(u), 4, np.int8)
    role[u.major_cell_type.eq('T/NK cell')] = 3
    role[u.sub_cell_type.str.startswith('NK_')] = 2
    role[u.sub_cell_type.eq('T_gdT_TRDV2')] = 1
    role[u.sub_cell_type.eq('T_gdT_TRDV1')] = 0
    rows = []
    quality = []
    gene_panel = {v: k for k, vv in PANELS.items() for v in vv}
    for label, g in u.groupby('sub_cell_type', sort=True):
        ix = g.index.to_numpy(int)
        count = n[ix]
        assert len(set(g.major_cell_type)) == 1
        total = int(count.sum())
        eligible = ix[count >= 20]
        quality.append(dict(source='GSE243013', original_author_type=label, major=g.major_cell_type.iloc[0], author_records=len(g), author_records_atleast20=len(eligible), native_cells=total, average_RNA=float(g.native_total_RNA.sum() / total), average_detected_genes=float(g.native_detected_gene_sum.sum() / total), core5_count1=int(rule[ix, 0].sum()), supported_count1=int(rule[ix, 1].sum()), core5_count2=int(rule[ix, 2].sum()), supported_count2=int(rule[ix, 3].sum()), core5_Frechet_lower_count=int(lower[ix].sum()), core5_Frechet_upper_count=int(upper[ix].sum()), reviewed_independent_patients=0))
        for j, marker in enumerate(MARKERS):
            detections = positive[ix, j]
            frac = positive[eligible, j] / n[eligible]
            means = cp[eligible, j]
            loo = (detections.sum() - detections)[total - count > 0] / (total - count)[total - count > 0]
            rows.append(dict(source='GSE243013', original_author_type=label, marker=marker, fixed_panel=gene_panel[marker], native_cells=total, author_units=len(ix), eligible_units=len(eligible), raw_count_sum=int(raw[ix, j].sum()), detected_cells=int(detections.sum()), detected_fraction=float(detections.sum() / total), cell_mean_CP10K=float(cp[ix, j] @ count / total), cell_mean_logCP10K=float(logcp[ix, j] @ count / total), pseudobulk_CP10K=float(raw[ix, j].sum() * 10000 / g.native_total_RNA.sum()), eligible_unit_detection_median=float(np.median(frac)) if len(frac) else None, eligible_unit_detection_q25=float(np.quantile(frac, .25)) if len(frac) else None, eligible_unit_detection_q75=float(np.quantile(frac, .75)) if len(frac) else None, eligible_unit_CP10K_median=float(np.median(means)) if len(means) else None, leave_one_author_record_detection_min=float(loo.min()) if len(loo) else None, leave_one_author_record_detection_max=float(loo.max()) if len(loo) else None))
    pd.DataFrame(rows).to_csv(out / 'private_all51_author_type_fixed25_marker_profiles.tsv', sep='\t', index=False)
    pd.DataFrame(quality).to_csv(out / 'private_all51_author_type_quality_and_coexpression_bounds.tsv', sep='\t', index=False)
    sample = pd.Categorical(u.sampleID).codes
    assert sample.max() + 1 == 243
    key = sample * 5 + role
    nn = np.bincount(key, weights=n, minlength=243 * 5).reshape(243, 5)
    pp = np.zeros((243 * 5, 25), np.int64)
    cc = np.zeros((243 * 5, 25), float)
    np.add.at(pp, key, positive)
    np.add.at(cc, key, cp * n[:, None])
    pp = pp.reshape(243, 5, 25)
    cc = cc.reshape(243, 5, 25)
    paired = []
    for a in [0, 1]:
        for b in [2, 3]:
            valid = (nn[:, a] >= 20) & (nn[:, b] >= 20)
            delta = pp[valid, a] / nn[valid, a, None] - pp[valid, b] / nn[valid, b, None]
            expr = cc[valid, a] / nn[valid, a, None] - cc[valid, b] / nn[valid, b, None]
            for j, marker in enumerate(MARKERS):
                paired.append(dict(source='GSE243013', reference_gamma_group=GROUPS[a], comparison_group=GROUPS[b], marker=marker, comparable_author_records=int(valid.sum()), independent_patients_reviewed=0, median_detection_difference=float(np.median(delta[:, j])) if valid.any() else None, mean_detection_difference=float(np.mean(delta[:, j])) if valid.any() else None, detection_difference_q25=float(np.quantile(delta[:, j], .25)) if valid.any() else None, detection_difference_q75=float(np.quantile(delta[:, j], .75)) if valid.any() else None, fraction_records_positive_difference=float(np.mean(delta[:, j] > 0)) if valid.any() else None, fraction_records_equal_difference=float(np.mean(delta[:, j] == 0)) if valid.any() else None, median_CP10K_difference=float(np.median(expr[:, j])) if valid.any() else None, inferential_test=False))
    assert len(paired) == 100
    pd.DataFrame(paired).to_csv(out / 'private_all100_within_author_record_marker_contrasts.tsv', sep='\t', index=False)
    schemas = [metadata_schema(v) for v in config['records']]
    save('private_all21_source_orthogonal_schema_inventory.json', schemas)
    metadata = root / 'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz'
    headers = pd.read_csv(metadata, nrows=0).columns.tolist()
    save('private_original_NSCLC_orthogonal_metadata_header.json', dict(source='GSE243013', metadata_sha256=sha(metadata), fields=headers, potential_TCR_fields=[v for v in headers if re.search(TCR, v, re.I)], potential_protein_fields=[v for v in headers if re.search(PROTEIN, v, re.I)], no_metadata_values_read=True))
    folders = {root} | {Path(v['path']).parent for v in config['records']}
    inventory = []
    for folder in sorted(folders):
        for path in sorted(folder.iterdir()):
            if path.is_file() and re.search(r'tcr|vdj|airr|contig|clonotyp|(^|_)adt|protein', path.name, re.I):
                info = dict(folder=str(folder), filename=path.name, bytes=path.stat().st_size, contents_read=False, header_fields=[])
                low = path.name.lower()
                if low.endswith(('.csv', '.csv.gz', '.tsv', '.tsv.gz')):
                    try:
                        info['header_fields'] = pd.read_csv(path, sep='\t' if '.tsv' in low else ',', nrows=0).columns.tolist()
                        info['header_status'] = 'HEADER_ONLY_READ_NO_VALUES'
                    except (pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeError, EOFError, OSError) as error:
                        info['header_status'] = 'HEADER_UNAVAILABLE_' + type(error).__name__
                inventory.append(info)
    save('private_complete_adjacent_orthogonal_file_opportunities.json', dict(one_level_folders=[str(v) for v in sorted(folders)], files=inventory, scope_excludes_other_server_locations=True, no_sequence_or_expression_values_read=True))
    save('private_cache_and_fixed_marker_slice_provenance.json', dict(full_current_cache_sha256_verified=cache_hashes, native_gene_axis_sha256=sha(root / 'GSE243013_genes.csv.gz'), native_unit_table_sha256=sha(upath), fixed_marker_native_columns=columns.tolist(), selected_cache_slices_sha256={name: __import__('hashlib').sha256(a.tobytes()).hexdigest() for name, a in arrays.items()}, parent_cache_full_independent_raw_replayed=False))
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file()) < 2 * 1024**2
    receipt = dict(status='COMPLETE_ALL51_FIXED25_GROUP_PROFILES_ALL100_PAIRED_DESCRIPTIONS_AND21_ORTHOGONAL_SCHEMAS_AUDIT_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), native_author_types=51, marker_positions=1275, paired_marker_positions=100, cached_source_entries=21, cache_hashes_fully_reverified=True, native_RNA_or_TCR_download=False, TCR_sequence_or_protein_values_read=False, clinical_test_or_formal_identity_adopted=False, resource_guard=guard, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
