"""Describe fixed panel RNA sources across every frozen major; arrays stay in RAM."""
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from scipy.sparse import csr_matrix


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def aggregate(raw, total, group, groups):
    """Two implementations of full grouped raw/CP10K accumulation."""
    assert (total > 0).all() and (raw >= 0).all()
    cp = raw * (10000.0 / total[:, None])
    n = np.bincount(group, minlength=groups)
    a = np.column_stack([np.bincount(group, weights=raw[:, j], minlength=groups) for j in range(raw.shape[1])])
    b = np.column_stack([np.bincount(group, weights=cp[:, j], minlength=groups) for j in range(raw.shape[1])])
    onehot = csr_matrix((np.ones(len(group)), (group, np.arange(len(group)))), shape=(groups, len(group)))
    np.testing.assert_array_equal(a, np.asarray(onehot @ raw))
    np.testing.assert_allclose(b, np.asarray(onehot @ cp), rtol=1e-11, atol=1e-8)
    np.testing.assert_array_equal(a.sum(axis=0), raw.sum(axis=0))
    return n, a, b


def shares(values, n):
    """Complete-panel shares only; empty scope/zero total gene stays unavailable."""
    if n.sum() == 0 or np.any(values.sum(axis=0) <= 0):
        return None
    weighted = (values / values.sum(axis=0)).mean(axis=1)
    mean = np.divide(values, n[:, None], out=np.zeros_like(values), where=n[:, None] > 0)
    equal = (mean / mean.sum(axis=0)).mean(axis=1)
    assert abs(weighted.sum() - 1) < 1e-10 and abs(equal.sum() - 1) < 1e-10
    return weighted, equal


def fixture():
    raw = np.array([[1, 2], [2, 0], [0, 3], [3, 4]], dtype=np.int64)
    total = np.array([3, 4, 6, 8])
    n, a, b = aggregate(raw, total, np.array([0, 0, 1, 1]), 3)
    np.testing.assert_array_equal(n, [2, 2, 0])
    np.testing.assert_array_equal(a, [[3, 2], [3, 7], [0, 0]])
    assert shares(b, n) is not None and shares(np.zeros((3, 2)), n) is None
    assert shares(a, np.zeros(3)) is None
    assert ((raw >= 1).sum(axis=1) >= 2).sum() == 2
    assert ((raw >= 2).sum(axis=1) >= 2).sum() == 1
    return dict(status='SYNTHETIC_GROUP_SUM_ZERO_EMPTY_AND_TWO_THRESHOLDS_ACCEPTED', real_inputs_used=False)


def main():
    cfg = json.load(sys.stdin)
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent == Path('/dev/shm') and not out.exists()
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    fixed_path = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    assert sha(fixed_path) == cfg['fixed_panel_contract_SHA256']
    assert sha(atlas / 'manifests/immune_authority_release_public.json') == cfg['release_SHA256']
    old = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    audit = root / 'immune_CM_fixed_panel_same_cell_RNA_audit_20261018_v1'
    accepted = json.loads((audit / 'private_independent_acceptance.json').read_text())
    assert accepted['status'].startswith('ACCEPTED_ALL_ORIGINAL_T_CELLS_')
    assert sha(old / 'private_terminal.json') == accepted['producer_terminal_sha256']
    old_manifest = json.loads((old / 'private_terminal.json').read_text())
    for name in ['private_fixed34_gene_T_counts_and_depth.npz', 'private_source_gene_measurement_contract.json']:
        assert sha(old / name) == old_manifest['files'][name]['sha256']
    cache = dict(np.load(old / 'private_fixed34_gene_T_counts_and_depth.npz', allow_pickle=False))
    panels = json.loads(fixed_path.read_text())['panels']
    genes = sorted({g for gg in panels.values() for g in gg})
    assert genes == cache['genes'].tolist()
    names = list(panels)
    pairs = cfg['fixed_panel_pairs']
    assert len(genes) == 34 and len(names) == 7
    z = zarr.open_group(str(atlas / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    counts = z['layers/counts']
    nc, nf = map(int, counts.attrs['shape'])
    assert [nc, nf] == cfg['counts_shape']
    md = {k: np.asarray(read_elem(z['obs'][k])).astype(str) for k in ['dataset', 'library_id', 'major']}
    sources = list(np.asarray(read_elem(z['uns/source_order'])).astype(str))
    majors = sorted(set(md['major']))
    ns, nm, ng = len(sources), len(majors), len(genes)
    source_code = pd.Index(sources).get_indexer(md['dataset'])
    major_code = pd.Index(majors).get_indexer(md['major'])
    assert (source_code >= 0).all() and (major_code >= 0).all()
    group = source_code * nm + major_code
    expected = np.bincount(group, minlength=ns * nm)
    independent = Counter(zip(md['dataset'], md['major']))
    assert all(expected[i * nm + j] == independent[s, m] for i, s in enumerate(sources) for j, m in enumerate(majors))
    cm_path = atlas / 'CM/T_whole_TME_primary_tumor/inputs/sample_metadata_private.tsv'
    cm = pd.read_csv(cm_path, sep='\t', keep_default_na=False)
    assert cm.sampleID.is_unique
    key = pd.Series(md['dataset']) + '::' + pd.Series(md['library_id'])
    assert set(cm.sampleID) <= set(key)
    masks = [np.ones(nc, bool), key.isin(cm.sampleID).to_numpy(), key.isin(cm.loc[cm.timepoint.eq('pre'), 'sampleID']).to_numpy()]
    assert np.array_equal(np.flatnonzero(md['major'] == 'T'), cache['global_original_T_row_index'])
    measured = np.asarray(read_elem(z['varm/measured_in_source']), bool)
    symbol = np.asarray(read_elem(z['var/gene_symbol'])).astype(str)
    canonical = np.asarray(read_elem(z['var/canonical_gene_identity'])).astype(str)
    status = np.asarray(read_elem(z['var/feature_identity_status'])).astype(str)
    aliases = np.asarray(read_elem(z['var/source_gene_aliases'])).astype(str)
    maps = []
    for si in range(ns):
        lookup = defaultdict(set)
        for col in np.flatnonzero(measured[:, si]):
            if canonical[col].startswith('AMBIGUOUS::') or status[col] == 'ambiguous_source_feature':
                continue
            labels = {symbol[col], canonical[col]}
            if status[col] == 'stable_Ensembl_ID':
                labels.update(aliases[col].split('|'))
            for name in labels:
                lookup[name].add(int(col))
        cols = [next(iter(lookup[g])) if len(lookup[g]) == 1 else -1 for g in genes]
        collision = {v for v in cols if v >= 0 and cols.count(v) > 1}
        maps.append(np.array([-1 if v in collision else v for v in cols]))
    old_maps = json.loads((old / 'private_source_gene_measurement_contract.json').read_text())
    assert len(old_maps) == ns
    for si, src in enumerate(sources):
        assert old_maps[si]['dataset'] == src
        np.testing.assert_array_equal(maps[si], old_maps[si]['physical_columns'])
    complete = np.array([[all(maps[si][genes.index(g)] >= 0 for g in panels[p]) for p in names] for si in range(ns)])
    out.mkdir(mode=0o700)
    (out / 'private_frozen_plan.json').write_text(json.dumps(dict(cfg, created_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), CM_material_metadata_SHA256=sha(cm_path), original_T_cache_SHA256=sha(old / 'private_fixed34_gene_T_counts_and_depth.npz')), indent=2) + '\n')
    shape = (3, ns * nm)
    n = np.zeros(shape, np.int64)
    rs = np.zeros((*shape, ng), np.float64)
    cs = np.zeros_like(rs)
    capture = np.zeros((*shape, 2, len(names)), np.int64)
    joint = np.zeros((*shape, 2, len(pairs)), np.int64)
    total_T = 0
    entries = 0
    for start in range(0, nc, cfg['chunk_rows']):
        stop = min(nc, start + cfg['chunk_rows'])
        ptr = np.asarray(counts['indptr'][start:stop + 1], np.int64)
        lo, hi = int(ptr[0]), int(ptr[-1])
        x = csr_matrix((np.asarray(counts['data'][lo:hi]), np.asarray(counts['indices'][lo:hi]), ptr - lo), shape=(stop - start, nf))
        x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
        assert np.isfinite(x.data).all() and (x.data > 0).all() and np.equal(x.data, np.floor(x.data)).all()
        ss = source_code[start:stop]
        assert measured[x.indices, np.repeat(ss, np.diff(x.indptr))].all()
        total = np.asarray(x.sum(axis=1)).ravel().astype(np.int64)
        assert (total > 0).all()
        raw = np.zeros((len(total), ng), np.int64)
        for si in np.unique(ss):
            ix = np.flatnonzero(ss == si)
            present = np.flatnonzero(maps[si] >= 0)
            raw[np.ix_(ix, present)] = x[ix][:, maps[si][present]].toarray()
        tt = md['major'][start:stop] == 'T'
        nt = int(tt.sum())
        if nt:
            old_raw = cache['raw_gene_counts'][total_T:total_T + nt]
            expected_raw = raw[tt].copy()
            for si in np.unique(ss[tt]):
                expected_raw[np.ix_(ss[tt] == si, maps[si] < 0)] = -1
            np.testing.assert_array_equal(expected_raw, old_raw)
            np.testing.assert_array_equal(total[tt], cache['full_native_total_counts'][total_T:total_T + nt])
            np.testing.assert_array_equal(np.diff(x.indptr)[tt], cache['full_native_detected_genes'][total_T:total_T + nt])
            total_T += nt
        positive = np.stack([np.column_stack([((raw[:, [genes.index(g) for g in panels[p]]] >= threshold).sum(axis=1) >= 2) & complete[ss, k] for k, p in enumerate(names)]) for threshold in [1, 2]], axis=1)
        for scope, mask in enumerate(masks):
            take = mask[start:stop]
            gg = group[start:stop][take]
            if not take.any():
                continue
            nn, aa, bb = aggregate(raw[take], total[take], gg, ns * nm)
            n[scope] += nn; rs[scope] += aa; cs[scope] += bb
            for ti in range(2):
                for pi in range(len(names)):
                    capture[scope, :, ti, pi] += np.bincount(gg[positive[take, ti, pi]], minlength=ns * nm)
                for qi, (left, right) in enumerate(pairs):
                    pp = positive[take, ti, names.index(left)] & positive[take, ti, names.index(right)]
                    joint[scope, :, ti, qi] += np.bincount(gg[pp], minlength=ns * nm)
        entries += hi - lo
        if start // cfg['chunk_rows'] % 16 == 0 or stop == nc:
            (out / 'private_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), rows_read=stop, total_rows=nc, stored_entries_read=entries, original_T_rows_exactly_replayed=total_T)) + '\n')
    assert total_T == len(cache['global_original_T_row_index'])
    np.testing.assert_array_equal(n[0], expected)
    for k, mask in enumerate(masks):
        np.testing.assert_array_equal(n[k], np.bincount(group[mask], minlength=ns * nm))
    descriptions = []
    overlaps = []
    for scope, scope_name in enumerate(cfg['scopes']):
        for si, source in enumerate(sources):
            sl = slice(si * nm, (si + 1) * nm)
            nums = n[scope, sl]
            for pi, panel_name in enumerate(names):
                gi = [genes.index(g) for g in panels[panel_name]]
                raw_shares = shares(rs[scope, sl][:, gi], nums) if complete[si, pi] else None
                cp_shares = shares(cs[scope, sl][:, gi], nums) if complete[si, pi] else None
                available = raw_shares is not None and cp_shares is not None
                reason = 'DESCRIPTIVE_RNA_CAPTURE_ONLY' if available else ('MISSING_OR_AMBIGUOUS_GENE' if not complete[si, pi] else 'EMPTY_SCOPE_OR_ZERO_PANEL_GENE')
                for mi, major in enumerate(majors):
                    rec = dict(source=source, scope=scope_name, major=major, panel=panel_name, original_cells=int(nums[mi]), source_scope_cells=int(nums.sum()), panel_status=reason, all_required_genes_measurable=bool(complete[si, pi]))
                    if available:
                        rec.update(raw_transcript_share=float(raw_shares[0][mi]), CP10K_cell_weighted_share=float(cp_shares[0][mi]), equal_major_mean_CP10K_profile_share=float(cp_shares[1][mi]))
                    for ti, threshold in enumerate([1, 2]):
                        rec['captured_cells_threshold' + str(threshold)] = int(capture[scope, si * nm + mi, ti, pi]) if complete[si, pi] else None
                    descriptions.append(rec)
            for qi, (left, right) in enumerate(pairs):
                ok = bool(complete[si, names.index(left)] and complete[si, names.index(right)])
                for mi, major in enumerate(majors):
                    for ti, threshold in enumerate([1, 2]):
                        overlaps.append(dict(source=source, scope=scope_name, major=major, panelA=left, panelB=right, threshold=threshold, original_cells=int(nums[mi]), both_panels_measurable=ok, both_captured_cells=int(joint[scope, si * nm + mi, ti, qi]) if ok else None))
    def save(name, value):
        (out / name).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    save('private_all_source_scope_major_panel_descriptions.json', descriptions)
    save('private_all_source_scope_major_two_arm_overlaps.json', overlaps)
    save('private_complete_numeric_proofs.json', dict(all_cells_and_original_T_native_values_replayed=True, every_chunk_grouped_raw_and_CP10K_sums_two_implementations=True, all_scope_metadata_counts_two_implementations=True, original_T_rows_replayed=total_T, all_cells_read=nc, stored_entries_read=entries, all_source_physical_maps_equal_prior=True, all_panel_shares_sum_to_one_when_available=True, missing_gene_and_empty_groups_retained=True, no_patient_or_library_rows_output=True))
    files = {p.name:dict(bytes=p.stat().st_size, SHA256=sha(p)) for p in out.iterdir() if p.is_file()}
    assert sum(v['bytes'] for v in files.values()) <= cfg['aggregate_output_cap_bytes']
    save('private_terminal.json', dict(status='COMPLETE_ALL_LINEAGE_FIXED_PANEL_RNA_SOURCE_DESCRIPTIONS_REVIEW_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), files=files, real_RNA_and_per_cell_arrays_RAM_only=True, original_project_read_only=True, new_clinical_P_or_labels=False, source_major_panel_positions=len(descriptions), source_major_pair_positions=len(overlaps)))
    print(json.dumps(dict(status='COMPLETE_REVIEW_PENDING', terminal_SHA256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
