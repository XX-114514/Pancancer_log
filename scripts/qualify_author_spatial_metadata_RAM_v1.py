"""Author label/coordinate qualification only. Individual arrays never leave RAM.

The caller supplies safe_label from the frozen biological vocabulary gate.
Unknown strings are not emitted; annotation names are not identity truth.
"""
import io
import json
from collections import Counter
import hashlib
import numpy as np

FIELDS = ('anno', 'anno2', 'bin50 structural annotation')
EXTRA_BIOLOGICAL_LABELS = frozenset(('MAIT', 'gdT', 'gammadelta_T', 'iTLSs', 'PCA', 'IEstroma'))


def disclosed(label, genes):
    return not has_identifying_token(label) and (safe_label(label, genes) or label in EXTRA_BIOLOGICAL_LABELS)


def text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def read_category(obj, n):
    import h5py
    if not isinstance(obj, h5py.Group) or not {'categories', 'codes'} <= set(obj.keys()):
        raise ValueError('EXPECTED_STRING_CATEGORICAL_GROUP')
    cats = obj['categories']
    assert isinstance(cats, h5py.Dataset) and cats.dtype.kind in 'OSU' and len(cats) <= 256
    labels = [text(v) for v in cats[:]]
    assert len(set(labels)) == len(labels) and all(len(v) <= 512 for v in labels)
    assert obj['codes'].shape == (n,)
    codes = np.empty(n, dtype=np.int16)
    for lo in range(0, n, 4096):
        raw = np.asarray(obj['codes'][lo:lo+4096])
        assert raw.dtype.kind in 'iu' and ((raw >= -1) & (raw < len(labels))).all()
        codes[lo:lo+len(raw)] = raw
    return labels, codes


def qualify(h):
    import h5py
    obs = h['obs']; key = text(obs.attrs.get('_index', '_index'))
    n = len(obs[key])  # Shape only; never read observation identifiers.
    assert 0 < n <= 2000000
    var = h['var']; vk = text(var.attrs.get('_index', '_index'))
    assert isinstance(var[vk], h5py.Dataset) and len(var[vk]) <= 100000
    genes = [text(v).strip().upper() for v in var[vk][:]]
    assert len(set(genes)) == len(genes)
    genes = set(genes)
    categories = {}; codes = {}; summary = {}
    for field in FIELDS:
        labels, values = read_category(obs[field], n)
        categories[field] = labels; codes[field] = values
        count = np.bincount(values[values >= 0], minlength=len(labels))
        safe = {v: int(c) for v, c in zip(labels, count) if disclosed(v, genes)}
        summary[field] = dict(qualified_counts=safe,
            unqualified_rows=int(sum(c for v, c in zip(labels, count) if not disclosed(v, genes))),
            missing_rows=int((values < 0).sum()),
            unqualified_categories=sum(not disclosed(v, genes) for v in labels))
        assert sum(safe.values()) + summary[field]['unqualified_rows'] + summary[field]['missing_rows'] == n
    ds = h['obsm']['spatial']
    assert isinstance(ds, h5py.Dataset) and ds.ndim == 2 and ds.shape[0] == n and 2 <= ds.shape[1] <= 3
    assert ds.dtype.kind in 'fiu'
    xy = np.empty((n, 2), np.float64)
    for lo in range(0, n, 4096):
        xy[lo:lo+4096] = ds[lo:lo+4096, :2]
    finite = np.isfinite(xy).all(axis=1)
    distinct, inverse, counts = np.unique(xy[finite], axis=0, return_inverse=True, return_counts=True)
    usable = finite.copy(); usable[finite] = counts[inverse] == 1
    # Independent tuple counter checks the complete exact-coordinate exclusion.
    replay_counts = Counter(map(tuple, xy[finite]))
    replay_usable = np.zeros(n, dtype=bool)
    for position in np.flatnonzero(finite): replay_usable[position] = replay_counts[tuple(xy[position])] == 1
    assert np.array_equal(usable, replay_usable)
    del replay_counts, replay_usable, distinct, inverse, counts
    fine = categories['anno2']; fc = codes['anno2']; coarse = categories['anno']; cc = codes['anno']
    t_ids = [i for i, label in enumerate(coarse) if label.lower().replace('_', '').replace(' ', '') in {'tcell', 'tcells'}]
    coarse_t = np.isin(cc, t_ids)
    cd4_ids = [i for i, label in enumerate(fine) if label.startswith('CD4_') and disclosed(label, genes)]
    cd8_ids = [i for i, label in enumerate(fine) if label.startswith('CD8_') and disclosed(label, genes)]
    cd4 = np.isin(fc, cd4_ids); cd8 = np.isin(fc, cd8_ids)
    tid = fine.index('CD4_Treg') if 'CD4_Treg' in fine else -2
    treg = (fc == tid) & coarse_t & usable
    qualified_cd4 = cd4 & coarse_t & usable; qualified_cd8 = cd8 & coarse_t & usable
    concordant = {label: int(((fc == i) & coarse_t & usable).sum())
                  for i, label in enumerate(fine) if i in cd4_ids + cd8_ids}
    replay = Counter(fine[int(fc[i])] for i in np.flatnonzero(coarse_t & usable)
                     if int(fc[i]) in cd4_ids + cd8_ids)
    assert all(replay[label] == count for label, count in concordant.items())
    # Each individual vector is hashed and retained only in source RAM, never serialized.
    digest = hashlib.sha256()
    for arr in (xy, cc, fc, codes[FIELDS[2]], usable): digest.update(np.ascontiguousarray(arr).tobytes())
    payload = dict(xy=xy, usable=usable, categories=categories, codes=codes,
                   coarse_t=coarse_t, cd4=qualified_cd4, cd8=qualified_cd8, treg=treg)
    result = dict(rows=n, coordinates_finite=int(finite.sum()), duplicate_xy_rows_excluded=int(finite.sum()-usable.sum()),
                  coordinates_usable=int(usable.sum()), fields=summary,
                  fine_T_not_coarse_T=int(((cd4 | cd8) & ~coarse_t).sum()),
                  coarse_T_not_typed_CD4_CD8=int((coarse_t & ~(cd4 | cd8)).sum()),
                  coarse_T=int(coarse_t.sum()), usable_CD4=int(qualified_cd4.sum()), usable_CD8=int(qualified_cd8.sum()),
                  usable_CD4_Treg=int(treg.sum()), usable_CD4_nonTreg=int((qualified_cd4 & ~treg).sum()),
                  concordant_fine_T_counts=concordant, source_RAM_array_SHA256=digest.hexdigest())
    assert result['usable_CD4'] == result['usable_CD4_Treg'] + result['usable_CD4_nonTreg']
    return payload, result


def aggregate(results):
    totals = {k: sum(v[k] for v in results) for k in ('rows', 'coordinates_finite', 'duplicate_xy_rows_excluded',
              'coordinates_usable', 'fine_T_not_coarse_T', 'coarse_T_not_typed_CD4_CD8', 'coarse_T',
              'usable_CD4', 'usable_CD8', 'usable_CD4_Treg', 'usable_CD4_nonTreg')}
    totals['images'] = len(results); totals['fields'] = {}
    for field in FIELDS:
        count = Counter(); coverage = Counter()
        for result in results:
            for label, value in result['fields'][field]['qualified_counts'].items():
                count[label] += value; coverage[label] += int(value > 0)
        totals['fields'][field] = dict(qualified_counts=dict(sorted(count.items())),
            images_with_label=dict(sorted(coverage.items())),
            unqualified_rows=sum(v['fields'][field]['unqualified_rows'] for v in results),
            missing_rows=sum(v['fields'][field]['missing_rows'] for v in results))
        f = totals['fields'][field]
        assert sum(f['qualified_counts'].values())+f['unqualified_rows']+f['missing_rows'] == totals['rows']
    counts = Counter(); coverage = Counter(); opportunity = Counter()
    for v in results:
        for label, count in v['concordant_fine_T_counts'].items():
            counts[label] += count; coverage[label] += int(count > 0)
            # Prospective coverage only; no distances or spatial effects computed.
            if label.startswith('CD8_'):
                opportunity[label] += int(v['usable_CD4_Treg'] >= 20 and v['usable_CD4_nonTreg'] >= 20
                                          and count >= 20 and v['usable_CD8']-count >= 20)
    totals.update(concordant_fine_T_counts=dict(sorted(counts.items())),
                  images_with_concordant_fine_T=dict(sorted(coverage.items())),
                  images_meeting_future_20each_fourgroup_gate=dict(sorted(opportunity.items())))
    totals.update(source='GSE328481', all_arrays_source_RAM_only=True, no_expression_or_identifier_values_read=True,
        no_spatial_effect_or_clinical_model=True, author_annotation_not_ground_truth=True,
        source_partial_read_does_not_reverify_whole_file_SHA=True,
        unknown_units_patient_phase_endpoint_not_imputed=True)
    totals['combined_source_RAM_array_SHA256'] = hashlib.sha256(''.join(v['source_RAM_array_SHA256'] for v in results).encode()).hexdigest()
    totals['all_exact_coordinate_and_fine_parent_independent_replays_accepted'] = True
    return totals


def synthetic():
    import h5py
    b = io.BytesIO()
    with h5py.File(b, 'w') as h:
        obs = h.create_group('obs'); obs.create_dataset('_index', shape=(7,), dtype='i8')
        var = h.create_group('var'); var.create_dataset('_index', data=np.array([b'FOXP3', b'LAG3']))
        def cat(name, labels, values):
            obj = obs.create_group(name); obj.create_dataset('categories', data=np.array(labels, dtype='S'))
            obj.create_dataset('codes', data=np.array(values, np.int8))
        cat('anno', ['Tcell', 'cancer_cell'], [0,0,0,0,1,0,0])
        cat('anno2', ['CD4_Treg', 'CD4_IL7R', 'CD8_Tex', 'CD4_P001'], [0,1,2,2,0,3,-1])
        cat('bin50 structural annotation', ['Tumor', 'private_subject'], [0,0,0,1,0,0,0])
        obsm=h.create_group('obsm'); obsm.create_dataset('spatial', data=np.array([[0,0],[1,0],[2,0],[2,0],[4,0],[5,0],[np.nan,0.]]))
        payload, result = qualify(h); out = aggregate([result])
        assert result['rows']==7 and result['coordinates_usable']==4 and result['duplicate_xy_rows_excluded']==2
        assert result['usable_CD4_Treg']==1 and result['fine_T_not_coarse_T']==1
        assert result['usable_CD8']==0 and result['usable_CD4_nonTreg']==1
        assert not any(v in json.dumps(out) for v in ('CD4_P001', 'private_subject'))
        assert payload['treg'].sum()==1
    for label in EXTRA_BIOLOGICAL_LABELS: assert disclosed(label, set())
    assert not disclosed('Treg-P001', {'Treg-P001'})


if __name__ == '__main__':
    from biological_annotation_label_gate_v2 import safe_label, has_identifying_token
    synthetic(); print('PASS_METADATA_GEOMETRY_HIERARCHY_AND_NONEMISSION')
