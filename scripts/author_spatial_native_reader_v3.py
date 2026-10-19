"""All-native-panel and source metadata parity; arrays remain in RAM."""
import gzip
import tarfile
import numpy as np
import pandas as pd
import anndata as ad
def native_integer(values):
    s = pd.Series(values).astype(str)
    a = pd.to_numeric(s, errors='raise').to_numpy()
    assert np.isfinite(a).all() and np.equal(a, np.floor(a)).all()
    return a.astype(np.int64).astype(str)


def frame_stream(tar, member, **kwargs):
    stream = tar.extractfile(member)
    assert stream is not None
    if member.endswith('.gz'):
        stream = gzip.GzipFile(fileobj=stream)
    return pd.read_csv(stream, **kwargs)


def load_native(path, cache, archive, target_genes):
    a = ad.read_h5ad(path)
    with np.load(cache, allow_pickle=False) as z:
        meta = {k: z[k] for k in ['author_labels', 'qc', 'fov', 'cell_id', 'coords', 'area', 'native_total', 'native_ngenes', 'negatives', 'negative_genes', 'background_q95']}
    n = a.n_obs
    assert len(meta['qc']) == n and a.var_names.is_unique
    mapping = dict(a.uns['archive_member_map'])
    fv, ci = native_integer(a.obs['fov']), native_integer(a.obs['cell_ID'])
    np.testing.assert_array_equal(fv, meta['fov'].astype(str))
    np.testing.assert_array_equal(ci, meta['cell_id'].astype(str))
    lookup = {key: i for i, key in enumerate(zip(fv, ci))}
    assert len(lookup) == n
    target_columns = []
    for gene in target_genes:
        hits = [i for i, name in enumerate(a.var_names) if gene == str(name) or ('gene_symbol' in a.var and gene == str(a.var.iloc[i].gene_symbol))]
        assert len(hits) == 1
        target_columns.append(hits[0])
    observed_target = a.X[:, target_columns].toarray()
    target = np.zeros_like(observed_target, dtype=np.int32)
    seen = np.zeros(n, bool)
    negatives = np.zeros_like(meta['negatives'])
    before = (archive.stat().st_size, archive.stat().st_mtime_ns)
    positions = 0
    with tarfile.open(archive, 'r:*') as tar:
        header = frame_stream(tar, mapping['expr'], nrows=0).columns.tolist()
        raw_columns = []
        for i, name in enumerate(a.var_names):
            aliases = {str(name)}
            if 'gene_symbol' in a.var:
                aliases.add(str(a.var.iloc[i].gene_symbol))
            matches = sorted(aliases & set(header))
            assert len(matches) == 1, 'Physical native feature lacks unique literal raw column'
            raw_columns.append(matches[0])
        assert len(set(raw_columns)) == len(raw_columns) == a.n_vars
        names = ['fov', 'cell_ID'] + raw_columns + meta['negative_genes'].astype(str).tolist()
        assert len(set(names)) == len(names)
        for frame in frame_stream(tar, mapping['expr'], usecols=names, chunksize=2048):
            ff, cc = native_integer(frame.fov), native_integer(frame.cell_ID)
            keep = cc != '0'
            frame, ff, cc = frame.loc[keep], ff[keep], cc[keep]
            ids = np.asarray([lookup[key] for key in zip(ff, cc)], int)
            assert len(set(ids)) == len(ids) and not seen[ids].any()
            raw = frame[raw_columns].to_numpy(float)
            assert np.isfinite(raw).all() and (raw >= 0).all() and np.array_equal(raw, np.floor(raw)) and raw.max() < 2**31
            np.testing.assert_array_equal(a.X[ids].toarray(), raw)
            if 'counts' in a.layers:
                np.testing.assert_array_equal(a.layers['counts'][ids].toarray(), raw)
            np.testing.assert_array_equal(raw.sum(1), meta['native_total'][ids])
            np.testing.assert_array_equal((raw > 0).sum(1), meta['native_ngenes'][ids])
            target[ids] = raw[:, target_columns]
            negatives[ids] = frame[meta['negative_genes'].astype(str)].to_numpy()
            seen[ids] = True
            positions += raw.size
        assert seen.all()
        np.testing.assert_array_equal(target, observed_target)
        np.testing.assert_array_equal(negatives, meta['negatives'])
        metadata = frame_stream(tar, mapping['metadata'], usecols=['fov', 'cell_ID', 'Area', 'CenterX_global_px', 'CenterY_global_px'])
        ff, cc = native_integer(metadata.fov), native_integer(metadata.cell_ID)
        keep = cc != '0'
        metadata, ff, cc = metadata.loc[keep], ff[keep], cc[keep]
        keys = list(zip(ff, cc))
        assert len(keys) == len(set(keys))
        assert set(lookup) <= set(keys), 'Native expression cell lacks unique source metadata'
        selected = np.array([k in lookup for k in keys], bool)
        extra_fov = ff[~selected]
        extra_counts = {field: int((extra_fov == field).sum()) for field in sorted(set(ff))}
        metadata_source_rows = len(keys)
        metadata = metadata.loc[selected]
        keys = [k for k, take in zip(keys, selected) if take]
        assert len(keys) == n and set(keys) == set(lookup)
        ids = np.array([lookup[k] for k in keys])
        np.testing.assert_allclose(metadata[['CenterX_global_px', 'CenterY_global_px']].to_numpy(float) * .12028, meta['coords'][ids], atol=1e-6, rtol=0)
        np.testing.assert_array_equal(metadata.Area.to_numpy(float), meta['area'][ids])
    assert before == (archive.stat().st_size, archive.stat().st_mtime_ns)
    for field in sorted(set(fv)):
        at = fv == field
        threshold = np.quantile(negatives[at].reshape(-1), .95)
        np.testing.assert_array_equal(meta['background_q95'][at], np.repeat(threshold, at.sum()))
    meta['metadata_only_counts_per_fov'] = extra_counts
    meta['counts'] = target
    del a
    return meta, dict(native_cells=n, native_features=len(raw_columns), all_native_RNA_positions_replayed=positions, target_positions=n * len(target_columns), full_native_X_and_optional_counts_exact=True, full_panel_RNA_and_detected_exact=True, all_negative_controls_and_FOV_thresholds_exact=True, raw_metadata_coordinates_and_area_exact=True, raw_metadata_unique_nonzero_ID_rows=metadata_source_rows, raw_metadata_keys_outside_frozen_expression=metadata_source_rows-n, metadata_only_counts_per_FOV=extra_counts, frozen_expression_universe_unchanged=True, unobserved_RNA_not_imputed=True, physical_units_independently_validated=False)
