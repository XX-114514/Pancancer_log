"""Freeze native-source archival membership without copying any expression."""
import argparse
import collections
import datetime
import json
from pathlib import Path

import numpy as np
from anndata._io.specs import read_elem

from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_archive_collection_v1 import row_axis_sha256, text_axis_sha256


def create_manifest(context_path):
    reader = FullContextAtlas(context_path)
    partitions, stores, cache = [], {}, {}
    def axis(z, field):
        key = (str(z.store.path), field)
        if key not in cache:cache[key] = np.asarray(read_elem(z['obs/' + field]), dtype=object)
        return cache[key]
    for dataset in reader.partitions:
        z = reader.store(dataset)
        legacy = reader.is_legacy(dataset)
        field = 'gse_id' if legacy else reader.primary.partitions[dataset].get('dataset_field', 'dataset')
        native_sources = axis(z, field)
        rows = np.flatnonzero(native_sources == dataset)
        key = reader.legacy_config['counts_key'] if legacy else reader.primary.partitions[dataset]['counts_key']
        if len(native_sources) != z[key].attrs['shape'][0]:raise ValueError('Native obs/count axes differ')
        context = reader.rows(dataset)
        if not np.isin(context, rows).all():raise ValueError('Frozen context not in archive')
        ids = axis(z, z['obs'].attrs.get('_index', '_index'))[rows]
        if (len(set(ids)) != len(ids) or any(str(v) in {'', 'None', 'nan'} for v in ids)):
            raise ValueError('Native cell identity is missing or repeated within source')
        unit_field = 'sample_uid' if legacy else reader.partitions[dataset].get('unit_field', 'library_id')
        units = axis(z, unit_field)[rows]
        partitions.append({'dataset': dataset, 'source_field': field, 'expression_unit_field': unit_field,
                           'n_archive_rows': len(rows), 'n_frozen_context_rows': len(context),
                           'n_expression_unit_keys': len(set(units)),
                           'native_rows_sha256': row_axis_sha256(rows),
                           'native_ids_sha256': text_axis_sha256(ids),
                           'native_units_sha256': text_axis_sha256(units),
                           'count_origin': key,
                           'cell_key_contract': 'ordered pair (source_key, native_obs_index)',
                           'unit_key_contract': 'ordered pair (source_key, expression_unit_key); not a patient assertion'})
        store_key = (str(z.store.path), field)
        if store_key not in stores:
            stores[store_key] = {'path': store_key[0], 'source_field': field,
                                'source_row_counts': dict(collections.Counter(str(v) for v in native_sources))}
    for value in stores.values():
        value['outside_selected_source_scope'] = {
            k: v for k, v in value['source_row_counts'].items() if k not in reader.partitions}
    return {'schema': 'selected_source_native_archive_v1',
            'created_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'frozen_context': {'path': str(Path(context_path).resolve()), 'sha256': sha256_file(context_path)},
            'partitions': partitions, 'native_store_source_scope': list(stores.values()),
            'archive_scope': 'All native rows in selected-source backed stores; earlier upstream QC removals are not reconstructed',
            'source_alias_or_patient_independence_accepted': False,
            'new_model_or_clinical_eligibility_accepted': False,
            'new_identity_accepted': False, 'raw_counts_copied_or_changed': False,
            'standalone_Zarr': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--context-manifest', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    with Path(args.output).open('x') as handle:
        json.dump(create_manifest(args.context_manifest), handle, indent=2, allow_nan=False)
        handle.write('\n')
