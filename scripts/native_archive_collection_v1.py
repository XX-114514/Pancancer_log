"""Selected-source native archive, distinct from the frozen analysis context.

This reader addresses the explicitly bound source stores, not upstream raw
droplets removed by earlier processing. It keeps all rows belonging to selected
sources in those stores, including analysis/QC exclusions. No native input is
written, no new identity is inferred, and no model eligibility is implied.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from anndata._io.specs import read_elem

from full_context_atlas_v1 import FullContextAtlas, sha256_file


def text_axis_sha256(values):
    """Length-framed UTF-8, insensitive to numpy object memory addresses."""
    digest = hashlib.sha256()
    for value in values:
        text = str(value).encode('utf-8')
        digest.update(len(text).to_bytes(8, 'little'));digest.update(text)
    return digest.hexdigest()


def row_axis_sha256(values):
    return hashlib.sha256(np.asarray(values, dtype='<i8').tobytes()).hexdigest()


class NativeArchiveCollection(FullContextAtlas):
    def __init__(self, manifest):
        path = Path(manifest)
        spec = json.loads(path.read_text())
        if spec.get('schema') != 'selected_source_native_archive_v1':
            raise ValueError('Unsupported native archive schema')
        bound = spec['frozen_context']
        if sha256_file(bound['path']) != bound['sha256']:
            raise ValueError('Frozen analysis context changed')
        super().__init__(bound['path'])
        self.archive_manifest_path, self.archive_manifest = path, spec
        self.archive_partitions = {p['dataset']: p for p in spec['partitions']}
        if (len(self.archive_partitions) != len(spec['partitions']) or
                set(self.archive_partitions) != set(self.partitions)):
            raise ValueError('Archive must preserve the frozen selected-source scope')
        self._archive_rows, self._native_arrays = {}, {}
        self._roster_index = np.full(self.legacy[self.legacy_config['counts_key']].attrs['shape'][0], -1, np.int64)
        self._roster_index[self._positions] = np.arange(len(self._positions))
        meta = pd.read_csv(self.legacy_config['sample_metadata_path'], sep='\t', dtype=str,
                           keep_default_na=False, usecols=['sampleID', 'gse_id', 'eligible_min100_cells'])
        if not meta.sampleID.is_unique or not set(meta.eligible_min100_cells) <= {'True', 'False'}:
            raise ValueError('Invalid frozen unit metadata')
        if not set(self._sample_categories) <= set(meta.sampleID):
            raise ValueError('Frozen roster has an unregistered expression unit')
        meta = meta.set_index('sampleID').loc[self._sample_categories]
        self._all_roster_sources = meta.gse_id.to_numpy()[self._sample_codes]
        self._all_roster_unit_eligible = meta.eligible_min100_cells.eq('True').to_numpy()[self._sample_codes]

    def _native_field(self, dataset, name):
        z = self.store(dataset)
        key = (str(z.store.path), name)
        if key not in self._native_arrays:
            if name not in z['obs']:raise KeyError(name)
            self._native_arrays[key] = np.asarray(read_elem(z['obs/' + name]), dtype=object)
        return self._native_arrays[key]

    def context_rows(self, dataset):
        """Original frozen selection/order; not expanded for the benchmark."""
        return FullContextAtlas.rows(self, dataset)

    def rows(self, dataset):
        if dataset not in self._archive_rows:
            entry = self.archive_partitions[dataset]
            expected_field = 'gse_id' if self.is_legacy(dataset) else self.primary.partitions[dataset].get('dataset_field', 'dataset')
            if entry['source_field'] != expected_field:
                raise ValueError('Archive selector must use native source identity')
            rows = np.flatnonzero(self._native_field(dataset, expected_field) == dataset)
            if len(rows) != entry['n_archive_rows'] or row_axis_sha256(rows) != entry['native_rows_sha256']:
                raise ValueError('Native archive row membership changed')
            # A partition cannot erase existing context rows or silently change
            # their order in the historical reader.
            if not np.isin(self.context_rows(dataset), rows).all():
                raise ValueError('Frozen context outside native archive')
            self._archive_rows[dataset] = rows
        return self._archive_rows[dataset]

    def native_ids(self, dataset):
        z = self.store(dataset)
        field = z['obs'].attrs.get('_index', '_index')
        values = self._native_field(dataset, field)[self.rows(dataset)]
        if text_axis_sha256(values) != self.archive_partitions[dataset]['native_ids_sha256']:
            raise ValueError('Native archive cell identity changed')
        return values

    def frozen_expression_units(self, dataset):
        """Native keys verified against all roster rows, including small units."""
        field = self.archive_partitions[dataset]['expression_unit_field']
        expected = 'sample_uid' if self.is_legacy(dataset) else self.partitions[dataset].get('unit_field', 'library_id')
        if field != expected:raise ValueError('Expression-unit namespace changed')
        units = self._native_field(dataset, field)[self.rows(dataset)]
        if text_axis_sha256(units) != self.archive_partitions[dataset]['native_units_sha256']:
            raise ValueError('Native archive unit identity changed')
        if self.is_legacy(dataset):
            index = self._roster_index[self.rows(dataset)]
            found = index >= 0
            expected_units = self._sample_categories[self._sample_codes[index[found]]]
            if (not np.array_equal(units[found], expected_units) or
                    not np.all(self._all_roster_sources[index[found]] == dataset)):
                raise ValueError('Native/frozen unit or source namespace mismatch')
        return units

    def membership(self, dataset):
        """Mutually exclusive context/low-unit/outside-roster provenance.

        False for in_legacy_roster in a newer source means not applicable;
        legacy_roster_applicable is provided so it is never an exclusion call.
        """
        rows = self.rows(dataset)
        context = np.isin(rows, self.context_rows(dataset))
        legacy = self.is_legacy(dataset)
        roster = self._roster_index[rows] if legacy else np.full(len(rows), -1)
        in_roster = roster >= 0
        low_unit = np.zeros(len(rows), bool)
        if legacy:
            low_unit[in_roster] = ~self._all_roster_unit_eligible[roster[in_roster]]
            if not np.array_equal(context, in_roster & ~low_unit):
                raise ValueError('Archive/context/roster qualification disagrees')
        status = np.full(len(rows), 'native_source_row_outside_frozen_context', dtype=object)
        status[context] = 'frozen_context_member'
        if legacy:
            status[low_unit] = 'frozen_roster_below_expression_unit_min100'
            status[~in_roster] = 'native_row_outside_frozen_roster'
        return {'in_frozen_context': context,
                'legacy_roster_applicable': np.full(len(rows), legacy, bool),
                'in_legacy_roster': in_roster,
                'below_legacy_expression_unit_gate': low_unit,
                'membership_status': status,
                'new_model_eligibility': 'unassessed_not_inferred_from_membership'}

    def frozen_labels(self, dataset):
        if not self.is_legacy(dataset):
            return FullContextAtlas.frozen_labels(self, dataset)
        index = self._roster_index[self.rows(dataset)]
        fine = np.full(len(index), 'unresolved_outside_frozen_roster', dtype=object)
        found = index >= 0
        fine[found] = self._label_categories[self._label_codes[index[found]]]
        return {'major': self.field(dataset, 'majorCluster'), 'fine': fine,
                'fine_available_from_frozen_roster': found,
                'status': 'inherited_identity_not_revalidated_outside_roster_unresolved'}

    def author_labels(self, dataset):
        if self.is_legacy(dataset):
            field = 'ev_author_label'
        else:
            annotation = self.primary.partitions[dataset]['annotation']
            field = annotation.get('field', 'source_author_label')
        if field not in self.field_names(dataset):
            return {'values': np.full(len(self.rows(dataset)), '', dtype=object),
                    'status': 'author_label_unavailable', 'field': None}
        return {'values': self.field(dataset, field), 'field': field,
                'status': 'carried_author_evidence_not_independent_truth'}

    def read_panel(self, dataset, positions, genes, layer='native_counts'):
        result = super().read_panel(dataset, positions, genes)
        if layer != 'native_counts':
            raise ValueError('Archive panels read native counts; applied correction remains a separately named source layer')
        result['view'] = 'selected_source_native_archive'
        result['archive_membership_does_not_imply_model_or_clinical_eligibility'] = True
        return result
