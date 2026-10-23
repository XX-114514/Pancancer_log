"""Cross-tab audited differentiation cues with native anatomical metadata.

Uses an explicit vocabulary of anatomical field values. Unreviewed free text
stays private and unclassified; cancer origin is not inferred to be biopsy site.
"""
import argparse
import collections
import datetime
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np

from annotation_representation_pilot_v1 import ROUTES, save
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256


FIELDS = ['tissue', 'tissue_site', 'biopsy_site', 'biopsysite', 'biopsysites', 'organ', 'anatomical_site']
ANATOMY = {
    'liver': 'liver', 'lung': 'lung', 'breast': 'breast', 'mammary gland': 'breast',
    'brain': 'brain', 'cerebral cortex': 'brain', 'cerebellum': 'brain',
    'kidney': 'kidney', 'colon': 'colon', 'rectum': 'rectum', 'small intestine': 'small_intestine',
    'pancreas': 'pancreas', 'skin': 'skin', 'bone marrow': 'bone_marrow',
    'lymph node': 'lymph_node', 'lymph nodes': 'lymph_node', 'blood': 'blood', 'peripheral blood': 'blood',
    'ovary': 'ovary', 'prostate': 'prostate', 'thyroid gland': 'thyroid', 'thyroid': 'thyroid',
    'adrenal gland': 'adrenal', 'spleen': 'spleen', 'tonsil': 'tonsil', 'peritoneum': 'peritoneum',
    'omentum': 'omentum', 'pleura': 'pleura', 'bone': 'bone', 'soft tissue': 'soft_tissue',
    'cervix': 'cervix', 'uterus': 'uterus', 'esophagus': 'esophagus', 'stomach': 'stomach',
    'nasopharynx': 'nasopharynx', 'oral cavity': 'oral_cavity', 'meninges': 'meninges'}
MISSING = {'', 'na', 'n/a', 'nan', 'none', '<na>', 'unknown', 'not provided', 'not available'}


def anatomy_values(values):
    result = []
    for raw in values:
        value = str(raw).strip().lower()
        result.append('missing_or_explicit_unknown' if value in MISSING else ANATOMY.get(value, 'unreviewed_value'))
    return np.asarray(result)


def context_counts(categories, support, panels):
    if support.shape != (len(categories), len(panels)):
        raise ValueError('Context and supplemental cue axes differ')
    result = []
    for label in np.unique(categories):
        mask = categories == label
        result.append({'context': str(label), 'n_archive_rows': int(mask.sum()),
            'cue_positive_counts': {name: int(support[mask, j].sum()) for j, name in enumerate(panels)}})
    return result


def main(root):
    start = time.monotonic()
    cfg = json.loads((root / 'config/context_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:
            raise ValueError('Frozen context input changed')
    for name in ['producer', 'audit']:
        status = json.loads((Path(cfg[name]) / 'provenance/supervisor_terminal_v1.json').read_text())
        if status['actual_child_exit_code'] != 0 or status['wall_timeout']:
            raise ValueError('Supplemental producer and independent audit must succeed')
    producer = Path(cfg['producer']); audit = Path(cfg['audit'])
    terminal = json.loads((producer / 'provenance/terminal_v1.json').read_text())
    checked = json.loads((audit / 'provenance/terminal_v1.json').read_text())
    if checked['producer_terminal_sha256'] != sha256_file(producer / 'provenance/terminal_v1.json'):
        raise ValueError('Independent audit is not bound to producer')
    test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_supplemental_native_context_v1.py')],
        capture_output=True, text=True, timeout=120)
    (root / 'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode != 0 or 'Ran 3 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise ValueError('Anatomical metadata counterexamples failed')
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': 0, 'n_tests': 3, 'no_skips': True})
    index_path = producer / 'outputs/supplemental_native_index_v1.json'
    if sha256_file(index_path) != terminal['index_sha256']:
        raise ValueError('Accepted supplemental index changed')
    index = json.loads(index_path.read_text()); reader = NativeArchiveCollection(cfg['archive_manifest'])
    if len(index) != len(reader.archive_partitions) or {r['source'] for r in index} != set(reader.archive_partitions):
        raise ValueError('Native context source scope differs')
    reports = []
    for record in index:
        ds = record['source']; rows = reader.rows(ds); n = len(rows)
        for field in ['arrays', 'metadata']:
            if sha256_file(producer / 'outputs' / record[field]) != record[field + '_sha256']:
                raise ValueError('Accepted supplemental artifact changed')
        meta = json.loads((producer / 'outputs' / record['metadata']).read_text())
        with np.load(producer / 'outputs' / record['arrays'], allow_pickle=False) as stored:
            np.testing.assert_array_equal(rows, stored['native_row'])
            support = stored['cue_support']
        if row_axis_sha256(rows) != meta['native_rows_sha256'] or text_axis_sha256(reader.native_ids(ds)) != meta['native_ids_sha256']:
            raise ValueError('Native identity axis differs')
        panels = meta['panel_order']; names = set(reader.field_names(ds))
        major = np.asarray([ROUTES.get(str(x), 'unresolved_context') for x in reader.frozen_labels(ds)['major']])
        report = {'source': ds, 'n_archive_rows': n,
            'inherited_major_context_not_new_identity': context_counts(major, support, panels),
            'anatomical_fields': {}, 'cancer_origin_used_as_biopsy_site': False,
            'organ_reference_eligibility_assigned': False}
        for field in FIELDS:
            if field not in names:
                report['anatomical_fields'][field] = {'status': 'field_absent_not_anatomy_absence'}
                continue
            categories = anatomy_values(reader.field(ds, field))
            counts = context_counts(categories, support, panels)
            if sum(c['n_archive_rows'] for c in counts) != n:
                raise ValueError('Anatomical context denominator differs')
            for j, panel in enumerate(panels):
                if sum(c['cue_positive_counts'][panel] for c in counts) != int(support[:, j].sum()):
                    raise ValueError('Anatomical cue cross-tab denominator differs')
            report['anatomical_fields'][field] = {'status': 'explicit_vocabulary_metadata_context_not_verified_biopsy_site',
                'counts': counts, 'unreviewed_free_text_not_exported': True}
        reports.append(report)
        save(root / 'checkpoints' / (ds + '.json'), {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'source': ds, 'completed_sources': len(reports)})
        reader._fields.clear()
    path = root / 'outputs/supplemental_native_context_v1.json'
    save(path, {'sources': reports, 'new_authoritative_labels': False, 'native_expression_read': False,
        'no_cue_or_tissue_metadata_automatically_certifies_identity_or_reference_scope': True,
        'individual_records_or_unreviewed_free_text_exported': False})
    if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > cfg['output_cap_bytes']:
        raise ValueError('Supplemental context output cap exceeded')
    save(root / 'provenance/terminal_v1.json', {'status': 'AUDITED_SUPPLEMENTAL_CUES_AND_NATIVE_CONTEXT_CROSSTABS_COMPLETE',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_sources': len(reports),
        'output_sha256': sha256_file(path), 'native_expression_read': False, 'new_authoritative_labels': False,
        'wall_seconds': time.monotonic() - start, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    root = Path(parser.parse_args().output)
    try:
        main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
