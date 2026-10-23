"""Qualify full-archive measured feature intersections without reading counts.

Feature coverage does not certify biological reference applicability. The gene
lists and source coverage stay private; no query gene absence becomes a zero.
"""
import argparse
import datetime
import json
from pathlib import Path
import resource
import time

from celltypist.models import Model
from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection
from native_annotation_review_v1 import native_measured_genes


def main(root):
    start = time.monotonic()
    cfg = json.loads((root / 'config/feature_coverage_spec_v1.json').read_text())
    for record in cfg['bindings'].values():
        if sha256_file(record['path']) != record['sha256']:
            raise ValueError('Frozen feature coverage input changed')
    reader = NativeArchiveCollection(cfg['archive_manifest'])
    old = json.loads(Path(cfg['pilot_feature_plan']).read_text())
    heldout = {g for panels in old['heldout_RNA_readouts'].values() for panel in panels.values() for g in panel}
    shared = json.loads(Path(cfg['pilot_selected_features']).read_text())['selected_feature_names']
    references = {}
    for name, entry in cfg['models'].items():
        if sha256_file(entry['path']) != entry['sha256']:
            raise ValueError('Cached reference model changed')
        references[name] = set(Model.load(entry['path']).features)
    common, union, records = None, set(), []
    for source in reader.archive_partitions:
        genes, columns = native_measured_genes(reader, source)
        if len(genes) != len(set(genes)) or len(columns) != len(set(columns)):
            raise ValueError('Canonical features do not map to unique physical columns')
        measured = set(genes)
        common = measured.copy() if common is None else common & measured
        union |= measured
        records.append({'source': source, 'n_unique_measured_canonical_genes': len(measured),
            'n_original_pilot_feature_candidates_measured': len(measured & set(old['common_model_feature_candidates'])),
            'n_original_selected_shared_features_measured': len(measured & set(shared)),
            'original_shared_feature_count': len(shared),
            'original_shared_reference_full_feature_axis_supported': set(shared) <= measured,
            'n_predeclared_readout_genes_measured': len(measured & heldout),
            'reference_coverage': {name: {'n_features': len(features), 'n_measured': len(measured & features),
                'fraction': len(measured & features) / len(features),
                'coverage_at_least_existing_threshold': len(measured & features) / len(features) >= cfg['reference_feature_fraction'],
                'biological_scope_qualified': False} for name, features in references.items()}})
        save(root / 'checkpoints' / (source + '.json'), {'source': source,
            'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'completed_sources': len(records)})
    common_nonheldout = sorted(common - heldout)
    payload = {'source_coverage': records, 'n_sources': len(records), 'union_canonical_gene_count': len(union),
        'common_measured_canonical_gene_count': len(common), 'common_nonheldout_feature_candidates': common_nonheldout,
        'heldout_readout_genes': sorted(heldout), 'counts_read': False, 'feature_selection_fitted': False,
        'source_model_eligibility_not_assigned': True, 'reference_scope_not_certified_by_feature_coverage': True,
        'no_gene_absence_filled_as_zero': True}
    output = root / 'outputs/full_archive_feature_coverage_v1.json'; save(output, payload)
    size = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
    if size > cfg['output_cap_bytes']:
        raise ValueError('Metadata output cap exceeded')
    save(root / 'provenance/terminal_v1.json', {'status': 'FULL_NATIVE_MEASURED_FEATURE_INTERSECTION_AND_REFERENCE_COVERAGE_COMPLETE',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'output_sha256': sha256_file(output),
        'n_sources': len(records), 'n_common_nonheldout_candidates': len(common_nonheldout),
        'n_sources_supporting_original_shared_feature_axis': sum(r['original_shared_reference_full_feature_axis_supported'] for r in records),
        'wall_seconds': time.monotonic() - start, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        'counts_read_or_copied': False, 'new_authoritative_labels': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    root = Path(parser.parse_args().output)
    try:
        main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
