"""Bind full-native local discovery fit pools and measured model candidates.

No graph, final identity, count copy, new QC call or downstream eligibility is
created here. Existing explicit exclusions constrain fitting; archive rows stay.
"""
import argparse
import datetime
import gc
import json
from pathlib import Path
import resource
import subprocess
import sys
import time
import numpy as np

from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256
from native_annotation_review_v1 import native_measured_genes
from native_local_fit_pool_v1 import fit_pool
from summarize_full_native_annotation_evidence_v1 import read_source


def main(root):
    start = time.monotonic()
    cfg = json.loads((root / 'config/local_discovery_input_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']: raise ValueError('Bound input changed')
    route_root, ref_root = Path(cfg['routing']), Path(cfg['reference'])
    for other in [route_root, ref_root]:
        terminal = json.loads((other / 'provenance/supervisor_terminal_v1.json').read_text())
        if terminal['actual_child_exit_code'] != 0 or terminal['wall_timeout']: raise ValueError('Input producer incomplete')
    test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_native_local_fit_pool_v1.py')],
                          capture_output=True, text=True, timeout=120)
    (root / 'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode or 'Ran 4 tests' not in test.stderr or 'skipped' in test.stderr.lower(): raise ValueError('Fit-pool counterexamples failed')
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': 0, 'n_tests': 4, 'no_skips': True})
    routes = json.loads((route_root / 'outputs/full_native_provisional_routing_index_v1.json').read_text())
    refs = {r['source']: r for r in json.loads((ref_root / 'outputs/full_archive_reference_index_v1.json').read_text())}
    reader = NativeArchiveCollection(cfg['archive_manifest'])
    if {r['source'] for r in routes} != set(refs) or set(refs) != set(reader.archive_partitions): raise ValueError('Source scopes differ')
    panels = json.loads(Path(cfg['pilot_feature_spec']).read_text())['heldout_RNA_readouts']
    heldout = sorted({gene for group in panels.values() for panel in group.values() for gene in panel})
    index = []
    for item in routes:
        ds = item['source']; meta, routed = read_source(route_root, item); rm, evidence = read_source(ref_root, refs[ds])
        rows = reader.rows(ds); n = len(rows)
        np.testing.assert_array_equal(rows, routed['native_row']); np.testing.assert_array_equal(rows, evidence['native_row'])
        np.testing.assert_array_equal(routed['positive_RNA'], evidence['positive_RNA'])
        identity_sha = text_axis_sha256(reader.native_ids(ds))
        if meta['native_ids_sha256'] != identity_sha or rm['native_ids_sha256'] != identity_sha: raise ValueError('Native identity differs')
        names = set(reader.field_names(ds))
        kwargs = {argument: reader.field(ds, field) if field in names else None
                  for argument, field in [('doublet', 'residual_doublet_flag'), ('excluded', 'annotation_exclude_from_primary')]}
        result = fit_pool(routed['positive_RNA'], reader.membership(ds)['in_frozen_context'], reader.frozen_expression_units(ds), **kwargs)
        result.update(native_row=rows, provisional_family_code=routed['provisional_family_code'])
        output = root / 'outputs' / (ds + '_local_fit_pool_v1.npz')
        with output.open('xb') as f: np.savez_compressed(f, **result)
        with np.load(output, allow_pickle=False) as saved:
            for name, value in result.items(): np.testing.assert_array_equal(saved[name], value)
        genes, columns = native_measured_genes(reader, ds)
        lookup = reader.feature_lookup(ds)
        heldout_columns = {int(column) for gene in heldout for column in lookup.get(gene, set())}
        allowed = np.asarray([int(column) not in heldout_columns and gene not in heldout for gene, column in zip(genes, columns)])
        candidate_path = root / 'outputs' / (ds + '_local_feature_candidates_v1.json')
        save(candidate_path, {'source': ds, 'genes': [g for g, ok in zip(genes, allowed) if ok],
            'native_columns': columns[allowed].tolist(), 'heldout_RNA_readouts': heldout,
            'heldout_physical_alias_columns_excluded': sorted(heldout_columns),
            'feature_selection_not_performed': True, 'native_all_gene_archive_unchanged': True})
        family_names = meta['family_code_dictionary']
        branches = []
        for code in np.unique(result['provisional_family_code']):
            member = result['provisional_family_code'] == code
            branches.append({'family': family_names[code] if code >= 0 else 'unresolved_context',
                'n_archive_rows': int(member.sum()), 'n_positive_RNA': int((member & result['positive_native_RNA']).sum()),
                'n_fit_pool_candidates': int((member & result['fit_pool_candidate']).sum())})
        counts = {name: int(value.sum()) for name, value in result.items() if value.dtype == bool}
        info = {'source': ds, 'n_archive_rows': n, 'native_rows_sha256': row_axis_sha256(rows), 'native_ids_sha256': identity_sha,
            'pool_flag_counts_overlapping_not_additive': counts, 'branches': branches, 'n_native_measured_features': len(genes),
            'n_model_feature_candidates': int(allowed.sum()), 'n_heldout_alias_columns': len(heldout_columns),
            'feature_candidates': candidate_path.name, 'feature_candidates_sha256': sha256_file(candidate_path),
            'fit_pool_is_provisional_not_complete_QC_acceptance': True,
            'preintegration_release_flags_not_used_as_automatic_cell_exclusions': True,
            'membership_and_explicit_exclusions_preserved_for_unfitted_projection_review': True,
            'patient_composition_clinical_eligibility_assigned': False}
        meta_path = output.with_suffix('.json'); save(meta_path, info)
        index.append({'source': ds, 'n_archive_rows': n, 'arrays': output.name, 'arrays_sha256': sha256_file(output),
                      'metadata': meta_path.name, 'metadata_sha256': sha256_file(meta_path)})
        save(root / 'checkpoints' / (ds + '.json'), {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                                  'source': ds, 'completed_sources': len(index)})
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > cfg['output_cap_bytes']: raise ValueError('Preparation cap exceeded')
        reader._fields.clear(); del result, routed, evidence, kwargs; gc.collect()
    path = root / 'outputs/full_native_local_discovery_input_index_v1.json'; save(path, index)
    save(root / 'provenance/terminal_v1.json', {'status': 'FULL_NATIVE_LOCAL_DISCOVERY_INPUTS_PREPARED_NO_GRAPH_FITTED',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_sources': len(index),
        'index_sha256': sha256_file(path), 'wall_seconds': time.monotonic() - start,
        'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        'native_expression_read': False, 'new_authoritative_labels': False})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--output', required=True); root = Path(p.parse_args().output)
    try: main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
