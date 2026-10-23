"""Version source modality evidence without rewriting frozen analysis inputs."""
import argparse
import datetime
import json
from pathlib import Path
import subprocess
import sys

from annotation_representation_pilot_v1 import save
from full_context_atlas_v1 import sha256_file


def modality_decision(native_counts, design_note, public_suspension):
    observed = {str(k).strip().lower() for k, value in native_counts.items() if value > 0}
    known = observed & {'cell', 'nucleus'}
    public = set(public_suspension)
    design_nuclear = design_note.startswith('snRNA_explicit_')
    if (known == {'cell'} and ('nucleus' in public or design_nuclear)) or (known == {'nucleus'} and public == {'cell'}):
        return 'conflicting_modality_evidence_requires_review'
    if known == {'cell', 'nucleus'} or public == {'cell', 'nucleus'}:
        return 'mixed_cell_and_nucleus_requires_per_record_branch'
    if observed == {'nucleus'}:
        return 'snRNA_native_and_public_dataset_concordant' if public == {'nucleus'} else 'snRNA_native_explicit_source_review_pending'
    if observed == {'cell'}:
        return 'scRNA_native_explicit_source_review_pending'
    if public == {'nucleus'} or design_nuclear:
        return 'snRNA_source_evidence_native_suspension_incomplete'
    return 'modality_unresolved_not_inferred_from_platform_name'


def main(root):
    cfg = json.loads((root / 'config/modality_overlay_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:
            raise ValueError('Source modality input changed')
    original = Path(cfg['source_role_run'])
    supervision = json.loads((original / 'provenance/supervisor_terminal_v1.json').read_text())
    if supervision['actual_child_exit_code'] != 0 or supervision['wall_timeout']:
        raise ValueError('Completed source-role registry required')
    test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_source_modality_overlay_v1.py')], capture_output=True, text=True, timeout=120)
    (root / 'provenance/test_stderr.txt').write_text(test.stderr)
    if test.returncode != 0 or 'Ran 3 tests' not in test.stderr or 'skipped' in test.stderr.lower():
        raise ValueError('Source modality counterexamples failed')
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': 0, 'n_tests': 3, 'no_skips': True})
    source_path = original / 'outputs/source_roles_and_coverage_v1.json'
    data = json.loads(source_path.read_text())
    evidence = json.loads(Path(cfg['public_dataset_evidence']).read_text())
    curated = {r['source']: r for r in evidence['rows']}
    if len(curated) != len(evidence['rows']) or not set(curated) <= {r['source'] for r in data['sources']}:
        raise ValueError('Public dataset evidence source scope differs')
    output = []
    for row in data['sources']:
        ds = row['source']; counts = row.get('native_status_cell_counts', {}).get('suspension_type', {}).get('counts', {})
        if counts and sum(counts.values()) != row['n_archive_cells']:
            raise ValueError('Native suspension denominator differs')
        e = curated.get(ds)
        if e is not None and ds != 'CXG_' + e['dataset_id'].replace('-', '_'):
            raise ValueError('Curated dataset identifier differs from native source')
        public = e['suspension_type'] if e else []
        decision = modality_decision(counts, row['design_decision'].get('assay_branch', ''), public)
        output.append({'source': ds, 'n_archive_rows': row['n_archive_cells'],
            'original_source_role_and_capture_decision': row['design_decision'],
            'original_analysis_roles': row['analysis_roles'],
            'native_suspension_counts': counts, 'public_dataset_evidence': e,
            'modality_branch_evidence': decision,
            'capture_composition_patient_and_alias_qualification_not_promoted': True,
            'source_id_match_does_not_prove_downloaded_asset_version': True,
            'nuclear_vs_cell_reference_transfer_not_yet_accepted': True})
    path = root / 'outputs/source_modality_overlay_v1.json'
    save(path, {'sources': output, 'original_source_registry_sha256': sha256_file(source_path),
        'frozen_benchmark_inputs_unchanged': True, 'new_authoritative_labels': False,
        'native_expression_read': False, 'new_downloads': False,
        'per_record_anatomy_required_for_multiorgan_source': True,
        'adult_reference_not_automatically_applicable_to_pediatric_source': True})
    if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > cfg['output_cap_bytes']:
        raise ValueError('Source modality metadata output cap exceeded')
    save(root / 'provenance/terminal_v1.json', {'status': 'SOURCE_MODALITY_EVIDENCE_OVERLAY_COMPLETE_NOT_MAPPING_ACCEPTANCE',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_sources': len(output),
        'output_sha256': sha256_file(path), 'native_expression_read': False, 'new_authoritative_labels': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    root = Path(parser.parse_args().output)
    try:
        main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
