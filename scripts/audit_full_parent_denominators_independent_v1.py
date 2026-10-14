"""Independently compare retained-sample parent counts with full node units.

Source-unit totals include libraries outside the retained panel; do not mistake
those whole-library exclusions for missing cells within a retained library.
No scores, effect estimates or historical results are changed.
"""
import csv, hashlib, json, os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import h5py
from anndata._io.specs import read_elem


def main():
    project = Path(os.environ['PANCANCER_PROJECT_ROOT'])
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    annotation = project / 'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1'
    out = root / 'full_parent_denominator_independent_audit_20261014_v1'
    assert not out.exists()
    catalogue_path = annotation / 'all_node_number_catalog.tsv'
    with catalogue_path.open() as f:
        catalogue = list(csv.DictReader(f, delimiter='\t'))
    rows = []
    for source in ['GSE161529', 'GSE182109', 'GSE183904']:
        units_path = annotation / 'by_gse' / source / 'units.tsv'
        full = Counter()
        with units_path.open() as f:
            for entry in csv.DictReader(f, delimiter='\t'):
                node = catalogue[int(entry['node_index'])]
                full[entry['sample_uid'], node['majorCluster']] += int(entry['n_cells'])
        native = root / 'state_singlecell_inputs_20261014_v1/private_native_sources' / (source + '.h5ad')
        with h5py.File(native) as f:
            obs = read_elem(f['obs'])
        selected = Counter(zip(obs.sample_uid.astype(str), obs.major.astype(str)))
        for major in ['T cell', 'Macrophage']:
            keys = [key for key in selected if key[1] == major]
            assert keys
            mismatch = [key for key in keys if full[key] != selected[key]]
            assert not mismatch, 'Retained-library parent denominator differs from full node units'
            rows.append(dict(source=source, major=major,
                retained_library_parent_keys=len(keys), all_keys_equal=True,
                source_units_sha256=hashlib.sha256(units_path.read_bytes()).hexdigest()))
    out.mkdir()
    report = dict(status='PASS_ALL_RETAINED_LIBRARY_PARENT_COUNTS_EQUAL_FULL_SOURCE_UNITS',
        finished_UTC=datetime.now(timezone.utc).isoformat(), comparisons=rows,
        interpretation='Source-wide cell totals differ because other libraries were excluded; retained-library parent denominators agree exactly.',
        scientific_results_changed=False, new_P_or_BH=False,
        catalogue_sha256=hashlib.sha256(catalogue_path.read_bytes()).hexdigest())
    (out / 'private_independent_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
