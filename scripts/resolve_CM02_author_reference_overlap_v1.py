"""Resolve overlapping author reference files without adopting cell identities.

All cell-level keys and labels stay in RAM. Export source/node/phase summaries
with explicit unresolved overlap and unannotated denominators only.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import json
import os

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    parent = Path(os.environ['PANCANCER_PARENT_OUTPUT'])
    audit = Path(os.environ['PANCANCER_BRIDGE_AUDIT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    terminal = json.loads((parent / 'private_terminal.json').read_text())
    accepted = json.loads(audit.read_text())
    assert accepted['all18_join_positions_and_source_crosscounts_replayed']
    assert accepted['parent_terminal_SHA256'] == sha(parent / 'private_terminal.json')
    source = parent.parent / 'probe_CM02_author_library_barcode_bridge_v1.py'
    assert sha(source) == terminal['source_SHA256']
    spec = importlib.util.spec_from_file_location('pinned_reader', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    tab = pd.read_csv(root / 'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv', sep='\t', dtype=str, keep_default_na=False,
                      usecols=['library_id'] + module.PATIENT_FIELDS + module.TIME_FIELDS)
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    z = zarr.open_group(str(atlas / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    mask = np.asarray(read_elem(z['obs/dataset'])).astype(str) == 'GSE286827'
    cur = pd.DataFrame({k: np.asarray(read_elem(z['obs'][k])).astype(str)[mask]
                        for k in ['library_id', 'timepoint', 'major', 'node_id', 'original_cell_id']})
    cur['position'] = np.arange(len(cur))
    cur['barcode'] = cur.original_cell_id.map(lambda v: module.barcode(v, 'within_reviewed_library_16nt'))
    assert cur.barcode.notna().all() and not cur.duplicated(['library_id', 'barcode']).any()
    opportunities = json.loads((parent / 'private_all18_library_barcode_opportunities.json').read_text())
    cross = json.loads((parent / 'private_all_source_node_state_crosscounts.json').read_text())
    mapped = []
    for name, field, expected in module.SOURCES[1:]:
        opportunity = next(r for r in opportunities if r['annotation_source'] == name)
        a = module.read_author(name, field, expected)
        lookup = {(module.norm(r[opportunity['patient_alias_field']]), module.norm(r[opportunity['time_alias_field']])): r['library_id']
                  for r in tab.to_dict('records')}
        a['library_id'] = [lookup.get((module.norm(p), module.norm(t))) for p, t in zip(a.patient, a.phase)]
        a['barcode'] = a.cell_key.map(lambda v: module.barcode(v, 'within_reviewed_library_16nt'))
        assert a.library_id.notna().all() and a.barcode.notna().all() and not a.duplicated(['library_id', 'barcode']).any()
        m = cur.merge(a[['library_id', 'barcode', 'label']], on=['library_id', 'barcode'], validate='one_to_one')
        got = {(p, g, n, l): int(v) for (p, g, n, l), v in m.groupby(['timepoint', 'major', 'node_id', 'label']).size().items()}
        old = {(r['phase'], r['current_major'], r['current_node'], r['author_label']): r['matched_cells']
               for r in cross if r['annotation_source'] == name and r['barcode_mode'] == 'within_reviewed_library_16nt' and r['current_cell_key_field'] == 'original_cell_id'}
        assert got == old
        m['annotation_source'] = name
        mapped.append(m)
    shared = set(mapped[0].position) & set(mapped[1].position)
    assert len(shared) == accepted['CD4_CD8_shared_current_cells']
    states = pd.concat(mapped, ignore_index=True)
    unambiguous = states[~states.position.isin(shared)].copy()
    assert unambiguous.position.is_unique
    overlaps = mapped[0].merge(mapped[1][['position', 'label']], on='position', suffixes=('_CD4', '_CD8'), validate='one_to_one')
    assert len(overlaps) == len(shared)
    denominators = []
    for (phase, major, node), part in cur.groupby(['timepoint', 'major', 'node_id']):
        positions = set(part.position)
        matches = positions & set(unambiguous.position)
        unresolved = positions & shared
        assert len(part) == len(matches) + len(unresolved) + len(positions - matches - unresolved)
        denominators.append(dict(phase=phase, current_major=major, current_node=node, all_current_cells=len(part),
                                 unique_fine_reference_cells=len(matches), overlapping_reference_cells=len(unresolved),
                                 no_fine_reference_cells=len(positions - matches - unresolved)))
    keys = ['timepoint', 'major', 'node_id', 'annotation_source', 'label']
    grouped = unambiguous.groupby(keys).size().to_dict()
    independent = Counter(unambiguous[keys].itertuples(index=False, name=None))
    assert grouped == independent and sum(grouped.values()) == len(unambiguous)
    rows = [dict(phase=p, current_major=g, current_node=n, annotation_source=s, author_label=l, matched_cells=int(v))
            for (p, g, n, s, l), v in sorted(grouped.items())]
    overlap_rows = [dict(phase=p, current_major=g, current_node=n, CD4_author_label=a, CD8_author_label=b, matched_cells=int(v))
                    for (p, g, n, a, b), v in overlaps.groupby(['timepoint', 'major', 'node_id', 'label_CD4', 'label_CD8']).size().items()]
    out.mkdir(mode=0o700)
    for name, value in [('private_all_node_phase_reference_denominators.json', denominators), ('private_unique_fine_state_crosscounts.json', rows),
                        ('private_all_overlapping_reference_label_pairs.json', overlap_rows)]:
        (out / name).write_text(json.dumps(value, indent=2) + '\n')
    result = dict(status='COMPLETE_SOURCE_REFERENCE_OVERLAP_EXPLICITLY_UNRESOLVED', completed_UTC=datetime.now(timezone.utc).isoformat(),
                  source_SHA256=sha(__file__), bridge_acceptance_SHA256=sha(audit), both_original_fine_reference_crosscounts_replayed=True,
                  reference_overlap_never_resolved_by_arbitrary_precedence=True, unique_reference_counts_replayed_with_Counter=True,
                  all_denominators_partitioned=True, files={p.name: dict(SHA256=sha(p), bytes=p.stat().st_size) for p in out.iterdir()},
                  no_cell_rows_or_arrays_saved=True, no_new_clinical_inference_or_authority_label_adoption=True)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], terminal_SHA256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
