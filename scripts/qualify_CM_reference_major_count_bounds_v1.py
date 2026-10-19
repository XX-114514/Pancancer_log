"""Necessary reference-support bounds from already reviewed source aggregates.

This checks original captured-cell counts only, not selected reference cells,
library support, RNA mass, identity validity or model performance. No individual
input is accepted. All scientific output belongs in the private workspace.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path


def qualify(records, sources, scope):
    grouped = defaultdict(list)
    for row in records:
        if row['scope'] == scope and row['source'] in sources:
            grouped[(row['source'], row['major'])].append(row)
    assert grouped and len(sources) == len(set(sources))
    majors = sorted({major for _, major in grouped})
    assert all((source, major) in grouped for source in sources for major in majors)
    rows, full = [], []
    for (source, major), group in sorted(grouped.items()):
        assert len(group) == 7 and len({r['panel'] for r in group}) == 7
        assert len({r['original_cells'] for r in group}) == 1
        assert len({r['source_scope_cells'] for r in group}) == 1
        n = group[0]['original_cells']
        assert isinstance(n, int) and not isinstance(n, bool) and n >= 0
        full.append(dict(source=source, major=major, cells=n,
                         total=group[0]['source_scope_cells']))
        if major != 'T':
            rows.append(dict(source=source, major=major, original_captured_cells=n,
                             necessary_cell_count_gate=n >= 20,
                             selected_cell_support_unknown=True,
                             library_support_unknown=True))
    totals = []
    for source in sources:
        group = [r for r in full if r['source'] == source]
        assert len({r['total'] for r in group}) == 1
        assert sum(r['cells'] for r in group) == group[0]['total']
        totals.append(dict(source=source, original_captured_cells=group[0]['total']))
    nonT = sorted({r['major'] for r in rows})
    contexts = []
    for held in sources:
        contexts.append(dict(
            source=held,
            nonT_types_below_heldout_20_cell_necessary_gate=[
                r['major'] for r in rows
                if r['source'] == held and not r['necessary_cell_count_gate']],
            nonT_types_with_fewer_than2_training_sources_at20cells=[
                major for major in nonT if sum(
                    r['necessary_cell_count_gate'] for r in rows
                    if r['source'] != held and r['major'] == major) < 2],
            count_feasibility_is_not_full_reference_acceptance=True))
    return dict(
        source_major_positions=rows, source_totals=totals, heldout_contexts=contexts,
        scope=scope, necessary_minimum_original_cells=20,
        required_training_sources=2, seven_panel_count_redundancy_exact=True,
        original_major_partition_conservation=True,
        no_class_removed_or_relabelled=True,
        rare_cell_count_not_RNA_mass_bound=True,
        zero_captured_not_biological_absence=True,
        not_selected_cell_or_library_support_acceptance=True,
        interpretation='Count-only necessary conditions under the stated reference policy; not a test of biological identifiability or clinical value.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-json', type=Path, required=True)
    parser.add_argument('--copy-receipt', type=Path, required=True)
    parser.add_argument('--scope', required=True)
    parser.add_argument('--sources', nargs='+', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    digest = hashlib.sha256(args.source_json.read_bytes()).hexdigest()
    receipt = json.loads(args.copy_receipt.read_text())
    proof = receipt['files'][args.source_json.name]
    assert digest == proof['SHA256'] and args.source_json.stat().st_size == proof['bytes']
    result = qualify(json.loads(args.source_json.read_text()), args.sources, args.scope)
    result.update(status='COMPLETE_SOURCE_AGGREGATE_REFERENCE_COUNT_BOUNDS',
                  completed_UTC=datetime.now(timezone.utc).isoformat(),
                  input_SHA256=digest,
                  source_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  no_new_RNA_or_individual_data_read=True)
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2)
        handle.write('\n')
    print(json.dumps(dict(status=result['status'],
                         output_SHA256=hashlib.sha256(args.output.read_bytes()).hexdigest())))


if __name__ == '__main__':
    main()
