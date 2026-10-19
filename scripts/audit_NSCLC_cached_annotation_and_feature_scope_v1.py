"""Re-read all cached annotation categories/counts and complete feature roles."""
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem
from external_functional_coverage_inputs_v1 import read_inputs
from audit_NSCLC_external_endpoints_metadata_v1 import classify
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 1024**2)
    src = root / 'NSCLC_cached_annotation_feature_qualification_20261019_v1'
    terminal = json.loads((src / 'private_terminal.json').read_text())
    assert terminal['status'] == 'COMPLETE_ALL61_FEATURE_AXES_AND21_CACHED_ANNOTATION_ENTRY_METADATA_AUDIT_PENDING'
    for name, d in terminal['files'].items(): assert sha(src / name) == d['sha256']
    out = root / 'NSCLC_cached_annotation_feature_audit_20261019_v1'
    assert not out.exists(); out.mkdir(mode=0o700)
    scope_path = root / 'NSCLC_external_endpoint_metadata_qualification_20261019_v1/private_followup_source_feature_and_annotation_scope_v1.json'
    assert sha(scope_path) == terminal['scope_sha256']
    scope = json.loads(scope_path.read_text())
    inputs = read_inputs(root)
    assert inputs['input_sha256'] == json.loads((src / 'private_native_axis_input_hashes.json').read_text())
    table = pd.read_csv(src / 'private_complete61_feature_role_qualification.tsv', sep='\t', keep_default_na=False)
    assert len(table) == 61 * sum(len(v) for v in scope['features_by_role'].values())
    seen = set()
    for row in table.itertuples():
        key = (row.arm, row.dataset, row.feature_role, row.feature)
        assert key not in seen; seen.add(key)
        assert row.feature in scope['features_by_role'][row.feature_role]
        positions = inputs['source_axes'][(row.arm, row.dataset)]['lookup'].get(row.feature, set())
        assert row.physical_matches == len(positions) and row.uniquely_measured == (len(positions) == 1)
        assert json.loads(row.native_columns) == sorted(positions)
        assert row.actual_expression_not_read
    contract = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    assert sha(contract) == terminal['source_contract_sha256']
    records = json.loads(contract.read_text())['records']; assert len(records) == 21
    fields = pd.read_csv(src / 'private_complete21_annotation_field_inventory.tsv', sep='\t', keep_default_na=False).set_index(['dataset', 'field'])
    categories = pd.read_csv(src / 'private_all_annotation_categories_counts.tsv.gz', sep='\t', keep_default_na=False)
    proofs = {v['dataset']: v for v in json.loads((src / 'private_all_source_metadata_read_proofs.json').read_text())}
    total_categories = 0; checked_fields = 0
    for record in records:
        path = Path(record['path']); rows_path = Path(record['row_positions_ref']); rows = np.load(rows_path, allow_pickle=False)
        p = proofs[record['dataset']]
        assert p['row_selector_sha256'] == sha(rows_path) and p['source_role'] == record['role'] and p['selected_rows'] == len(rows)
        h = h5py.File(path, 'r') if path.suffix == '.h5ad' else zarr.open_group(str(path), mode='r')
        try:
            obs = h['obs']; n = len(read_elem(obs[str(obs.attrs['_index'])])); assert n == p['raw_metadata_rows']
            expected_fields = []
            for field in obs.keys():
                s = str(field).casefold()
                disallowed = any(x in s for x in ['patient', 'donor', 'sample', 'barcode', 'library', 'identifier']) or re.search(r'(^|_)id($|_)', s)
                matches = re.search(r'cell.?type|cell.?subtype|sub.?cell.?type|annot|^clust|sub.?clust|leiden|louvain|^state', s)
                if not disallowed and matches: expected_fields.append(field)
            assert len(expected_fields) == p['candidate_annotation_fields']
            recorded_fields = set(fields.xs(record['dataset']).index) if record['dataset'] in fields.index.get_level_values(0) else set()
            assert set(expected_fields) == recorded_fields
            for field in expected_fields:
                row = fields.loc[(record['dataset'], field)]
                v = pd.Categorical(read_elem(obs[field])); assert len(v) == n and len(v.categories) == row.categories
                assert row.source_role == record['role']
                if len(v.categories) > 2000:
                    assert row.status == 'TOO_MANY_CATEGORIES_NOT_INTERPRETED'; continue
                names = [str(x) for x in v.categories]
                selected = pd.Series(v).astype('string').iloc[rows]
                counts = selected.value_counts(dropna=True)
                recorded = categories[categories.dataset.eq(record['dataset']) & categories.field.eq(field)].set_index('author_label')
                assert len(recorded) == len(names) and set(recorded.index) == set(names)
                summary = {1: 0, 2: 0, 3: 0}
                for name in names:
                    c = classify(name); amount = int(counts.get(name, 0)); r = recorded.loc[name]
                    assert r.text_class == c and r.selected_cells == amount and r.source_role == record['role']
                    if c: summary[c] += amount
                assert row.explicit_gamma_selected == summary[1] and row.marker_gamma_selected == summary[3] and row.mixed_gamma_selected == summary[2]
                total_categories += len(names); checked_fields += 1
        finally:
            if isinstance(h, h5py.File): h.close()
    assert total_categories == len(categories)
    restored = json.loads((src / 'private_restored_author_annotation_category_summary.json').read_text())
    side = root / 'GSE176021_exact_original_author_annotation_sidecar_20261019_v2/private_author_label_lookup_and_coverage.json'
    proof = root / 'GSE176021_exact_original_author_annotation_audit_20261019_v2/private_terminal.json'
    assert sha(side) == restored['label_lookup_sha256'] and sha(proof) == restored['author_audit_sha256']
    original = json.loads(side.read_text())
    for b in restored['branches']:
        labels = original[b['branch']]['labels']
        assert b['categories'] == len(labels) and b['explicit_gamma_categories'] == sum(classify(x) == 1 for x in labels) and b['marker_gamma_categories'] == sum(classify(x) == 3 for x in labels)
    receipt = dict(status='ACCEPTED_ALL61_NATIVE_FEATURE_ROLES_ALL21_CACHED_ANNOTATION_FIELDS_AND_RESTORED_AUTHOR_CATEGORY_SCOPE', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), feature_role_positions_checked=len(table), cached_source_entries_checked=len(records), annotation_fields_checked=checked_fields, annotation_categories_checked=total_categories, resource_guard=guard, expression_read=False, clinical_or_identity_adoption=False, scientific_outputs_public=False)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file()) < 1024**2
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
