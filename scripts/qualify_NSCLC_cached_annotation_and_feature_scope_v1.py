"""Metadata only: every cached contract entry and every existing native axis.

Learned molecular endpoints are read from a private frozen request, never from
stdout or this public source. Quarantined/reference roles stay quarantined or
reference; metadata availability is not endpoint or clinical acceptance.
"""
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
from NSCLC_external_endpoint_text_rules_v1 import gamma_text
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def annotation_field(value):
    s = str(value).casefold()
    if re.search(r'patient|donor|sample|barcode|library|identifier|(^|_)id($|_)', s):
        return False
    return bool(re.search(r'cell.?type|cell.?subtype|sub.?cell.?type|annot|^clust|sub.?clust|leiden|louvain|^state', s))


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 2 * 1024**2)
    request = root / 'NSCLC_external_endpoint_metadata_qualification_20261019_v1/private_followup_source_feature_and_annotation_scope_v1.json'
    assert sha(request) == '895f85e0fc9e86527cfa178fa535994db1146dbf1c8a7875a0c4707330de6d16'
    scope = json.loads(request.read_text())
    audit = root / 'NSCLC_external_endpoint_metadata_audit_20261019_v1/private_terminal.json'
    assert sha(audit) == scope['source_metadata_audit_sha256']
    out = root / 'NSCLC_cached_annotation_feature_qualification_20261019_v1'
    assert not out.exists(); out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    contract = root.parents[2] / 'Project_v5/runs/20261013_immunotherapy_atlas_v1/config/source_contract_v2_private.json'
    config = json.loads(contract.read_text()); records = config['records']
    assert len(records) == 21 and len({v['dataset'] for v in records}) == 21
    save('plan_before_cached_annotation_feature_reads.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), private_scope_sha256=sha(request), request_filesystem_mtime_UTC=datetime.fromtimestamp(request.stat().st_mtime, timezone.utc).isoformat(), request_timestamp_is_manual_not_authoritative=True, resource_guard=guard, source_contract_sha256=sha(contract), scope='All21 existing contract entries including quarantine/reference/supplemental, and all61 old atlas-arm native axes. Observe inherited row selectors, no cell-level output or expression matrix. No selected source or field may be promoted to clinical/identity acceptance.', field_rule='Lexical annotation fields excluding patient/donor/sample/library/barcode/ID. Up to2000 categories per field; more is explicitly unqualified. Gamma/marker text is opportunity only; no new cell annotation.'))
    inputs = read_inputs(root); axes = inputs['source_axes']; assert len(axes) == 61
    features = []; feature_rows = []
    for role, genes in scope['features_by_role'].items():
        features.extend((role, gene) for gene in genes)
    for (arm, dataset), axis in sorted(axes.items()):
        for role, gene in features:
            positions = axis['lookup'].get(gene, set())
            feature_rows.append(dict(arm=arm, dataset=dataset, feature_role=role, feature=gene, physical_matches=len(positions), uniquely_measured=len(positions) == 1, native_columns=sorted(int(v) for v in positions), actual_expression_not_read=True))
    pd.DataFrame(feature_rows).to_csv(out / 'private_complete61_feature_role_qualification.tsv', sep='\t', index=False)
    field_rows = []; category_rows = []; sources = []
    for record in records:
        path = Path(record['path']); rows_path = Path(record['row_positions_ref'])
        rows = np.load(rows_path, allow_pickle=False)
        assert rows.ndim == 1 and np.issubdtype(rows.dtype, np.integer) and len(np.unique(rows)) == len(rows)
        h = h5py.File(path, 'r') if path.suffix == '.h5ad' else zarr.open_group(str(path), mode='r')
        try:
            obs = h['obs']; n = len(read_elem(obs[str(obs.attrs['_index'])]))
            assert not len(rows) or (rows.min() >= 0 and rows.max() < n)
            fields = sorted(k for k in obs.keys() if annotation_field(k))
            source = dict(dataset=record['dataset'], source_role=record['role'], raw_metadata_rows=n, selected_rows=len(rows), row_selector_sha256=sha(rows_path), candidate_annotation_fields=len(fields), annotation_field_value_hashes={})
            for field in fields:
                values = read_elem(obs[field]); assert len(values) == n
                values = pd.Categorical(values)
                if len(values.categories) > 2000:
                    field_rows.append(dict(dataset=record['dataset'], source_role=record['role'], field=field, categories=len(values.categories), status='TOO_MANY_CATEGORIES_NOT_INTERPRETED', explicit_gamma_selected=0, marker_gamma_selected=0, mixed_gamma_selected=0))
                    continue
                labels = [str(v) for v in values.categories]
                codes = values.codes[rows]
                counts = np.bincount(codes[codes >= 0], minlength=len(labels))
                classes = np.array([gamma_text(v) for v in labels], np.int8)
                import hashlib
                digest = hashlib.sha256(json.dumps(labels, ensure_ascii=False).encode() + np.asarray(codes, '<i4').tobytes()).hexdigest()
                source['annotation_field_value_hashes'][field] = digest
                for i, label in enumerate(labels):
                    category_rows.append(dict(dataset=record['dataset'], source_role=record['role'], field=field, author_label=label, text_class=int(classes[i]), selected_cells=int(counts[i])))
                field_rows.append(dict(dataset=record['dataset'], source_role=record['role'], field=field, categories=len(labels), status='METADATA_ONLY_NOT_IDENTITY_ACCEPTED', explicit_gamma_selected=int(counts[classes == 1].sum()), marker_gamma_selected=int(counts[classes == 3].sum()), mixed_gamma_selected=int(counts[classes == 2].sum())))
            sources.append(source)
        finally:
            if isinstance(h, h5py.File): h.close()
        save('private_progress.json', dict(stage='CACHED_SOURCE_ANNOTATION_FIELDS_ONLY', completed_source_entries=len(sources), observed_UTC=datetime.now(timezone.utc).isoformat()))
    pd.DataFrame(field_rows).to_csv(out / 'private_complete21_annotation_field_inventory.tsv', sep='\t', index=False)
    pd.DataFrame(category_rows).to_csv(out / 'private_all_annotation_categories_counts.tsv.gz', sep='\t', index=False, compression='gzip')
    save('private_all_source_metadata_read_proofs.json', sources)
    save('private_native_axis_input_hashes.json', inputs['input_sha256'])
    side = root / 'GSE176021_exact_original_author_annotation_sidecar_20261019_v2/private_author_label_lookup_and_coverage.json'
    proof = root / 'GSE176021_exact_original_author_annotation_audit_20261019_v2/private_terminal.json'
    assert sha(proof) == 'd8d21a1c27465d9b8248857c0028e5106509e43426cbdfef9f44ccb1d71c4169'
    restored = json.loads(side.read_text())
    save('private_restored_author_annotation_category_summary.json', dict(source='GSE176021', author_audit_sha256=sha(proof), label_lookup_sha256=sha(side), branches=[dict(branch=k, categories=len(v['labels']), explicit_gamma_categories=sum(gamma_text(s) == 1 for s in v['labels']), marker_gamma_categories=sum(gamma_text(s) == 3 for s in v['labels'])) for k, v in restored.items()], all_source_fields_missing_is_not_cell_absence=True))
    files = {p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}
    assert sum(v['bytes'] for v in files.values()) < 1024**2
    receipt = dict(status='COMPLETE_ALL61_FEATURE_AXES_AND21_CACHED_ANNOTATION_ENTRY_METADATA_AUDIT_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), scope_sha256=sha(request), source_contract_sha256=sha(contract), files=files, gene_axes=61, cached_source_entries=21, expression_read=False, clinical_adoption=False, identity_adoption=False, scientific_outputs_public=False)
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
