"""All spatial inputs: accepted author labels, technical coverage and raw headers."""
import csv
import gzip
import json
import os
import re
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata._io.specs import read_elem

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource

T_LABELS = ['CD4+ T cells', 'CD8+ T cells', 'T cell proliferating', 'T cells', 'Treg']
MPS_LABELS = ['Alveolar macrophages', 'Alveolar mph proliferating', 'Classical monocytes',
              'DC2', 'Interstitial mph prevascular', 'Macrophages', 'Monocytes',
              'Nonclassical monocytes', 'Plasmacytoid DC']
UNKNOWN = ['AUTHOR_NOT_RETRIEVED', 'NOT_IN_AUTHOR_QC_OBJECT', 'Unknown']
MISSING = {'', 'NA', 'NAN', 'NONE', 'NULL', '<NA>'}


def header_roles(stream, name):
    line = stream.readline(16385)
    if len(line) >= 16385:
        return dict(status='HEADER_EXCEEDS_16KiB_NOT_TRUNCATED_INTO_DATA', fields=None, role_fields=[])
    text = line.decode('utf-8-sig')
    header = next(csv.reader([text], delimiter='\t' if name.endswith(('.tsv', '.tsv.gz')) else ','))
    roles = [v for v in header if re.search(r'annot|cell.?type|lineage|subtype|cluster|barcode|entity|fov|cell.?id|patient|donor|subject|core', v, re.I)]
    celltype = [v for v in header if re.search(r'annot|cell.?type|lineage|subtype', v, re.I)]
    return dict(status='HEADER_ONLY_NO_ROWS_READ', fields=len(header), role_fields=roles, possible_identity_fields=celltype, identity_values_or_exact_join_reviewed=False)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 1024**2)
    parent = root / 'complementary_spatial_interface_scope_20261019_v2'
    assert sha(parent / 'private_terminal.json') == '4ff95eb65ce082b27476574512ffe47b2158e59e6a995efb9ef52d9ad5d3e266'
    for name, info in json.loads((parent / 'private_terminal.json').read_text())['files'].items():
        assert sha(parent / name) == info['sha256']
    priority = root / 'gamma_receiver_interface_sensitivity_review_20261019_v1/private_actual_visual_and_priority_acceptance.json'
    assert sha(priority) == 'f90587b180d374d4f92e5cf43d6cf92d1556adb127f55a793669144c2a468ed1'
    manifest = root / 'immune_same_name_family_source_review_v1/spatial_native_input_repair_v2/analysis_input_manifest.tsv'
    assert sha(manifest) == '6fda649c8b2bf4abe010ba8d53d50e495bae29ce8d30648675d9be9bf91dcf83'
    files = pd.read_csv(manifest, sep='\t', keep_default_na=False)
    assay = pd.read_csv(parent / 'private_complete14_interface_molecular_and_geometry_scope.tsv', sep='\t')
    assert len(files) == 14 and len(assay) == 14 and set(files.sample_id) == set(assay.sample_id)
    assay = assay.set_index('sample_id')
    old = root / 'immune_same_name_family_source_review_v1/spatial_contact_native_preflight_v3'
    accepted = json.loads((old / 'independent_audit.json').read_text())
    assert sha(old / 'inventory.json') == accepted['producer_inventory_sha256']
    inventory = json.loads((old / 'inventory.json').read_text())
    old_inputs = {v['sample_id']: v for v in inventory['inputs']}
    state = Path(json.loads((root / 'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    reference_path = state / 'state_spatial_proxy_native_cells_20261016_v1/inventory.json'
    reference = json.loads(reference_path.read_text())
    assert reference['status'] == 'ALL_FROZEN_NONDISCOVERY_CALIBRATION_CELLS_NATIVE_COUNTS_AND_ACCEPTED_REFERENCE_ALIGNED'
    out = root / 'all14_author_spatial_lineage_qualification_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_all14_author_lineage_counts.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), accepted_scope_sha256=sha(parent / 'private_terminal.json'), scientific_priority_sha256=sha(priority), source_inputs=14, T_author_labels=T_LABELS, macrophage_monocyte_DC_author_labels=MPS_LABELS, missing_author_label_not_negative=True, minimum_lineage_cells_per_FOV=20, metadata_header_cap_bytes=16384, no_new_gene_values_or_spatial_effects=True, no_classifier_fitted=True, reference_inventory_sha256=sha(reference_path), reference_label_note='Original18 major labels include overlapping unspecified lymphoid/myeloid categories; harmonization and T-versus-NK challenge remain pending, not a ready gold-standard classifier.', resource_guard=guard))
    records, fov_rows, cache_proofs = [], [], []
    for source in files.itertuples():
        path = Path(source.analysis_path)
        expected = assay.loc[source.sample_id]
        # Parent scope binds each full current H5AD SHA to the frozen panel axis.
        authority = next(v for v in json.loads((parent / 'private_all14_coordinate_input_authority.json').read_text()) if v['sample_id'] == source.sample_id)
        assert sha(path) == authority['current_input_sha256']
        row = dict(dataset_id=source.dataset_id, sample_id=source.sample_id, cancer=source.cancer, native_cells=int(expected.native_cells), accepted_author_cache=False, author_labels_available=False, known_author_QC_cells=np.nan, author_T_QC_cells=np.nan, author_Mac_monocyte_DC_QC_cells=np.nan, known_FOVs_with_T_and_MPS_atleast20=np.nan, gene_proxy_paired_FOVs=0, FOV_lineage_pair_status='AUTHOR_LABELS_UNAVAILABLE_NOT_CELL_ABSENCE', original_resource_gene_proxy_measurable=bool(expected.original_resource_endpoints_measurable), specified_assembly_measurable=bool(expected.strict_ligand_assembly_measurable), selected_geometry_has_source_proof=bool(expected.selected_geometry_has_source_proof), new_target_gene_raw_parity_accepted_by_this_stage=False, fine_gamma_identity_validated=False, absolute_physical_units_validated=False, independent_patient_mapping_reviewed=False, source_input_sha256=authority['current_input_sha256'])
        if source.sample_id in old_inputs:
            rel = 'private_native_cell_caches/' + source.sample_id + '.npz'
            cache = old / rel
            assert sha(cache) == inventory['outputs'][rel]
            assert authority['current_input_sha256'] == old_inputs[source.sample_id]['sha256']
            with np.load(cache, allow_pickle=False) as z:
                labels, qc, fov, xy = z['author_labels'].astype(str), z['qc'].astype(bool), z['fov'].astype(str), z['coords']
                assert len(labels) == len(qc) == len(fov) == row['native_cells'] and xy.shape == (len(labels), 2) and np.isfinite(xy).all()
                known = ~np.isin(labels, UNKNOWN)
                T, M = np.isin(labels, T_LABELS) & qc, np.isin(labels, MPS_LABELS) & qc
                assert not (T & M).any()
                valid_fov = ~np.isin(np.char.upper(np.char.strip(fov)), sorted(MISSING))
                for field in sorted(set(fov[valid_fov])):
                    at = fov == field
                    nt, nm = int((at & T).sum()), int((at & M).sum())
                    fov_rows.append(dict(dataset_id=source.dataset_id, sample_id=source.sample_id, cancer=source.cancer, FOV=field, native_cells=int(at.sum()), QC_cells=int((at & qc).sum()), author_T_QC_cells=nt if known.any() else np.nan, author_MPS_QC_cells=nm if known.any() else np.nan, lineages_both_atleast20=nt >= 20 and nm >= 20 if known.any() else None, original_resource_gene_proxy_measurable=row['original_resource_gene_proxy_measurable'], fine_gamma_identity_validated=False, independent_patient_mapping_reviewed=False))
                count = sum(v['lineages_both_atleast20'] is True for v in fov_rows if v['sample_id'] == source.sample_id)
                row.update(accepted_author_cache=True, author_labels_available=bool(known.any()), known_author_QC_cells=int((known & qc).sum()), author_T_QC_cells=int(T.sum()) if known.any() else np.nan, author_Mac_monocyte_DC_QC_cells=int(M.sum()) if known.any() else np.nan, known_FOVs_with_T_and_MPS_atleast20=count if known.any() else np.nan, gene_proxy_paired_FOVs=count if row['original_resource_gene_proxy_measurable'] else 0, FOV_lineage_pair_status='AUTHOR_LABELS_UNAVAILABLE_NOT_CELL_ABSENCE' if not known.any() else 'AUTHOR_LINEAGE_TECHNICAL_FIELD_OPPORTUNITY' if count else 'NO_PAIRED_FIELD_WITH20_PER_AUTHOR_LINEAGE')
                cache_proofs.append(dict(sample_id=source.sample_id, cache_relative_path=rel, cache_sha256=sha(cache), used_metadata_fields=['author_labels', 'qc', 'fov', 'coords'], expression_arrays_read=False, exact_original_RDS_label_mapping_inherited=True))
        row['broad_resource_proxy_field_opportunity'] = bool(row['author_labels_available'] and row['original_resource_gene_proxy_measurable'] and row['selected_geometry_has_source_proof'] and row['known_FOVs_with_T_and_MPS_atleast20'] > 0)
        records.append(row)
    assert len(records) == 14
    pd.DataFrame(records).to_csv(out / 'private_complete14_author_lineage_and_assay_scope.tsv', sep='\t', index=False)
    pd.DataFrame(fov_rows).to_csv(out / 'private_complete_old5_technical_field_lineage_counts.tsv.gz', sep='\t', index=False, compression='gzip')
    save('private_all5_author_cache_proofs.json', cache_proofs)
    headers = []
    for gse in sorted(files.dataset_id.unique()):
        archive = root.parents[1] / 'data_expansion_20260909/raw/spatial' / gse / (gse + '_RAW.tar')
        if not archive.exists():
            headers.append(dict(dataset_id=gse, archive_present=False, source_label_absence_claimed=False))
            continue
        stat = (archive.stat().st_size, archive.stat().st_mtime_ns)
        with tarfile.open(archive, 'r:*') as tar:
            candidates = [v for v in tar.getmembers() if v.isfile() and re.search(r'metadata|annotation|cell.?type|label|cell.?meta|cluster|(?:^|/)cells\.', v.name, re.I) and re.search(r'\.(csv|tsv|txt)(\.gz)?$', v.name, re.I)]
            assert len(candidates) <= 128, 'Metadata discovery scope exceeds predeclared cap; do not silently truncate'
            for member in candidates:
                raw = tar.extractfile(member)
                stream = gzip.GzipFile(fileobj=raw) if member.name.endswith('.gz') else raw
                headers.append(dict(dataset_id=gse, archive=str(archive), archive_stat_bytes=stat[0], archive_stat_mtime_ns=stat[1], member=member.name, **header_roles(stream, member.name)))
        assert stat == (archive.stat().st_size, archive.stat().st_mtime_ns)
    save('private_complete_bounded_raw_metadata_headers.json', headers)
    # Catalog existing annotation-like sidecars in bounded spatial parents only.
    # An unread RDS/cluster file is a review lead, not accepted cell identity.
    search_roots = {Path(v).parent for v in files.analysis_path}
    search_roots.update(root.parents[1] / 'data_expansion_20260909/raw/spatial' / gse for gse in files.dataset_id.unique())
    candidates = set()
    for folder in sorted(search_roots):
        if not folder.is_dir():
            continue
        children = list(folder.iterdir())
        for sub in [v for v in children if v.is_dir() and re.search(r'annot|metadata|label|cell.?type', v.name, re.I)]:
            children.extend(sub.iterdir())
        for path in children:
            if path.is_file() and re.search(r'annot|metadata|label|cell.?type|cluster', path.name, re.I) and re.search(r'\.(rds|rda|rdata|csv|tsv|txt)(\.gz)?$', path.name, re.I):
                candidates.add(path)
    assert len(candidates) <= 256, 'Bounded sidecar catalog exceeds cap; do not silently truncate'
    save('private_bounded_existing_annotation_sidecar_catalog.json', [dict(path=str(path), bytes=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns, content_or_exact_join_reviewed=False, identity_availability_not_inferred_from_filename=True) for path in sorted(candidates)])
    ref = pd.DataFrame(reference['private_native_files'])
    classes = ref.groupby('major').agg(cells=('cells', 'sum'), sources=('gse_id', 'nunique'), files=('relative_file', 'count')).reset_index()
    assert ref.cells.sum() == reference['cells'] == 91447 and ref.gse_id.nunique() == 52 and len(classes) == 18
    classes.to_csv(out / 'private_complete18_reference_class_scope.tsv', sep='\t', index=False)
    aggregate = pd.DataFrame(records).groupby(['dataset_id', 'cancer'], sort=True).agg(source_inputs=('sample_id', 'size'), author_available_inputs=('author_labels_available', 'sum'), measured_original_resource_inputs=('original_resource_gene_proxy_measurable', 'sum'), broad_resource_proxy_field_opportunity_inputs=('broad_resource_proxy_field_opportunity', 'sum'), T_author_QC_cells=('author_T_QC_cells', lambda x: x.sum(min_count=1)), MPS_author_QC_cells=('author_Mac_monocyte_DC_QC_cells', lambda x: x.sum(min_count=1)), paired_technical_FOVs=('known_FOVs_with_T_and_MPS_atleast20', lambda x: x.sum(min_count=1)), gene_proxy_paired_FOVs=('gene_proxy_paired_FOVs', 'sum')).reset_index()
    aggregate.to_csv(out / 'private_complete7_source_lineage_coverage.tsv', sep='\t', index=False)
    receipt = dict(status='COMPLETE_ALL14_AUTHOR_LINEAGE_AND_SOURCE_HEADER_SCOPE_NO_NEW_EXPRESSION_OR_SPATIAL_TEST', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), resource_guard=guard, spatial_input_positions=14, original_cache_metadata_sources_checked=len(cache_proofs), reference_class_positions=18, scientific_values_public=False, new_expression_values_or_P=False, original_cell_labels_or_coordinates_changed=False, fine_gamma_or_binding_validation=False, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in receipt['files'].values()) < 1024**2
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
