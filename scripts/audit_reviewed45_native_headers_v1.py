"""Independently verify every native header and existing coverage summary.

No producer import, no expression count reads, and no patient/sample exports.
File/header acceptance is explicitly weaker than count or clinical acceptance.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import h5py
import pandas as pd


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8388608), b''):
            h.update(block)
    return h.hexdigest()


def main():
    started = datetime.now(timezone.utc).isoformat()
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    annotation = Path(os.environ['PANCANCER_ANNOTATION_ROOT'])
    parent = root / 'immune_same_name_family_source_review_v1'
    out = root / 'reviewed45_native_header_recon_v1'
    assert not (out / 'independent_audit.json').exists()
    manifest = json.loads((out / 'inventory.json').read_text())
    for name, value in manifest['outputs'].items():
        assert sha(out / name) == value
    membership = pd.read_csv(parent / 'membership_review_v2/source_node_analysis_membership.tsv', sep='\t')
    coverage = pd.read_csv(parent / 'membership_review_v2/all5160_reviewed_family_sample_coverage.tsv', sep='\t')
    pairs = pd.read_csv(parent / 'membership_review_v2/all12_reviewed_context_pair_coverage.tsv', sep='\t')
    scope = pd.read_csv(parent / 'frozen1290_scope.tsv', sep='\t')
    for key, path in {'membership': parent / 'membership_review_v2/source_node_analysis_membership.tsv',
                      'family_coverage': parent / 'membership_review_v2/all5160_reviewed_family_sample_coverage.tsv',
                      'pair_coverage': parent / 'membership_review_v2/all12_reviewed_context_pair_coverage.tsv',
                      'scope': parent / 'frozen1290_scope.tsv'}.items():
        assert sha(path) == manifest['input_sha256'][key]
    member = membership[membership.analysis_family_main_membership]
    ids = sorted(member.gse_id.unique())
    sources = pd.read_csv(out / 'all45_native_header_source_summary.tsv', sep='\t', keep_default_na=False)
    families = pd.read_csv(out / 'all180_reviewed_family_unit_coverage_summary.tsv', sep='\t')
    assert len(ids) == len(sources) == 45 and sources.gse_id.is_unique
    assert set(sources.gse_id) == set(ids) and member.subCluster.nunique() == 13
    assert len(families) == 180 and not families[['gse_id', 'family']].duplicated().any()
    expected_family_keys = {(sid, family) for sid in ids for family in coverage.family.unique()}
    assert set(zip(families.gse_id, families.family)) == expected_family_keys
    assert len(scope) == 1290 and not scope.sampleID.duplicated().any()
    source_checks = []
    for row in sources.itertuples():
        sid = row.gse_id
        native = annotation / '20260923_core_anndata_zarr_v1/objects/harmonized_per_gse' / (sid + '_symbol_counts.h5ad')
        with h5py.File(native, 'r') as handle:
            shape = tuple(map(int, handle['X'].attrs['shape']))
            genes = [str(v) for v in handle['var'][handle['var'].attrs['_index']].asstr()[:]]
            assert shape == (int(row.native_source_cells), int(row.native_genes))
            assert len(genes) == shape[1] == len(set(genes))
            assert hashlib.sha256(json.dumps(genes, ensure_ascii=False, separators=(',', ':')).encode()).hexdigest() == row.gene_axis_json_sha256
            assert str(handle['X/data'].dtype) == row.count_storage_dtype
            assert ';'.join(sorted(handle['obs'].keys())) == row.obs_schema
        assert native.stat().st_size == int(row.native_bytes) and native.stat().st_mtime_ns == int(row.native_mtime_ns)
        m = member[member.gse_id.eq(sid)]
        s = scope[scope.gse_id.eq(sid)]
        f = coverage[coverage.gse_id.eq(sid)]
        p = pairs[pairs.gse_id.eq(sid) & pairs['mode'].eq('marker_supported_main')]
        assert ';'.join(sorted(m.subCluster)) == row.reviewed_node_labels
        assert len(m) == int(row.reviewed_original_nodes) and m.family.nunique() == int(row.reviewed_families)
        assert len(s) == int(row.frozen_scope_libraries) and int(s.n_cells.sum()) == int(row.frozen_scope_cells)
        assert ';'.join(sorted(set(s.cancerType))) == row.frozen_scope_cancer_labels
        assert int(s.unit_kind.eq('reviewed_patient').sum()) == int(row.reviewed_patient_scope_libraries)
        assert json.loads(row.unit_kind_counts_json) == s.unit_kind.value_counts().to_dict()
        assert int(f.marker_supported_n_cells.ge(20).sum()) == int(row.family_library_units_ge20)
        assert int(p.n_samples.ge(20).sum()) == int(row.pair_strata_records_ge20_samples)
        assert int(p.n_samples.ge(10).sum()) == int(row.pair_strata_records_ge10_samples)
        assert row.status == 'HEADER_AVAILABLE_NOT_FULL_COUNT_OR_CLINICAL_ACCEPTANCE'
        for fr in families[families.gse_id.eq(sid)].itertuples():
            subset = f[f.family.eq(fr.family)]
            assert len(subset) == fr.frozen_library_units
            assert int(subset.marker_supported_n_cells.ge(20).sum()) == fr.libraries_marker_supported_ge20
            assert int(subset.marker_supported_n_cells.sum()) == fr.marker_supported_cells
        source_checks.append(dict(gse_id=sid, status='PASS_NATIVE_HEADER_AND_EXISTING_UNIT_COVERAGE'))
    assert int(sources.native_source_cells.sum()) == manifest['native_source_cells_sum']
    assert int(sources.frozen_scope_libraries.sum()) == manifest['scope_libraries']
    assert int(sources.family_library_units_ge20.sum()) == manifest['family_library_units_ge20']
    assert int(sources.pair_strata_records_ge20_samples.gt(0).sum()) == manifest['sources_with_ge20_pair_stratum_records']
    assert int(sources.pair_strata_records_ge10_samples.gt(0).sum()) == manifest['sources_with_ge10_pair_stratum_records']
    result = dict(status='PASS_ALL45_HEADERS_ALL180_EXISTING_FAMILY_UNIT_SUMMARIES',
                  started_UTC=started, finished_UTC=datetime.now(timezone.utc).isoformat(),
                  auditor_sha256=sha(__file__), producer_inventory_sha256=sha(out / 'inventory.json'),
                  source_headers=45, source_family_summaries=180, checks=source_checks,
                  count_values_or_full_file_sha_verified=False, clinical_units_accepted=False,
                  formal_states=False, scientific_release=False, axis_validation_complete=False)
    (out / 'independent_audit.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
