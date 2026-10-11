"""Read native headers and existing reviewed coverage; no cell values/models.

Outputs contain source/family summaries, never sample or patient rows.
An existing output directory is rejected. File presence is not acceptance.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os

import h5py
import pandas as pd


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(8388608), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    start = datetime.now(timezone.utc).isoformat()
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    annotation = Path(os.environ['PANCANCER_ANNOTATION_ROOT'])
    parent = root / 'immune_same_name_family_source_review_v1'
    output = root / 'reviewed45_native_header_recon_v1'
    assert not output.exists()
    inputs = {
        'membership': parent / 'membership_review_v2/source_node_analysis_membership.tsv',
        'family_coverage': parent / 'membership_review_v2/all5160_reviewed_family_sample_coverage.tsv',
        'pair_coverage': parent / 'membership_review_v2/all12_reviewed_context_pair_coverage.tsv',
        'scope': parent / 'frozen1290_scope.tsv',
    }
    hashes = {key: sha(path) for key, path in inputs.items()}
    member = pd.read_csv(inputs['membership'], sep='\t')
    accepted = member[member.analysis_family_main_membership.eq(True)].copy()
    source_ids = sorted(set(accepted.gse_id))
    assert len(source_ids) == 45 and accepted.subCluster.nunique() == 13
    assert not set(source_ids) & {'GSE169609', 'GSE299623'}
    family = pd.read_csv(inputs['family_coverage'], sep='\t')
    pair = pd.read_csv(inputs['pair_coverage'], sep='\t')
    scope = pd.read_csv(inputs['scope'], sep='\t')
    assert len(family) == 5160 and len(scope) == 1290
    assert not family[['sampleID', 'family']].duplicated().any()
    assert not scope.sampleID.duplicated().any()
    native = annotation / '20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'
    rows = []
    for source in source_ids:
        path = native / (source + '_symbol_counts.h5ad')
        members = accepted[accepted.gse_id.eq(source)]
        units = scope[scope.gse_id.eq(source)]
        source_families = family[family.gse_id.eq(source)]
        pairs = pair[pair.gse_id.eq(source) & pair['mode'].eq('marker_supported_main')]
        status = 'HEADER_AVAILABLE_NOT_FULL_COUNT_OR_CLINICAL_ACCEPTANCE'
        with h5py.File(path, 'r') as handle:
            genes = handle['var'][handle['var'].attrs['_index']].asstr()[:].tolist()
            shape = tuple(map(int, handle['X'].attrs['shape']))
            gene_hash = hashlib.sha256(json.dumps(genes, ensure_ascii=False,
                                                 separators=(',', ':')).encode()).hexdigest()
            assert len(genes) == shape[1] and len(set(genes)) == len(genes)
            encoding = handle['X'].attrs['encoding-type']
            if isinstance(encoding, bytes):
                encoding = encoding.decode()
            assert encoding == 'csr_matrix'
            obs_fields = ';'.join(sorted(handle['obs'].keys()))
            dtype = str(handle['X/data'].dtype)
        rows.append(dict(
            gse_id=source, native_source_cells=shape[0], native_genes=shape[1],
            gene_axis_unique=True, gene_axis_json_sha256=gene_hash,
            native_bytes=path.stat().st_size, native_mtime_ns=path.stat().st_mtime_ns,
            count_storage_dtype=dtype, matrix_encoding=encoding, obs_schema=obs_fields,
            reviewed_original_nodes=len(members), reviewed_families=members.family.nunique(),
            reviewed_node_labels=';'.join(sorted(members.subCluster)),
            frozen_scope_libraries=len(units), frozen_scope_cancer_labels=';'.join(sorted(set(units.cancerType))),
            frozen_scope_cells=int(units.n_cells.sum()),
            reviewed_patient_scope_libraries=int(units.unit_kind.eq('reviewed_patient').sum()),
            unit_kind_counts_json=json.dumps(units.unit_kind.value_counts().to_dict(), sort_keys=True),
            family_library_units_ge20=int(source_families.marker_supported_n_cells.ge(20).sum()),
            pair_strata_records_ge20_samples=int(pairs.n_samples.ge(20).sum()),
            pair_strata_records_ge10_samples=int(pairs.n_samples.ge(10).sum()),
            pair_records_are_not_independent_evidence=True, status=status))
    output.mkdir()
    table = pd.DataFrame(rows)
    table.to_csv(output / 'all45_native_header_source_summary.tsv', sep='\t', index=False)
    family_rows = family[family.gse_id.isin(source_ids)].groupby(['gse_id', 'family']).agg(
        frozen_library_units=('sampleID', 'size'),
        libraries_marker_supported_ge20=('marker_supported_n_cells', lambda value: int(value.ge(20).sum())),
        marker_supported_cells=('marker_supported_n_cells', 'sum')).reset_index()
    family_rows.to_csv(output / 'all180_reviewed_family_unit_coverage_summary.tsv', sep='\t', index=False)
    assert len(family_rows) == 180
    manifest = dict(
        schema='reviewed45_native_header_recon_v1', status='COMPLETE_HEADERS_AND_EXISTING_COVERAGE_ONLY',
        started_UTC=start, finished_UTC=datetime.now(timezone.utc).isoformat(),
        code_sha256=sha(__file__), input_sha256=hashes, sources=45,
        reviewed_original_node_union=13, native_source_cells_sum=int(table.native_source_cells.sum()),
        native_genes_min=int(table.native_genes.min()), native_genes_max=int(table.native_genes.max()),
        scope_libraries=int(table.frozen_scope_libraries.sum()),
        family_library_units_ge20=int(table.family_library_units_ge20.sum()),
        sources_with_ge20_pair_stratum_records=int(table.pair_strata_records_ge20_samples.gt(0).sum()),
        sources_with_ge10_pair_stratum_records=int(table.pair_strata_records_ge10_samples.gt(0).sum()),
        original_node_union_labels=sorted(set(accepted.subCluster)),
        count_values_or_full_file_sha_verified=False, new_models_fitted=False,
        clinical_units_accepted=False, formal_states=False, scientific_release=False,
        axis_validation_complete=False,
        limits=['Native file/header availability is not full count/identity/clinical acceptance',
                'Existing directed pair coverage records and source/cancer strata are not independent biological replicates',
                'Exact clinical patient block and discovery/validation overlap require separate review',
                'Original numbered subtypes remain separate; no formal state assignment or canonical merge'],
        outputs={path.name: sha(path) for path in output.iterdir() if path.is_file()})
    (output / 'inventory.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest), flush=True)


if __name__ == '__main__':
    main()
