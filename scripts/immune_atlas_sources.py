#!/usr/bin/env python3
"""Freeze existing sources and exact cell selection for an immune therapy atlas.

All row-level registries and physical paths stay in the private server run.
The public report contains only study-level counts and conservative roles.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

import h5py
import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem


NEW_ROLES = {
    'GSE286827': ('core', 'HNSCC', 'whole_tumor', 'ICI'),
    'GSE256326': ('core', 'RMC', 'whole_tumor', 'ICI_and_naive_comparator'),
    'GSE294109': ('core', 'ccRCC', 'whole_tumor', 'ICB_TKI_mixed_cohort'),
    'GSE292046': ('core', 'GBM', 'T_enriched', 'neoantigen_vaccine'),
    'CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea': ('CAR_product', 'LBCL', 'CAR_T_infusion_product', 'CAR_T'),
    'CXG_34f5307e_7b4d_4a48_b68f_2ba844c6414b': ('mechanism_reference', 'Breast', 'metastatic_tumor', 'not_ICI_verified'),
    'CXG_67b6b9ac_b841_42ce_b2d4_ed2ec12327bd': ('mechanism_reference', 'LBCL', 'nonmalignant_LME_subset', 'not_ICI_verified'),
    'CXG_843c1459_e852_429c_9339_17b5760e7716': ('mechanism_reference', 'HGG', 'brain_tumor', 'not_ICI_verified'),
    'CXG_b252b015_b488_4d5c_b16e_968c13e48a2c': ('mechanism_reference', 'HGSOC', 'untreated_tumor', 'untreated'),
}
MISSING = {'', 'nan', 'None', '<NA>', 'NA', '__SOURCE_VALUE_MISSING__', 'unknown', 'Unknown'}


def text_column(frame, choices):
    result = pd.Series('', index=frame.index, dtype=str)
    for column in choices:
        if column not in frame:
            continue
        value = frame[column].astype(str)
        valid = ~value.isin(MISSING)
        use = result.eq('') & valid
        result.loc[use] = value.loc[use]
    return result


def write_new(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)


def read_frame(path, key):
    if path.suffix == '.h5ad':
        with h5py.File(path) as handle:
            return read_elem(handle[key])
    return read_elem(zarr.open_group(str(path), mode='r')[key])


def build(project, output, code_git_sha):
    output.mkdir(parents=True, exist_ok=False, mode=0o700)
    for folder in ['config', 'manifests', 'logs', 'objects', 'figures', 'tables', 'scripts']:
        (output / folder).mkdir(mode=0o700)
    shutil.copyfile(__file__, output/'scripts/immune_atlas_sources.py')
    v5 = project / 'Project_v5'
    collection = v5 / 'runs/20261010_clinical_scrna_expansion_v1'
    parent = v5 / 'runs/20261012_public_preprocessed_zarr_v2'
    audit = json.loads((parent / 'manifests/independent_final_audit.json').read_text())
    assert audit['all_checks_pass'] and audit['n_cells'] == 3180973
    new_zarr = parent / 'outputs/public_preprocessed_union_v1.zarr'
    new_obs = read_frame(new_zarr, 'obs')
    records = []
    metadata = []
    public = []

    def register(name, path, frame, row_positions, counts_key, var_key, role,
                 cancer, capture, therapy, qc, origin, sample_columns, patient_columns,
                 author_columns, count_denominator_column=None):
        selected = frame.iloc[row_positions].copy()
        library = text_column(selected, sample_columns)
        if library.eq('').any():
            raise ValueError(f'Missing original library key: {name}')
        canonical = pd.DataFrame(index=pd.Index(name + '::' + selected.index.astype(str), name='cell_id'))
        canonical['parent_cell_id'] = selected.index.astype(str).to_numpy()
        canonical['original_cell_id'] = text_column(selected, ['original_cell_id', 'barcode']).where(
            text_column(selected, ['original_cell_id', 'barcode']).ne(''), selected.index.astype(str)).to_numpy()
        canonical['source_row_position'] = np.asarray(row_positions, dtype=np.int64)
        canonical['dataset'] = name
        canonical['library_id'] = library.to_numpy()
        canonical['sampleID'] = (name + '::' + library).to_numpy()
        canonical['patient_id'] = text_column(selected, patient_columns).to_numpy()
        canonical['patient_id_source_status'] = 'source_key_requires_review'
        canonical['assay'] = text_column(selected, ['assay', 'modality', 'suspension_type']).to_numpy()
        canonical['tissue'] = text_column(selected, ['tissue', 'tissue_label', 'sample_type']).to_numpy()
        canonical['timepoint'] = text_column(selected, ['timepoint', 'time_label', 'pre_post', 'Timepoint_3', 'Timepoint_2']).to_numpy()
        canonical['treatment'] = text_column(selected, ['treatment', 'Treatment', 'Treatment.or.Mode.of.Action', 'author_therapy']).to_numpy()
        canonical['source_response'] = text_column(selected, ['response', 'response_literature', 'outcome']).to_numpy()
        canonical['source_author_label'] = text_column(selected, author_columns).to_numpy()
        canonical['cancer_type'] = text_column(selected, ['Cancer_type_update']).where(
            text_column(selected, ['Cancer_type_update']).ne(''), cancer).to_numpy()
        canonical['capture_scope'] = capture
        canonical['analysis_role'] = role
        canonical['source_origin'] = origin
        canonical['prior_QC_status'] = qc
        canonical['clinical_validation_eligible'] = False
        canonical['independent_validation_eligible'] = False
        canonical['clinical_endpoint_definition_status'] = 'not_frozen'
        for column in ['patient_identity_status', 'response_source', 'response_definition', 'treatment_source',
                       'residual_doublet_flag', 'ambient_status', 'ambient_high_loss_review',
                       'analysis_correction_applied', 'analysis_total_counts', 'source_final_AC']:
            if column in selected:
                canonical[column] = selected[column].to_numpy()
        metadata.append(canonical)
        positions_path = output / 'manifests' / (name + '_source_rows.npy')
        np.save(positions_path, np.asarray(row_positions, dtype=np.int64))
        record = dict(dataset=name, path=str(path), counts_key=counts_key, var_key=var_key,
                      row_positions_ref=str(positions_path), n_cells=len(canonical),
                      n_libraries=canonical.sampleID.nunique(), role=role, cancer_type=cancer,
                      capture_scope=capture, therapy_background=therapy, prior_QC_status=qc,
                      count_denominator_column=count_denominator_column,
                      clinical_validation_eligible=False, independent_validation_eligible=False)
        if path == new_zarr:
            record['new_parent_dataset'] = name
            record['normalized_key'] = 'X'
            record['ambient_removed_key'] = 'layers/ambient_removed_counts'
        records.append(record)
        public.append({k:record[k] for k in ['dataset', 'n_cells', 'n_libraries', 'role', 'cancer_type',
                                            'capture_scope', 'therapy_background', 'prior_QC_status',
                                            'clinical_validation_eligible', 'independent_validation_eligible']})

    for name, (role, cancer, capture, therapy) in NEW_ROLES.items():
        positions = np.flatnonzero(new_obs['dataset'].astype(str).eq(name).to_numpy())
        assert len(positions) == audit['source_cell_counts'][name]
        register(name, new_zarr, new_obs, positions, 'layers/counts', 'var', role, cancer, capture, therapy,
                 'full_native_QC_DD_ambient_assessment_audited', 'audited_new_public_parent',
                 ['original_library_id'], ['patient_uid', 'patient_id', 'donor_id', 'Patient'],
                 ['author_cell_type', 'CellType', 'merged_cellType', 'cell_type', 'Moreno_Atlas_Projected_ID'])

    old = project / 'Project_v2/Downstream/0716_Data/immunotherapy_4cohort_analysis/20260803_v1'
    old_config = json.loads((old / 'config/pipeline_config.json').read_text())
    raw_root = project / 'Project_v2/Downstream/0716_Data/raw_geo_supplements'
    for name in ['GSE169246', 'GSE176021', 'GSE316195']:
        raw_path = raw_root / name / (name + '_rawdata.h5ad')
        qc_path = project / 'Project_v2' / old_config['inputs'][name]['h5ad']
        raw_obs = read_frame(raw_path, 'obs')
        qc_obs = read_frame(qc_path, 'obs')
        calls_path = old / 'outputs/scsn/doublet' / (name + '_scrublet.tsv.gz')
        calls = pd.read_csv(calls_path, sep='\t', dtype={'original_obs_name':str})
        high = calls.set_index('original_obs_name')['high_confidence_doublet'].astype(str).str.lower().isin(['true', '1'])
        accepted = set(qc_obs.index.astype(str)) - set(high.index[high])
        positions = np.flatnonzero(raw_obs.index.astype(str).isin(accepted))
        if len(positions) != len(accepted):
            raise ValueError(f'QC-to-native cell-ID mapping incomplete: {name}')
        spec = old_config['inputs'][name]
        with h5py.File(raw_path) as handle:
            counts_key = 'layers/counts' if 'counts' in handle['layers'] else 'X'
        register(name, raw_path, raw_obs, positions, counts_key, 'var', 'core', spec['cancer_type'],
                 spec['assay_scope'], 'ICI', 'existing_sample_QC_and_Scrublet_high_confidence_reused',
                 'previously_used_immunotherapy_cohort', ['gsm_id', 'sample_label'] if name == 'GSE176021' else ['sample_label', 'gsm_id'],
                 ['patient_id'], ['cd3_CellType', 'cd8_CellType'])

    icb_path = v5 / 'runs/20261011_ICB_external_collection_inventory_v1/expression_asset_v1/c1d02236-6709-4feb-b7d6-30e2a203e8d1.h5ad'
    icb_obs = read_frame(icb_path, 'obs')
    icb_obs['native_study_library_key'] = icb_obs['PMID_donor_id'].astype(str) + '::' + icb_obs['pre_post'].astype(str)
    for study in sorted(icb_obs.Study_name.astype(str).unique()):
        mask = icb_obs.Study_name.astype(str).eq(study).to_numpy()
        role = 'overlap_quarantine' if study in ['Bi', 'Tirosh', 'Jerby_Arnon'] else 'core_source_review'
        register('ICB_' + study, icb_path, icb_obs, np.flatnonzero(mask), 'raw/X', 'raw/var', role,
                 '|'.join(sorted(icb_obs.loc[mask,'Cancer_type_update'].astype(str).unique())), 'mixed_tumor_original_author_scope',
                 'author_ICB_context_endpoint_definition_review', 'author_processed_count_audit_PASS_QC_history_requires_review',
                 'multi_study_author_compilation', ['native_study_library_key'], ['PMID_donor_id'],
                 ['author_cell_type_update', 'author_cell_type'], count_denominator_column='nCount_RNA')

    additional = project / 'Project_v1/Data/data/raw/Breast/GSE305078/GSE305078_rawdata.h5ad'
    if additional.is_file():
        frame = read_frame(additional, 'obs')
        register('GSE305078', additional, frame, np.arange(len(frame)), 'X', 'var',
                 'additional_QC_pending', 'TNBC', 'tumor_source_scope_review', 'ICB',
                 'native_initial_QC_doublet_ambient_pending', 'existing_public_native_matrix', ['batch'], [], [])
    combined = pd.concat(metadata)
    assert combined.index.is_unique
    combined.to_csv(output / 'manifests/all_source_cell_registry_private.tsv.gz', sep='\t')
    write_new(output / 'config/source_contract_private.json', dict(schema='immune_atlas_source_contract_v1',
              code_git_sha=code_git_sha,
              project_root=str(project), run_root=str(output), records=records,
              old_curated_metadata_ref=str(old / 'metadata/curated/sample_metadata_enriched_v2.tsv'),
              clinical_sidecar_ref=str(collection / 'metadata/ICI_clinical_verified_sidecar.tsv'),
              new_parent_audit_sha256=hashlib.sha256((parent/'manifests/independent_final_audit.json').read_bytes()).hexdigest(),
              all_rows_ref=str(output/'manifests/all_source_cell_registry_private.tsv.gz')))
    report = dict(schema='immune_atlas_source_preflight_public_v1', observed_utc=datetime.now(timezone.utc).isoformat(),
                  code_git_sha=code_git_sha, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  sources=public, n_sources=len(public), n_registry_cells=len(combined),
                  n_core_prior_QC_cells=sum(x['n_cells'] for x in public if x['role']=='core'),
                  n_core_source_review_cells=sum(x['n_cells'] for x in public if x['role']=='core_source_review'),
                  n_extra_QC_pending_cells=sum(x['n_cells'] for x in public if x['role']=='additional_QC_pending'),
                  all_registry_cell_IDs_unique=True, raw_clinical_rows_published=False,
                  independent_validation_frozen=False, authoritative_release=False,
                  stage='source_inventory_complete_counts_and_biological_execution_pending')
    write_new(output / 'manifests/source_preflight_public.json', report)
    print(json.dumps(report), flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--project-root', required=True, type=Path)
    parser.add_argument('--run-root', required=True, type=Path)
    parser.add_argument('--code-git-sha', required=True)
    args=parser.parse_args()
    build(args.project_root, args.run_root, args.code_git_sha)
