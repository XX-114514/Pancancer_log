"""Exact existing public GEO metadata join. All derived files stay on server.

No expression download, patient output, inferred sex or inferred MMR status.
This is an input/coverage audit, not a clinical association test.
"""
import gzip
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_soft(path):
    records, record = [], None
    with gzip.open(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('^SAMPLE = '):
                if record is not None:
                    records.append(record)
                record = {'gsm': line.strip().split(' = ', 1)[1]}
            elif record is not None and line.startswith('!Sample_characteristics_ch1 = '):
                key, sep, value = line.strip().split(' = ', 1)[1].partition(':')
                if sep:
                    key = key.strip().lower()
                    assert key not in record, 'Duplicate characteristic key'
                    record[key] = value.strip()
    if record is not None:
        records.append(record)
    return pd.DataFrame(records)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    result = root / 'CRC_fine_pair_DIALOGUE_20261014_v1'
    parent = root / 'CM_interface_cancer_association_v1/source_unit_review_v1'
    out = result / 'private_MMR_metadata_audit_v1'
    assert not out.exists(), 'Preserve existing version'
    paths = {
        'accepted_parent_audit': parent / 'independent_audit.json',
        'accepted_sidecar': parent / 'public_library_patient_tissue_processing_sidecar.tsv',
        'cached_public_GEO': parent / 'GSE178341_family.soft.gz',
        'frozen_partition': result / 'private_inputs/patient_library_partition.tsv',
    }
    proof = json.loads(paths['accepted_parent_audit'].read_text())
    assert proof['status'] == 'PASS_ALL265_PUBLIC_METADATA_SIDE_KEYS_EXACT_GSM_LIBRARY_PATIENT_SPECIMEN_TISSUE_PROCESSING_AND_MMR_VALUES'
    side = pd.read_csv(paths['accepted_sidecar'], sep='\t', keep_default_na=False)
    side = side[side.gse_id.eq('GSE178341') & side.mapping_status.eq('exact_public_library_match')]
    assert side.sampleID.is_unique
    chosen = pd.read_csv(paths['frozen_partition'], sep='\t', keep_default_na=False)
    assert len(chosen) == 38 and chosen.analysis_block.is_unique and chosen.sample_uid.is_unique
    joined = chosen.merge(side, left_on='sample_uid', right_on='sampleID', how='left', validate='one_to_one', indicator=True)
    assert joined['_merge'].eq('both').all(), 'Unresolved chosen library'
    assert joined.analysis_block.eq(joined.reviewed_analysis_block).all(), 'Patient block mismatch'
    author = parse_soft(paths['cached_public_GEO'])
    assert author.gsm.is_unique
    fields = ['mmrstatus', 'mmr_ihc', 'patienttypeid', 'specimen_type', 'processing_type',
              'tissuesitesimple', 'tissuesite_detailed', 'tumorstage', 'nodestatussimple',
              'metastasisstatus', 'histologicgradesimple']
    assert set(fields).issubset(author.columns)
    joined = joined.drop(columns='_merge').merge(author[['gsm'] + fields], on='gsm', how='left',
                                               validate='many_to_one', suffixes=('_accepted', '_author'), indicator=True)
    assert joined['_merge'].eq('both').all()
    assert joined.mmrstatus_accepted.eq(joined.mmrstatus_author).all(), 'Independent MMR mismatch'
    assert joined.reviewed_patient_id.eq(joined.patienttypeid).all(), 'Independent author patient mismatch'
    assert joined.processing_type_raw.eq(joined.processing_type).all(), 'Independent processing mismatch'
    invalid = {'', '-', 'na', 'nan', 'n/a', 'unknown', 'not available', 'none'}
    summary = []
    for partition, group in joined.groupby('partition', sort=True):
        for field in fields:
            col = 'mmrstatus_author' if field == 'mmrstatus' else field
            summary.append({'partition': partition, 'field': field, 'patients': len(group),
                            'nonmissing': int((~group[col].str.lower().str.strip().isin(invalid)).sum()),
                            'distinct_values': int(group[col].nunique()),
                            'counts_private': json.dumps(group[col].value_counts().to_dict(), sort_keys=True)})
    out.mkdir()
    joined.drop(columns='_merge').to_csv(out / 'all38_exact_private_author_metadata.tsv', sep='\t', index=False)
    pd.DataFrame(summary).to_csv(out / 'private_covariate_coverage.tsv', sep='\t', index=False)
    diagnostic = pd.read_csv(result / 'private_program_interpretation_v1/all36_private_pathway_nuisance_diagnostics.tsv', sep='\t')
    assert len(diagnostic) == 36
    # The decision aid stays private, including reference enrichment significance.
    decisions = diagnostic[['cell_type', 'direction', 'reference', 'exploratory_BH']].copy()
    decisions['interpretation_only_not_mechanism'] = True
    decisions.to_csv(out / 'private_program_nuisance_decision_aid.tsv', sep='\t', index=False)
    inventory = {'status': 'EXACT_AUTHOR_MMR_INPUT_AUDIT_COMPLETE_RESULTS_PRIVATE',
                 'finished_UTC': datetime.now(timezone.utc).isoformat(),
                 'inputs_sha256': {key: sha(path) for key, path in paths.items()},
                 'producer_sha256': sha(__file__), 'clinical_tests_executed': False,
                 'sex_covariate': 'not_present_in_cached_GEO_characteristics_do_not_infer',
                 'MMR_is_not_immunotherapy_response': True,
                 'holdout_not_independent_of_original_CM_discovery': True,
                 'publication_scope_pending': True,
                 'outputs_sha256': {p.name: sha(p) for p in out.iterdir() if p.is_file()}}
    (out / 'private_inventory.json').write_text(json.dumps(inventory, indent=2) + '\n')
    print('EXACT_AUTHOR_MMR_INPUT_AUDIT_COMPLETE_RESULTS_REMAIN_SERVER')


if __name__ == '__main__':
    main()
