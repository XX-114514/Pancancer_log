"""Full public processed axes/metadata/geometry/count qualification, private I/O."""
import argparse
import csv
import gzip
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def main(a):
    proof = json.loads((a.source / 'download_technical_public.json').read_text())
    assert proof['status'] in ['ALL8_PUBLIC_PROCESSED_FILES_OFFICIAL_MD5_LENGTH_ACCEPTED', 'ALL8_SERVER_PUBLIC_FILES_LENGTH_RANGE_MD5_AND_ACCEPTED_LOCAL_SHA_VERIFIED']
    for f in proof['files']:
        p = a.source / f['relative_path']
        assert p.stat().st_size == f['bytes'] and sha(p) == f['sha256']
    assert not a.output.exists()
    a.output.mkdir(mode=0o700)
    root = a.source / 'single_cell_portal'
    metadata = pd.read_csv(root / 'metadata.txt.gz', sep='\t', skiprows=[1], dtype=str, keep_default_na=False)
    assert metadata.NAME.is_unique
    names = pd.read_csv(root / 'counts/barcodes.tsv.gz', header=None, dtype=str)[0]
    genes = pd.read_csv(root / 'counts/features.tsv.gz', header=None, dtype=str)[0]
    assert names.is_unique and genes.is_unique
    assert len(names) == len(metadata) and set(names) == set(metadata.NAME)
    metadata = metadata.set_index('NAME').loc[names].reset_index()
    metadata['matrix_column'] = np.arange(len(names))
    workbook = openpyxl.load_workbook(a.supplement, read_only=True, data_only=True)
    panelrows = list(workbook['Q_MERFISH_Genes'].values)
    official = [str(x[0]).strip() for x in panelrows[1:] if x[0] is not None]
    companion = [str(x[2]).strip() for x in panelrows[1:] if x[2] is not None]
    assert len(official) == len(set(official)) == 479 and set(official) == set(genes)
    coordinates = []
    scopes = []
    for patient in ['Pt38', 'Pt43', 'Pt67', 'Pt73']:
        p = root / ('spatial_' + patient + '.txt.gz')
        d = pd.read_csv(p, sep='\t', skiprows=[1], dtype={'NAME': str})
        assert d.NAME.is_unique and np.isfinite(d[['X', 'Y']]).all().all()
        assert set(d.NAME) <= set(metadata.NAME)
        match = metadata.set_index('NAME').loc[d.NAME]
        assert match.donor_id.eq(patient).all()
        assert set(d.NAME) == set(metadata.loc[metadata.donor_id.eq(patient), 'NAME'])
        d['donor_id_coordinate_source'] = patient
        coordinates.append(d)
        scopes.append(dict(sample=patient, author_cells=len(d),
                           author_global_ID_coordinate_metadata_bijection=True,
                           coordinate_units='AUTHOR_DISTANCE_UNIT_PHYSICAL_LINK_PENDING',
                           boundaries_or_nuclear_centroids_provided=False))
    coords = pd.concat(coordinates, ignore_index=True)
    assert coords.NAME.is_unique and set(coords.NAME) == set(metadata.NAME)
    metadata = metadata.merge(coords, on='NAME', how='left', validate='one_to_one', sort=False)
    assert np.array_equal(metadata.matrix_column, np.arange(len(names)))
    assert metadata.donor_id.eq(metadata.donor_id_coordinate_source).all()
    totals = np.zeros(len(names), np.int64)
    detected = np.zeros(len(names), np.int64)
    parsed = 0
    with gzip.open(root / 'counts/counts.mtx.gz', 'rt') as stream:
        header = stream.readline().strip()
        assert header.startswith('%%MatrixMarket matrix coordinate') and header.endswith('general')
        while True:
            line = stream.readline()
            if not line.startswith('%'):
                break
        dimensions = tuple(map(int, line.split()))
        assert dimensions[:2] == (len(genes), len(names))
        for chunk in pd.read_csv(stream, sep=r'\s+', header=None, names=['gene', 'cell', 'count'],
                                 dtype={'gene': np.int32, 'cell': np.int32, 'count': np.float64}, chunksize=1000000):
            x = chunk['count'].to_numpy()
            assert np.isfinite(x).all() and (x > 0).all() and np.equal(x, np.rint(x)).all()
            gi, ci = chunk.gene.to_numpy() - 1, chunk.cell.to_numpy() - 1
            assert gi.min() >= 0 and gi.max() < len(genes) and ci.min() >= 0 and ci.max() < len(names)
            totals += np.bincount(ci, weights=x, minlength=len(names)).astype(np.int64)
            detected += np.bincount(ci, minlength=len(names))
            parsed += len(chunk)
            if parsed % 10000000 == 0:
                print(json.dumps(dict(phase='all_raw_MatrixMarket_values', entries=parsed)), flush=True)
    assert parsed == dimensions[2]
    ntx = pd.to_numeric(metadata.ntx).to_numpy()
    ngene = pd.to_numeric(metadata.ngene).to_numpy()
    count_match = bool(np.array_equal(totals, ntx))
    detection_match = bool(np.array_equal(detected, ngene))
    metadata['raw_panel_RNA'] = totals
    metadata['raw_panel_detected_genes'] = detected
    metadata['coordinate_kind'] = 'author_published_cell_point_not_provided_boundary'
    metadata['global_segmentation_key_source'] = 'SCP2510_raw_barcodes_metadata_spatial_exact_NAME'
    metadata['cell_key'] = 'SCP2510::' + metadata.NAME
    metadata['legacy_FOV_Baysor_fragment_key_used'] = False
    metadata.to_csv(a.output / 'private_all_author_cells_registry.tsv.gz', sep='\t', index=False, compression='gzip')
    clinical_rows = list(workbook['B_Patient_Metadata'].values)
    clinical = {str(row[0]).strip(): dict(zip(clinical_rows[0], row)) for row in clinical_rows[1:]
                if row[0] is not None and re.fullmatch(r'\d+', str(row[0]).strip())}
    h = list(workbook['H_MER_Tile_Cell_Counts'].values)
    i = list(workbook['I_MER_Tile_Cell_Comp'].values)
    assert {str(r[1]) for r in h[1:] if r[1]} == {'Patient 38', 'Patient 43', 'Patient 67', 'Patient 73'}
    assert {str(r[0]) for r in i[1:] if r[0]} == {'Pt38', 'Pt43', 'Pt67', 'Pt73'}
    patient_rows = []
    for sample in ['Pt38', 'Pt43', 'Pt67', 'Pt73']:
        row = clinical.get(sample[2:], {})
        patient_rows.append(dict(study='SCP2510', donor_id=sample,
                                 patient_number_sheet_B=sample[2:],
                                 exact_published_H_I_donor_alias_and_B_number=True,
                                 clinical_record_present=bool(row),
                                 RECIST=row.get('RECIST'),
                                 response_evaluable=row.get('RECIST') in ['CR', 'PR', 'SD', 'PD'],
                                 tissue_site=row.get('Tissue Site'), histology=row.get('Histology'),
                                 PFS_days=row.get('PFS (days)'),
                                 PFS_event_source_0_censored=row.get('PFS censored (0=yes)'),
                                 OS_days=row.get('OS (days)'),
                                 OS_event_source_0_notdead=row.get('OS Censored (0=not dead)'),
                                 source='official Supplementary Table1 sheets B/H/I and SCP donor_id',
                                 specimen_phase='PRE_ICI_STUDY_DESIGN_SPECIMEN_SPECIFIC_REVIEW_PENDING'))
    pd.DataFrame(patient_rows).to_csv(a.output / 'private_four_case_source_linkage.tsv', sep='\t', index=False)
    pd.DataFrame(scopes).to_csv(a.output / 'private_all4_segmentation_geometry_qualification.tsv', sep='\t', index=False)
    (a.output / 'private_official_panel_and_companion_channels.json').write_text(json.dumps(dict(
        MERFISH_genes=official, companion_smFISH=companion, features_in_count_matrix=list(genes)), indent=2) + '\n')
    quarantine = [dict(object='legacy_FOV_plus_Baysor_fragments', status='NOT_AUTHOR_GLOBAL_CELL_KEYS_NOT_USED'),
                  dict(object='cell_boundaries_nuclear_centroids', status='NOT_IN_DOWNLOADED_PUBLIC_PROCESSED_FILES'),
                  dict(object='physical_units_and_specimen_phase', status='REQUIRES_AUTHOR_METHOD_EXPORT_LINKAGE')]
    pd.DataFrame(quarantine).to_csv(a.output / 'private_isolation_and_unresolved_registry.tsv', sep='\t', index=False)
    terminal = dict(status='ALL_AUTHOR_PROCESSED_INPUTS_FULL_NUMERIC_AND_EXACT_KEY_REVIEW_COMPLETE',
                    finished_utc=datetime.now(timezone.utc).isoformat(),
                    original_cells=len(names), genes=len(genes), all_numeric_entries=parsed,
                    all4_metadata_coordinates_counts_keys_bijective=True,
                    all_count_values_positive_integers=True,
                    metadata_ntx_matches_raw_counts=count_match,
                    metadata_ngene_matches_stored_entries=detection_match,
                    note='nnz/detected equality remains conditional on unique COO pairs; independent audit checks duplicates',
                    all4_published_donor_case_keys_found=all(x['clinical_record_present'] for x in patient_rows),
                    response_evaluable_cases=sum(x['response_evaluable'] for x in patient_rows),
                    no_spatial_effects_or_clinical_model=True,
                    author_processed_ID_authority_not_raw_transcript_per_cell_reconstruction=True,
                    no_boundary_contact_claim=True,
                    producer_sha256=sha(Path(__file__)),
                    input_download_receipt_sha256=sha(a.source / 'download_technical_public.json'),
                    supplement_sha256=sha(a.supplement),
                    files={p.name: sha(p) for p in a.output.iterdir() if p.is_file()})
    (a.output / 'qualification_technical_public.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps({k:v for k,v in terminal.items() if k not in ['files']}), flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for field in ['source', 'supplement', 'output']:
        p.add_argument('--' + field, type=Path, required=True)
    main(p.parse_args())
