"""Alternate full COO parser, duplicate bitmap, all cell/clinical key replay."""
import argparse
import gzip
import hashlib
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2), b''):
            h.update(b)
    return h.hexdigest()


def main(a):
    proof = json.loads((a.producer/'qualification_technical_public.json').read_text())
    for name, digest in proof['files'].items():
        assert sha(a.producer/name) == digest
    assert not a.output.exists()
    a.output.mkdir(mode=0o700)
    base = a.source/'single_cell_portal'
    with gzip.open(base/'counts/barcodes.tsv.gz','rt') as f:
        names = [s.rstrip('\n\r') for s in f]
    with gzip.open(base/'counts/features.tsv.gz','rt') as f:
        genes = [s.rstrip('\n\r') for s in f]
    registry = pd.read_csv(a.producer/'private_all_author_cells_registry.tsv.gz',sep='\t',dtype={'NAME':str})
    assert registry.NAME.tolist() == names and registry.cell_key.is_unique
    assert registry.matrix_column.tolist() == list(range(len(names)))
    totals = np.zeros(len(names),np.int64)
    detected = np.zeros(len(names),np.int64)
    bitmap = np.zeros((len(names)*len(genes)+7)//8,np.uint8)
    entries = 0
    with gzip.open(base/'counts/counts.mtx.gz','rt') as f:
        header = f.readline()
        while True:
            line = f.readline()
            if not line.startswith('%'):
                dimensions = list(map(int,line.split()));break
        assert dimensions[:2] == [len(genes),len(names)]
        while True:
            lines = list(itertools.islice(f,200000))
            if not lines:
                break
            block = np.fromstring(''.join(lines),sep=' ',dtype=np.float64).reshape(-1,3)
            assert len(block) == len(lines) and np.isfinite(block).all()
            assert np.equal(block,np.rint(block)).all()
            gi,ci,values = block.astype(np.int64).T
            gi -= 1;ci -= 1
            assert gi.min()>=0 and gi.max()<len(genes) and ci.min()>=0 and ci.max()<len(names) and values.min()>0
            keys = ci*len(genes)+gi
            assert len(np.unique(keys)) == len(keys), 'Within-chunk duplicate COO pairs'
            byte = keys//8
            bit = np.left_shift(np.uint8(1),(keys%8).astype(np.uint8))
            assert not np.any(np.bitwise_and(bitmap[byte],bit)), 'Cross-chunk duplicate COO pairs'
            np.bitwise_or.at(bitmap,byte,bit)
            np.add.at(totals,ci,values)
            np.add.at(detected,ci,1)
            entries += len(block)
    assert entries == dimensions[2] == proof['all_numeric_entries']
    assert np.array_equal(totals,registry.raw_panel_RNA)
    assert np.array_equal(detected,registry.raw_panel_detected_genes)
    assert np.array_equal(totals,pd.to_numeric(registry.ntx))
    assert np.array_equal(detected,pd.to_numeric(registry.ngene))
    original_meta = pd.read_csv(base/'metadata.txt.gz',sep='\t',skiprows=[1],dtype=str,keep_default_na=False).set_index('NAME')
    assert set(original_meta.index) == set(names)
    for col in ['donor_id','library','biosample_id','hubID','hubType','type_lvl1','type_lvl2']:
        assert registry[col].fillna('').astype(str).tolist() == original_meta.loc[names,col].tolist()
    seen = set()
    for case in ['Pt38','Pt43','Pt67','Pt73']:
        d = pd.read_csv(base/('spatial_'+case+'.txt.gz'),sep='\t',skiprows=[1],dtype={'NAME':str})
        assert not seen.intersection(d.NAME);seen.update(d.NAME)
        aligned = registry.set_index('NAME').loc[d.NAME]
        assert aligned.donor_id.eq(case).all()
        np.testing.assert_array_equal(aligned[['X','Y']].to_numpy(),d[['X','Y']].to_numpy())
    assert seen == set(names)
    book = openpyxl.load_workbook(a.supplement,read_only=True,data_only=True)
    clinical = list(book['B_Patient_Metadata'].values)
    source = {str(row[0]):dict(zip(clinical[0],row)) for row in clinical[1:] if row[0] is not None}
    linked = pd.read_csv(a.producer/'private_four_case_source_linkage.tsv',sep='\t',dtype={'patient_number_sheet_B':str})
    for row in linked.itertuples():
        assert row.RECIST == source[row.patient_number_sheet_B]['RECIST']
    result = dict(status='ALL_PUBLIC_PROCESSED_COO_PAIRS_COUNTS_GEOMETRY_METADATA_AND_CASE_LINKS_INDEPENDENTLY_ACCEPTED',
                  finished_utc=datetime.now(timezone.utc).isoformat(),all_cells=len(names),all_genes=len(genes),
                  all_COO_entries=entries,duplicate_COO_pairs=0,all4_original_author_key_geometry_replayed=True,
                  original_metadata_ntx_ngene_exact=True,
                  producer_terminal_sha256=sha(a.producer/'qualification_technical_public.json'),
                  auditor_sha256=sha(Path(__file__)),no_raw_transcript_to_processed_cell_claim=True,
                  no_physical_boundary_contact_or_clinical_effect_claim=True)
    (a.output/'independent_audit_technical_public.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    for field in ['source','producer','supplement','output']:
        p.add_argument('--'+field,type=Path,required=True)
    main(p.parse_args())
