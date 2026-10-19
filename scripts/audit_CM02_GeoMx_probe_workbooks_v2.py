"""Independent DCC aggregation and complete workbook numeric/axis replay."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import tarfile
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import openpyxl
import pandas as pd

def main(a):
    with gzip.open(a.source/'GSE292098_Hs_R_NGS_WTA_v1.0.pkc.gz','rt') as f:p=json.load(f)
    mapping={}
    for target in p['Targets']:
        for probe in target['Probes']:
            if target['CodeClass'].startswith('Endogenous'):
                assert probe['RTS_ID'] not in mapping or mapping[probe['RTS_ID']]==target['DisplayName']
                mapping[probe['RTS_ID']]=target['DisplayName']
    result=np.load(a.producer/'private_candidate_mapped_DCC_gene_counts.npz');genes=result['gene'].tolist();lookup={g:i for i,g in enumerate(genes)}
    registry=pd.read_csv(a.producer/'private_all586_DCC_gene_mapping_registry.tsv',sep='\t');arrays=result['counts'];assert len(arrays)==586
    count=0;records=0
    with tarfile.open(a.source/'GSE271689_RAW.tar') as t:
        for member in t:
            if not member.isfile():continue
            assert registry.source_file.iloc[count]==member.name
            b=t.extractfile(member).read();b=gzip.decompress(b) if b[:2]==b'\x1f\x8b' else b
            section=b.decode().split('<Code_Summary>')[1].split('</Code_Summary>')[0].strip()
            vector=np.zeros(len(genes),np.int64)
            for row in csv.reader(io.StringIO(section)):
                if row[0]=='RTS_ID':continue
                assert len(row)==2
                value=int(row[1]);assert value>=0
                if row[0] in mapping:vector[lookup[mapping[row[0]]]]+=value
                records+=1
            assert np.array_equal(vector,arrays[count]);count+=1
    cells=0
    for path in a.source.glob('*_0_filter_Greek.xlsx'):
        workbook=openpyxl.load_workbook(path,read_only=True,data_only=True);sheet=workbook['TargetCountMatrix'];sheet.reset_dimensions();rows=sheet.iter_rows(values_only=True);header=next(rows)
        stored=np.load(a.producer/(path.stem+'_private_matrix.npz'))
        assert list(header[1:])==stored['segment'].tolist()
        for i,row in enumerate(rows):
            assert str(row[0])==stored['target'][i]
            assert np.array_equal(np.asarray(row[1:],float),stored['values'][i]);cells+=len(row)-1
        assert i+1==len(stored['target'])
    out=a.output;assert not out.exists();out.mkdir(mode=0o700)
    proof=dict(status='ALL586_DCC_CANDIDATE_GENE_SUMS_AND_ALL_GREEK_MATRIX_VALUES_AXES_ACCEPTED',finished_utc=datetime.now(timezone.utc).isoformat(),
               replayed_probe_records=records,replayed_workbook_values=cells,
               exact_Yale_PKC_version_and_patient_ROI_scale_still_pending=True,
               producer_receipt_sha256=hashlib.sha256((a.producer/'probe_workbook_technical_public.json').read_bytes()).hexdigest(),
               auditor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),no_new_effects=True)
    (out/'independent_GeoMx_audit_technical_public.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--producer',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
