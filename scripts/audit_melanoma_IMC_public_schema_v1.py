"""Re-read public schema and reconcile protein aliases; no individual rows."""
from datetime import datetime,timezone
from pathlib import Path
import csv
import hashlib
import io
import json
import os
import re
import zipfile

import preflight_melanoma_IMC_public_archive_v2 as source


def sha(b):return hashlib.sha256(b).hexdigest()


def norm(s):return re.sub('[^a-z0-9]','',s.lower())


def main():
    prior=Path(os.environ['PANCANCER_PUBLIC_PREFLIGHT_SOURCE'])
    out=Path(os.environ['PANCANCER_PUBLIC_PREFLIGHT_AUDIT']);assert not out.exists();out.mkdir(mode=0o700)
    terminal=json.loads((prior/'private_terminal.json').read_text())
    assert terminal['source_SHA256']==sha(Path(source.__file__).read_bytes())
    for n,p in terminal['files'].items():assert sha((prior/n).read_bytes())==p['SHA256']
    contract=json.loads((prior/'private_official_panel_and_file_contract.json').read_text())
    d=json.loads(source.fetch(source.BASE+f'api/records/{source.RECORD}',2*1024**2))
    files={x['key']:x for x in d['files']}
    assert all(files[n]['checksum']==f['checksum'] and files[n]['size']==f['size'] for n,f in contract['files'].items())
    raw=source.fetch(files['Melanoma_panel.csv']['links']['self'],1024**2)
    assert sha(raw)==contract['panel_SHA256']
    panel=list(csv.DictReader(io.StringIO(raw.decode('utf-8-sig'))))
    aliases={'CD278':['ICOS'],'CD274':['PDL1'],'CD134':['OX40'],'CD223':['LAG3'],
        'CD279':['PD1'],'CD366':['TIM3'],'CD127':['IL7R'],'FoxP3':['FOXP3'],
        'Ki67':['Ki-67','MKI67'],'E_Cadherin':['E-Cadherin'],'Collagen_I':['CollagenI'],
        'CA_IX':['CAIX'],'HLA_DR':['HLA-DR'],'Alpha_SMA':['aSMA','alphaSMA','SMA']}
    f=files['Melanoma_code_data.zip'];rf=source.RangeFile(f['links']['self'],f['size'])
    schemas=[];logical_metadata=[]
    with zipfile.ZipFile(rf) as z:
        entries=z.infolist()
        signature=sha(json.dumps([(e.filename,e.file_size,e.compress_size,e.CRC) for e in entries],separators=(',',':')).encode())
        assert signature==json.loads((prior/'private_complete_directory_schema.json').read_text())['directory_identity_SHA256']
        targets=[e for e in entries if e.filename.endswith('/input/sc_data.csv') or e.filename=='input/sc_data.csv']
        assert len(targets)==1
        for e in targets:
            with z.open(e) as stream:header=stream.readline(65537)
            assert len(header)<=65536
            cols=next(csv.reader([header.decode('utf-8-sig').strip()]))
            assert sha(header)==json.loads((prior/'private_selected_header_schema.json').read_text())[0]['header_SHA256']
            assert len(cols)==len(set(cols)) and 'sample_id' in cols
            mapping={}
            for r in panel:
                possible={r['marker'],r['label'],r['Metal'],r['Metal']+'Di'}|set(aliases.get(r['marker'],[]))
                hits=[c for c in cols if norm(c) in {norm(p) for p in possible}]
                assert len(hits)<=1
                mapping[r['marker']]=hits[0] if hits else None
            schemas.append(dict(logical_role='processed_single_cell_protein_table',header_SHA256=sha(header),
                column_names_schema_only=cols,panel_to_processed_protein_aliases=mapping,
                declared_uncompressed_bytes=e.file_size,declared_compressed_bytes=e.compress_size,
                no_values_or_clinical_records_read=True,full_member_CRC_not_checked=True))
        # Directory-only metadata candidates: never return coded file names.
        for e in entries:
            if e.is_dir() or not e.filename.lower().endswith(('.csv','.tsv','.txt')):continue
            name=e.filename.lower()
            if '/rnaseq' in name or '/figs/' in name:continue
            if '/clidata/' in name:role='bulk_validation_clinical_files_not_imaging_bridge_accepted'
            elif any(k in name for k in ['meta','anno','sample','spatial']):role='possible_imaging_link_or_spatial_schema'
            else:continue
            logical_metadata.append(dict(logical_role=role,extension=Path(e.filename).suffix,
                bytes=e.file_size,compressed_bytes=e.compress_size,path_SHA256=sha(e.filename.encode())))
    proof=dict(status='COMPLETE_INDEPENDENT_PUBLISHED_PANEL_DIRECTORY_HEADER_AND_ALIAS_RECONCILIATION',
        completed_UTC=datetime.now(timezone.utc).isoformat(),audit_source_SHA256=sha(Path(__file__).read_bytes()),
        producer_terminal_SHA256=sha((prior/'private_terminal.json').read_bytes()),
        complete_directory_schema_replayed=True,source_panel_physical_markers_and_exact_processed_aliases=schemas,
        directory_only_candidate_roles=logical_metadata,
        clinical_patient_phase_response_bridge_not_accepted=True,protein_positivity_thresholds_not_accepted=True,
        no_RNA_recipe_equivalence_or_new_cell_identity_claim=True,
        no_cell_clinical_rows_or_raw_images_or_learned_arrays_read_or_saved=True,
        no_author_code_executed=True,no_effect_tests=True,range_bytes=rf.transferred)
    target=out/'private_complete_schema_alias_audit.json';target.write_text(json.dumps(proof,indent=2)+'\n')
    print(json.dumps(dict(status=proof['status'],receipt_SHA256=sha(target.read_bytes()),schemas=schemas)))


if __name__=='__main__':main()
