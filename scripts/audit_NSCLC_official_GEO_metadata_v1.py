"""Independent original SOFT replay; exact deposit bridge, no patient aliases."""
import gzip
import hashlib
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,2*1024**2);src=root/'NSCLC_official_GEO_metadata_qualification_20261019_v1';out=root/'NSCLC_official_GEO_metadata_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    t=json.loads((src/'private_terminal.json').read_text());assert t['status']=='COMPLETE_OFFICIAL_GEO_METADATA_SCHEMA_AND_EXACT_TOKEN_DIAGNOSTIC_NO_ALIAS_ADOPTED';assert t['source_sha256']=='6e03375d2674b322508086ae1e60391104c361935e2ce7e7ae13777722966679'
    for name,r in t['files'].items():assert sha(src/name)==r['sha256']
    payload=json.loads((src/'private_all_original_family_metadata.json').read_text());expected=payload['sample_records'];records=[];series={};active=None
    with gzip.open(src/'private_original_GSE243013_family.soft.gz','rt') as stream:
        for line in stream:
            line=line.rstrip('\r\n')
            if line.startswith('^SAMPLE = '):active={'GSM':line[len('^SAMPLE = '):],'fields':{}};records.append(active)
            elif line.startswith('^'):active=None
            elif ' = ' in line:
                k,v=line.split(' = ',1)
                if active is not None and k.startswith('!Sample_'):active['fields'].setdefault(k,[]).append(v)
                elif k.startswith('!Series_'):series.setdefault(k,[]).append(v)
    assert records==expected and series==payload['series_fields'] and len(records)==243
    clinical=pd.read_csv(root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv',sep='\t',dtype=str,keep_default_na=False);native=set(clinical.sampleID);assert len(native)==243
    assert set(r['fields']['!Sample_title'][0] for r in records)==native and len({r['GSM'] for r in records})==243
    mapping=[];characteristics=set();protocols=[]
    for r in records:
        f=r['fields'];assert len(f['!Sample_title'])==1 and f['!Sample_geo_accession']==[r['GSM']]
        title=f['!Sample_title'][0];mapping.append(dict(GSM=r['GSM'],original_author_record=title,bridge='EXACT_ORIGINAL_GEO_TITLE',reviewed_independent_patient=False))
        protocols.append(f['!Sample_extract_protocol_ch1'])
        for value in f['!Sample_characteristics_ch1']:characteristics.add(value.split(':',1)[0].strip())
    assert characteristics=={'tissue','disease state'} and all(p==protocols[0] for p in protocols) and len(protocols[0])==2
    assert '7AAD−CD235a−' in protocols[0][0] and 'FACS' in protocols[0][0] and 'CD45' not in '\n'.join(protocols[0])
    assert any('samples(n=234)' in text and '234' in text for text in series['!Series_overall_design'])
    pd.DataFrame(mapping).to_csv(out/'private_all243_exact_GEO_author_record_bridge.tsv',sep='\t',index=False)
    shared={'public_protocol_sections':protocols[0],'meaning':'Two shared protocol paragraphs for all243records, not two alternative protocols. Documented live/nonerythrocytic7AAD-negative/CD235a-negative FACS; CD45 enrichment not stated. Immune-only analysis scope comes from supplied processed B/T/My object.','official_design_vs_deposit':'Original design says234samples/234patients;243GSMs/exact native author records. No explicit patient/subject/donor characteristic or survival field. Counts alone do not establish duplicate cases, map9controls or upgrade author records to independent patients.','no_running_clinical_change':True}
    (out/'private_shared_protocol_and_identity_limits.json').write_text(json.dumps(shared,indent=2)+'\n')
    receipt={'status':'ACCEPTED_ALL243_ORIGINAL_GEO_RECORDS_EXACT_AUTHOR_TITLES_SHARED_PROTOCOL_NO_PATIENT_ALIAS','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(src/'private_terminal.json'),'resource_guard':guard,'original_GEO_records':243,'exact_native_author_bridges':243,'shared_protocol_paragraphs':2,'explicit_patient_identity_fields':0,'adopted_patient_aliases':0,'running_clinical_input_changed':False,'files':{p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()}}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
