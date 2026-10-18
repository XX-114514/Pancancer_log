"""Qualify existing primary-paper XML/methods and failed clinical-table transfers."""
import hashlib
import json
import os
import xml.etree.ElementTree as ET
from datetime import datetime,timezone
from pathlib import Path
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'GSE176021_published_metadata_access_20261019_v1';out=root/'GSE176021_published_protocol_qualification_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    path=source/'private_original_response.bin';response=json.loads((source/'private_metadata_access_receipt.json').read_text());assert sha(path)==response['sha256'] and response['bytes']==path.stat().st_size==201479
    article=ET.fromstring(path.read_bytes());dois=[x.text for x in article.iter('article-id') if x.attrib.get('pub-id-type')=='doi'];assert dois==['10.1038/s41586-021-03752-4']
    sections={}
    for sec in article.iter('sec'):
        title=sec.find('title');title=' '.join(title.itertext()) if title is not None else ''
        if title in ['Patients and biospecimens','scRNA-seq–TCR-seq']:
            sections[title]=' '.join(' '.join(sec.itertext()).split())
    assert len(sections)==2 and 'two doses of anti-PD-1 (nivolumab) before surgical resection' in sections['Patients and biospecimens'] and 'no more than 10% residual viable tumour cells' in sections['Patients and biospecimens']
    assert 'Cryobanked T cells were thawed' in sections['scRNA-seq–TCR-seq'] and 'highly viable CD3' in sections['scRNA-seq–TCR-seq']
    supp=[]
    for x in article.iter('supplementary-material'):
        href=[v for y in x.iter() for k,v in y.attrib.items() if k.endswith('href')]
        if any(s.endswith('MOESM8_ESM.xlsx') for s in href):supp.append({'id':x.attrib['id'],'href':href,'caption':' '.join(' '.join(x.itertext()).split())})
    assert len(supp)==1
    failures=[]
    for name in ['private_supplementary_tables_access_20261019_v1.json','private_EPMC_supplement_archive_access_20261019_v1.json','private_supplementary_tables_clash_route_access_20261019_v1.json']:
        p=source/name;d=json.loads(p.read_text());assert d['status']=='RESPONSE_NOT_ACCEPTED_NO_LABEL_ADOPTED';failures.append({'receipt':name,'sha256':sha(p),'URL':d['URL'],'error_type':d['error_type'],'route':d.get('route','direct_server_request')})
    assert not (source/'private_original_Supplementary_Tables_1_10.xlsx').exists()
    d={'status':'PRIMARY_ARTICLE_PROTOCOL_AND_3_FAILED_TABLE_ACCESS_ATTEMPTS_QUALIFIED_NO_CASE_MPR_MAP','qualified_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'paper_XML_sha256':sha(path),'paper_original_access_receipt_sha256':sha(source/'private_metadata_access_receipt.json'),'DOI':dois[0],'PMID':'34290408','PMCID':'PMC8338555','methods':sections,'clinical_definition':'Primary tumour MPR: no more than10% residual viable tumour at surgery after neoadjuvant nivolumab. This definition does not provide individual source/case labels.','processing':'Cryobanked highly viable CD3-sorted T cells. Distinguish tumour/normal/LN tissues from longitudinal blood and separate antigen-stimulation experiments; no whole-TME claim or blanket claim that all source cells are cultured.','supplementary_table_link_from_primary_XML':supp,'failed_access_receipts':failures,'case_MPR_map_accepted':False,'clinical_association_computed':False,'gene_or_author_label_reassignment':False,'next':'Clinical branch requires complete original table and exact case/material bridge. Meanwhile use accepted author-reference sidecar for frozen-view identity/composition diagnostics, without retesting closed finite families.'}
    p=out/'private_terminal.json';p.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'status':d['status'],'terminal_sha256':sha(p)}))


if __name__=='__main__':main()
