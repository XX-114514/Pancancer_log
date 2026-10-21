"""Independent stdlib XML replay of source clinical/modality coverage, RAM only."""
from collections import Counter
from decimal import Decimal,InvalidOperation
from pathlib import Path
import hashlib
import posixpath
import re
import xml.etree.ElementTree as ET
import zipfile

NS='{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'


def audit(path,expected):
    assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==expected['original_workbook_SHA256']
    with zipfile.ZipFile(path) as z:
        strings=[]
        if 'xl/sharedStrings.xml' in z.namelist():
            strings=[''.join(t.text or '' for t in si.iter(NS+'t')) for si in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        rel={r.attrib['Id']:r.attrib['Target'] for r in ET.fromstring(z.read('xl/_rels/workbook.xml.rels'))}
        sheets={}
        for s in ET.fromstring(z.read('xl/workbook.xml')).find(NS+'sheets'):
            target=rel[s.attrib['{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id']]
            sheets[s.attrib['name']]=target.lstrip('/') if target.startswith('/') else posixpath.normpath('xl/'+target)
        def value(c):
            v=c.find(NS+'v');text=v.text if v is not None else None
            if c.attrib.get('t')=='s':return strings[int(text)]
            if c.attrib.get('t')=='inlineStr':return ''.join(t.text or '' for t in c.iter(NS+'t'))
            return text
        def getrows(sheet,key_only=False):
            ans=[]
            for r in ET.fromstring(z.read(sheets[sheet])).find(NS+'sheetData'):
                cells={re.sub(r'\d','',c.attrib['r']):value(c) for c in r if not key_only or re.sub(r'\d','',c.attrib['r'])=='A'}
                ans.append(cells)
            return ans
        def ident(v):
            if v is None or not re.fullmatch(r'[1-9][0-9]*(?:\.0+)?',v):return None
            return str(int(Decimal(v)))
        b=getrows('B_Patient_Metadata');headers=b[0];registry={}
        for r in b[1:]:
            k=ident(r.get('A'))
            if k is not None:
                assert k not in registry
                registry[k]={headers[a]:v for a,v in r.items()}
        def counts(ids):
            rows=[registry[k] for k in ids]
            responses=Counter(r.get('RECIST') if r.get('RECIST') in ['CR','PR','SD','PD'] else 'MISSING_OR_UNRECOGNIZED' for r in rows)
            outcomes={}
            for name,duration,event in [('PFS','PFS (days)','PFS censored (0=yes)'),('OS','OS (days)','OS Censored (0=not dead)')]:
                valid=[]
                for r in rows:
                    try:t=Decimal(r.get(duration) or 'NaN');e=Decimal(r.get(event) or 'NaN')
                    except InvalidOperation:continue
                    if t.is_finite() and e.is_finite() and t>0 and e in [0,1]:valid.append(e)
                outcomes[name]=(len(valid),sum(e==1 for e in valid))
            return len(rows),dict(responses),outcomes
        def check(ids,ref):
            n,res,outcomes=counts(ids)
            assert n==ref['patients'] and res==ref['RECIST_counts']
            for name,(n,e) in outcomes.items():
                rr=ref['reported_survival_metadata'][name]
                assert n==rr['positive_reported_duration_and_binary_flag'] and e==rr['events_code1']
        check(set(registry),expected['all_author_clinical'])
        groups=[]
        for sheet,label in [('C_ISH_Panel_Freq_of_DAPI','co_registered_RNA'),('D_Ab_Panel_Freq_of_DAPI','co_registered_protein')]:
            raw=getrows(sheet,True)[1:];values=[ident(r.get('A')) for r in raw if r.get('A') is not None]
            ids={v for v in values if v is not None};ref=expected['modalities'][label]
            assert len(values)==ref['published_patient_rows'] and len(ids)==ref['strictly_parsed_unique_patient_keys']
            assert len(ids-set(registry))==ref['missing_author_clinical_key']
            check(ids&set(registry),ref['clinical']);groups.append(ids)
        assert len(groups[0]&groups[1])==expected['RNA_protein_shared_published_patient_keys']
        check(groups[0]&groups[1]&set(registry),expected['both_modalities_clinical'])
    return dict(status='PASS_INDEPENDENT_XML_PATIENT_KEYS_RECIST_AND_SURVIVAL_METADATA',
        separate_parser_and_decimal_numeric_rules=True,no_individual_keys_or_arrays_returned=True,
        no_new_frequency_or_expression_or_effect_values_selected=True)
