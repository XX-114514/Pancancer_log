"""Bounded official GEO metadata only; diagnose explicit identity evidence privately.

No inferred aliases adopted and no running clinical input is changed.
"""
import gzip
import hashlib
import json
import os
import urllib.request
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);out=root/'NSCLC_official_GEO_metadata_qualification_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    url='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243013/soft/GSE243013_family.soft.gz'
    plan={'created_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'resource_guard':guard,'official_url':url,'scope':'Official family SOFT metadata only; compressed cap2MiB, decompressed cap16MiB; ignore sample/platform data tables. No expression retrieval. Preserve all original metadata rows privately. Diagnose explicit patient/subject/donor fields and exact existing author-case tokens only; no approximate namespace bridges, no assumed patient aliases, no effect tests or changes to running frozen clinical family.','motivation':'Original PubMed abstract reports234post-neoadjuvant patients while native deposit has243authorrecords. Count difference alone neither proves duplicate patients nor resolves9NoICI/2ambiguous treatment records. Publisher fulltext unavailable; official supplement HEAD probes403 did not download any files. Seek exact original identity evidence rather than count-based inference.'}
    (out/'plan_before_official_metadata_only.json').write_text(json.dumps(plan,indent=2)+'\n')
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'public-GEO-metadata-qualification/1.0'}),timeout=35) as f:raw=f.read(2*1024**2+1)
        assert len(raw)<=2*1024**2 and raw.startswith(b'\x1f\x8b')
        import io
        with gzip.GzipFile(fileobj=io.BytesIO(raw),mode='rb') as f:data=f.read(16*1024**2+1)
        assert len(data)<=16*1024**2;records=[];current=None;table=False;series_fields={}
        for line in data.decode('utf-8').splitlines():
            if line.startswith('^'):
                if line.startswith('^SAMPLE = '):
                    current={'GSM':line.split(' = ',1)[1],'fields':{}};records.append(current)
                else:current=None
                table=False
            elif line.startswith('!sample_table_begin') or line.startswith('!platform_table_begin'):table=True
            elif line.startswith('!sample_table_end') or line.startswith('!platform_table_end'):table=False
            elif not table and ' = ' in line:
                key,value=line.split(' = ',1)
                if current is not None and key.startswith('!Sample_'):current['fields'].setdefault(key,[]).append(value)
                elif key.startswith('!Series_'):series_fields.setdefault(key,[]).append(value)
        assert records and len({r['GSM'] for r in records})==len(records)
        clinical=pd.read_csv(root/'NSCLC_immune_RNA_MCP_clinical_scope_20261019_v1/private_all243_current_RNA_clinical_scope.tsv',sep='\t',dtype=str,keep_default_na=False);native=set(clinical.sampleID);assert len(native)==243
        identity_keys={'patient','patient id','patient_id','patientid','subject','subject id','subject_id','donor','donor id','donor_id'};counts={};id_values={};intersections=[];schema={}
        for r in records:
            fields=r['fields'];exact=set();explicit={}
            for key,values in fields.items():
                schema[key]=schema.get(key,0)+1
                for value in values:
                    if value in native:exact.add(value)
                    if key.startswith('!Sample_characteristics') and ':' in value:
                        characteristic,text=value.split(':',1);characteristic=characteristic.strip().lower();text=text.strip()
                        if text in native:exact.add(text)
                        if characteristic in identity_keys:
                            explicit.setdefault(characteristic,[]).append(text);counts[characteristic]=counts.get(characteristic,0)+1;id_values.setdefault(characteristic,set()).add(text)
            intersections.append({'GSM':r['GSM'],'exact_native_author_tokens':sorted(exact),'explicit_original_identity_fields':explicit})
        payload=json.dumps({'series_fields':series_fields,'sample_records':records},indent=2)+'\n';matching=json.dumps(intersections,indent=2)+'\n';assert len(raw)+len(payload.encode())+len(matching.encode())<7*1024**2
        (out/'private_original_GSE243013_family.soft.gz').write_bytes(raw);(out/'private_all_original_family_metadata.json').write_text(payload);(out/'private_all_original_identity_intersections.json').write_text(matching)
        proof={'status':'COMPLETE_OFFICIAL_GEO_METADATA_SCHEMA_AND_EXACT_TOKEN_DIAGNOSTIC_NO_ALIAS_ADOPTED','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'official_url':url,'compressed_bytes':len(raw),'decompressed_bytes':len(data),'official_gzip_sha256':hashlib.sha256(raw).hexdigest(),'original_GEO_records':len(records),'explicit_identity_field_record_positions':counts,'explicit_identity_distinct_values':{k:len(v) for k,v in id_values.items()},'records_with_exact_native_author_token':sum(bool(r['exact_native_author_tokens']) for r in intersections),'unique_exact_native_author_tokens':len({x for r in intersections for x in r['exact_native_author_tokens']}),'metadata_field_schema_record_counts':schema,'adopted_patient_aliases':0,'running_clinical_input_changed':False,'files':{p.name:{'sha256':sha(p),'bytes':p.stat().st_size} for p in out.iterdir() if p.is_file()}}
    except Exception as e:
        proof={'status':'OFFICIAL_GEO_METADATA_UNAVAILABLE_NO_ALIAS_ADOPTED','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'error_class':type(e).__name__,'http_code':getattr(e,'code',None),'running_clinical_input_changed':False,'adopted_patient_aliases':0}
    (out/'private_terminal.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps({'status':proof['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
