"""Public GEO processed inventory v2: accept numbered supplementary fields.

No controlled sequencing reads, clonotype strings, private rows or effect tests.
"""
import gzip
import hashlib
import json
import os
import re
from datetime import datetime,timezone
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import urlparse

from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_native_MFA_interpretation_inputs_v1 import sha

URL='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE328nnn/GSE328481/soft/GSE328481_family.soft.gz'
OUT='GSE328481_public_processed_spatial_metadata_20261018_v2'


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    out=root/OUT;assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        question='Assess independent LUAD spatial measurement/annotation/scale and sample identity feasibility for fixed immune CM functional axes and supported ligand-state candidates. Metadata only, no expression or effect selection.',
        source=URL,source_access='Official public GEO processed matrix metadata only. Do not access controlled GSA FASTQ, send messages or retrieve clonotype sequences.',
        maximum_compressed_bytes=8*1024**2,maximum_uncompressed_bytes=32*1024**2,
        sample_scope='All samples in series; all supplementary links preserved. No gene/effect-based source choice.',
        future_expression='Not authorized by this script; require explicit bounded versioned schema/measurement plan and actual byte quota before download. No assumptions about raw counts from archive names.')
    (out/'plan_before_public_metadata_request.json').write_text(json.dumps(plan,indent=2)+'\n')
    req=Request(URL,headers={'User-Agent':'PanCancer-public-metadata-research/1.0'})
    with urlopen(req,timeout=40) as response:
        content=response.read(8*1024**2+1);assert len(content)<=8*1024**2
        transport=dict(final_url=response.url,status=response.status,content_length=response.headers.get('Content-Length'),last_modified=response.headers.get('Last-Modified'))
    with gzip.GzipFile(fileobj=__import__('io').BytesIO(content)) as stream:
        raw=stream.read(32*1024**2+1);assert len(raw)<=32*1024**2
    assert raw.startswith(b'^SERIES') or b'^SERIES = GSE328481' in raw
    text=raw.decode('utf-8');(out/'public_series_family.soft.gz').write_bytes(content)
    samples=[];series={};current=series
    permitted_series={'Series_title','Series_summary','Series_overall_design','Series_supplementary_file','Series_pubmed_id','Series_sample_id','Series_relation'}
    permitted_sample={'Sample_title','Sample_source_name_ch1','Sample_organism_ch1','Sample_characteristics_ch1','Sample_data_processing','Sample_supplementary_file','Sample_description','Sample_library_strategy','Sample_type','Sample_relation'}
    for line in text.splitlines():
        if line.startswith('^SAMPLE = '):
            current={'GSM':line.split(' = ',1)[1]};samples.append(current)
        elif line.startswith('^'):current=series if line.startswith('^SERIES = ') else {}
        elif line.startswith('!') and ' = ' in line:
            key,value=line[1:].split(' = ',1)
            key=re.sub(r'^(Sample|Series)_supplementary_file_\d+$',r'\1_supplementary_file',key)
            if (key in permitted_series and current is series) or (key in permitted_sample and current is not series):current.setdefault(key,[]).append(value)
    assert samples and len({s['GSM'] for s in samples})==len(samples)
    # Only supplementary processed-file HEAD requests; never controlled-access relation URLs.
    assert any(s.get('Sample_supplementary_file') for s in samples), 'No sample links parsed; inspect actual SOFT field names'
    files=[]
    for sample in samples:
        for url in sample.get('Sample_supplementary_file',[]):
            if url=='NONE':continue
            parsed=urlparse(url);assert parsed.hostname=='ftp.ncbi.nlm.nih.gov'
            endpoint='https://'+parsed.netloc+parsed.path
            row=dict(GSM=sample['GSM'],url=endpoint,filename=Path(parsed.path).name)
            try:
                with urlopen(Request(endpoint,method='HEAD'),timeout=25) as response:
                    row.update(status=response.status,content_length=response.headers.get('Content-Length'),last_modified=response.headers.get('Last-Modified'),etag=response.headers.get('ETag'))
            except Exception as e:row.update(status='HEAD_FAILED',reason=type(e).__name__+': '+str(e))
            files.append(row)
    (out/'private_all_series_sample_metadata.json').write_text(json.dumps(dict(series=series,samples=samples),indent=2)+'\n')
    (out/'public_all_processed_file_inventory.json').write_text(json.dumps(files,indent=2)+'\n')
    outputs=[p for p in out.iterdir() if p.is_file()]
    terminal=dict(status='COMPLETE_PUBLIC_GEO_METADATA_ONLY_EXPRESSION_SCHEMA_UNVERIFIED',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),transport=transport,source_gzip_SHA256=hashlib.sha256(content).hexdigest(),samples=len(samples),processed_links=len(files),controlled_or_expression_downloads=False,
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in outputs})
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],completed_UTC=terminal['completed_UTC'],samples=len(samples),processed_links=len(files),terminal_sha256=sha(out/'private_terminal.json'))))

if __name__=='__main__':main()
