"""Bounded metadata/panel/ZIP schema qualification, no clinical effects.

The archive directory and selected code/header bytes stay in RAM. No images,
individual records or learned matrices are saved. Author code is never run.
"""
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import csv
import hashlib
import io
import json
import os
import re
import urllib.request
import zipfile

BASE='https://zenodo.org/'
RECORD=6838169


def sha(b):return hashlib.sha256(b).hexdigest()


def opener():return urllib.request.build_opener(urllib.request.ProxyHandler({}))


def fetch(url,limit):
    assert url.startswith(BASE)
    with opener().open(url,timeout=35) as response:
        assert response.status==200
        raw=response.read(limit+1);assert len(raw)<=limit
        return raw


class RangeFile:
    def __init__(self,url,size):
        assert url.startswith(BASE)
        self.url=url;self.size=size;self.position=0;self.transferred=0
        self.requests=0;self.identity=None;self.directory_reads=[]
    def tell(self):return self.position
    def seek(self,offset,whence=0):
        self.position=offset if whence==0 else self.position+offset if whence==1 else self.size+offset
        assert 0<=self.position<=self.size
        return self.position
    def seekable(self):return True
    def read(self,n=-1):
        if n<0:n=self.size-self.position
        n=min(n,self.size-self.position)
        if n==0:return b''
        assert n<=16*1024**2 and self.transferred+n<=64*1024**2
        start=self.position;end=start+n-1
        request=urllib.request.Request(self.url,headers={'Range':f'bytes={start}-{end}'})
        with opener().open(request,timeout=35) as response:
            assert response.status==206,'Range_not_supported_no_full_payload_fallback'
            assert response.headers['Content-Range']==f'bytes {start}-{end}/{self.size}'
            identity=(response.headers.get('ETag'),response.headers.get('Last-Modified'))
            assert any(identity)
            if self.identity is None:self.identity=identity
            else:assert self.identity==identity
            raw=response.read(n+1);assert len(raw)==n
        self.position+=n;self.transferred+=n;self.requests+=1
        return raw


def main():
    out=Path(os.environ['PANCANCER_PUBLIC_PREFLIGHT_OUTPUT']);assert not out.exists();out.mkdir(mode=0o700)
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    save('plan_before_any_payload.json',dict(frozen_UTC=datetime.now(timezone.utc).isoformat(),
        source_SHA256=sha(Path(__file__).read_bytes()),record=RECORD,
        scope='Official API/panel, complete ZIP central directory, generic author code and selected CSV headers only.',
        budget='API<=2MiB, panel<=1MiB, each range<=16MiB, total selected bytes<=64MiB. No full archive fallback.',
        no_author_code_execution=True,no_individual_rows_images_RNA_deconvolution_or_learned_weights_selected=True,
        no_new_clinical_or_spatial_tests=True))
    raw=fetch(BASE+f'api/records/{RECORD}',2*1024**2);record=json.loads(raw)
    assert record['id']==RECORD and record['metadata']['access_right']=='open'
    assert record['metadata']['license']['id']=='cc-by-4.0'
    files={f['key']:f for f in record['files']}
    assert set(files)=={'Melanoma_panel.csv','Melanoma_code_data.zip','melanoma_IMC_data.zip'}
    panel=files['Melanoma_panel.csv'];b=fetch(panel['links']['self'],1024**2)
    assert len(b)==panel['size'] and hashlib.md5(b).hexdigest()==panel['checksum'].split(':')[1]
    rows=list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))))
    markers=[r['marker'] for r in rows];assert len(markers)==len(set(markers))
    save('private_official_panel_and_file_contract.json',dict(API_SHA256=sha(raw),panel_SHA256=sha(b),
        source_access='open',license=record['metadata']['license']['id'],source_markers=markers,
        files={n:dict(size=f['size'],checksum=f['checksum']) for n,f in files.items()},
        required_protein_axis_presence={g:g in markers for g in ['CD3','CD4','CD8a','FoxP3','Ki67','CD223','CD279','CD366','CD127','CD137','CD25','CD152','TCF1','TOX']},
        original_RNA_recipe_not_assumed_equivalent_to_protein_panel=True))
    item=files['Melanoma_code_data.zip'];assert item['size']==668227028
    rf=RangeFile(item['links']['self'],item['size'])
    with zipfile.ZipFile(rf) as archive:
        entries=archive.infolist();names=[e.filename for e in entries]
        assert len(names)==len(set(names))
        assert all(not n.startswith('/') and '..' not in Path(n).parts for n in names)
        signature=sha(json.dumps([(e.filename,e.file_size,e.compress_size,e.CRC) for e in entries],separators=(',',':')).encode())
        # Only generic file names are exported; coded image/case names remain RAM.
        generic=[e for e in entries if not e.is_dir() and e.filename.endswith(('.Rmd','.R'))]
        generic_names=[e.filename for e in generic]
        safe_generic=[n for n in generic_names if re.fullmatch(r'[A-Za-z_./ -]+\.(?:R|Rmd)',n)]
        header_candidates=[e for e in entries if e.filename.endswith('/input/sc_data.csv') or e.filename=='input/sc_data.csv']
        directory=dict(total_entries=len(entries),directory_identity_SHA256=signature,
            file_extensions=dict(Counter(Path(e.filename).suffix.lower() for e in entries if not e.is_dir())),
            generic_code_paths=safe_generic,sc_data_header_opportunities=len(header_candidates),
            full_archive_bytes=item['size'],provided_archive_checksum=item['checksum'],full_archive_checksum_not_verified=True)
        save('private_complete_directory_schema.json',directory)
        schemas=[]
        allowed=set(markers)|{'x','y','X','Y','sample','Sample','sample_id','SampleID','roi','ROI','patient','Patient','response','Response','celltype','cell_type','CellType','cluster','Cluster','cellID','cell_id','major','subtype','type','Type','anno','annotation','group','Group','region','Region','IM','CT','Pos_X','Pos_Y','X_position','Y_position'}
        for entry in header_candidates:
            assert entry.file_size<=4*1024**3 and entry.compress_size<=2*1024**3
            with archive.open(entry) as stream:header=stream.readline(65537)
            assert len(header)<=65536 and header.endswith(b'\n')
            columns=next(csv.reader([header.decode('utf-8-sig').strip()]))
            schemas.append(dict(logical_role='author_sc_data',header_SHA256=sha(header),column_count=len(columns),
                safe_whitelisted_column_names=[s for s in columns if s in allowed],
                measured_channel_column_presence={row['marker']:bool(set(columns)&{row['marker'],row['label'],row['Metal'],row['Metal']+'Di'}) for row in rows},
                declared_member_bytes=entry.file_size,declared_compressed_member_bytes=entry.compress_size,
                only_header_parsed_no_clinical_or_cell_rows=True,full_member_CRC_not_verified=True))
        save('private_selected_header_schema.json',schemas)
        code_summaries=[]
        for entry in generic:
            if Path(entry.filename).name not in ['functions.R','Clustering.Rmd','Spatial_analysis.Rmd','generate_Figs.Rmd']:continue
            assert entry.file_size<=2*1024**2
            code=archive.read(entry);assert len(code)==entry.file_size
            text=code.decode('utf-8-sig');fields=set(re.findall(r'\$([A-Za-z_][A-Za-z0-9_.]*)',text))
            code_summaries.append(dict(generic_file=Path(entry.filename).name,bytes=len(code),SHA256=sha(code),full_member_CRC_checked_by_zipfile=True,
                safe_field_names=sorted(fields&allowed),
                keywords={k:k.lower() in text.lower() for k in ['positive','decoy','response','clinical','coord','centroid','Ki67','FoxP3','norm','FDR','phenograph']},
                author_code_not_executed=True,no_individual_literals_saved=True))
        save('private_generic_code_schema_proofs.json',code_summaries)
    save('private_terminal.json',dict(status='COMPLETE_OFFICIAL_PANEL_DIRECTORY_AND_SELECTED_HEADER_CODE_QUALIFICATION_ONLY',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__).read_bytes()),
        range_requests=rf.requests,range_bytes=rf.transferred,range_identity=rf.identity,
        no_full_archive_or_raw_images_read=True,no_clinical_or_cell_rows_parsed=True,no_effects_or_identity_assigned=True,
        no_individual_or_learned_arrays_saved=True,
        files={p.name:dict(SHA256=sha(p.read_bytes()),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_OFFICIAL_PANEL_DIRECTORY_AND_SELECTED_HEADER_CODE_QUALIFICATION_ONLY',terminal_SHA256=sha((out/'private_terminal.json').read_bytes()))))


if __name__=='__main__':main()
