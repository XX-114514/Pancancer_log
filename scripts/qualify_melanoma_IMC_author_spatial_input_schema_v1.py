"""Read only the author's explicit spatial-input headers, not graph rows.

No author code execution, coordinate values, cell IDs or learned matrices.
Directory/prefetch buffers remain RAM. Header-only is not full member CRC.
"""
from pathlib import Path
from datetime import datetime,timezone
import csv,hashlib,importlib.util,json,os,re,traceback,zipfile


def sha(b):return hashlib.sha256(b).hexdigest()


def load(name,path,expected):
    assert sha(path.read_bytes())==expected
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def main():
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert not out.exists();out.mkdir(mode=0o700)
    def save(n,v):(out/n).write_text(json.dumps(v,indent=2)+'\n')
    save('plan_before_author_spatial_input_headers.json',dict(frozen_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__).read_bytes()),
        scope='Fixed author Spatial_analysis.Rmd plus at most two explicitly read input CSV headers, max16KiB each.',
        no_cell_or_graph_rows_parsed_or_persisted=True,no_author_code_execution=True,
        generic_header_whitelist_only_coded_names_suppressed=True,no_physical_coordinate_or_graph_semantics_accepted_from_header_alone=True))
    helper=Path(os.environ['PANCANCER_IMC_TRANSPORT_HELPER'])
    source=load('published_archive_transport',helper,'9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2')
    transport=load('bounded_archive_transport',Path(os.environ['PANCANCER_IMC_BOUNDED_TRANSPORT']),'709f9b5ced875fea663557bad936f32d004e12d5daf291048b7de476064b0d68')
    d=json.loads(source.fetch(source.BASE+'api/records/6838169',2*1024**2));f=next(x for x in d['files'] if x['key']=='Melanoma_code_data.zip')
    assert f['size']==668227028 and f['checksum']=='md5:da74b7228763b76dcd2d9830c1a10cd9'
    rf=transport.RangeFile(source,f['links']['self'],f['size']);results=[]
    safe={'','X','Y','x','y','CellID','cell_id','cellID','CellType','celltype','type','cluster','Cluster','ROI','sample_id','sample','Sample',
        'ImageNumber','ObjectNumber','Position','CellClass','Group','class','ID','from','to','source','target','cell1','cell2',
        'ObjectNumber1','ObjectNumber2','neighbor','neighbour','neighbors','neighbours','distance','weight',
        'Location_Center_X','Location_Center_Y','Centroid_X','Centroid_Y','centroid_x','centroid_y'}
    with zipfile.ZipFile(rf) as z:
        b=z.read('Rscript/Spatial_analysis.Rmd');assert sha(b)=='d8faddad268f2900a83f8cfa3f8352dd697128d408a528181fb997c5773dba56'
        text=b.decode('utf-8-sig')
        for role in ['sc_data','n_dat']:
            line=next(s for s in text.splitlines() if re.match(r'\s*'+role+r'\s*<-\s*read.csv',s))
            match=re.search(r'''read.csv\(\s*(["'])([^"'\r\n]+)\1''',line);assert match
            literal=match.group(2)
            while literal.startswith('../') or literal.startswith('./'):literal=literal.split('/',1)[1]
            assert literal in z.namelist()
            entry=z.getinfo(literal)
            with z.open(entry) as stream:header=stream.readline(16385)
            assert len(header)<=16384 and header.endswith(b'\n')
            columns=next(csv.reader([header.decode('utf-8-sig').strip()]))
            results.append(dict(author_variable=role,source_member_path_SHA256=sha(literal.encode()),
                uncompressed_member_bytes=entry.file_size,compressed_member_bytes=entry.compress_size,
                header_SHA256=sha(header),source_column_count=len(columns),
                safe_generic_column_names=[c for c in columns if c in safe],
                other_unexported_column_count=sum(c not in safe for c in columns),
                explicit_coordinate_axis_names_present=any(c in columns for c in ['Location_Center_X','Centroid_X','centroid_x','x']) and any(c in columns for c in ['Location_Center_Y','Centroid_Y','centroid_y','y']),
                only_complete_header_parsed=True,full_member_CRC_not_verified=True))
    last=json.loads(source.fetch(source.BASE+'api/records/6838169',2*1024**2));other=next(x for x in last['files'] if x['key']==f['key'])
    assert all(f.get(k)==other.get(k) for k in ['id','size','checksum','links']);rf.close()
    save('private_complete_author_spatial_input_schema.json',dict(status='COMPLETE_AUTHOR_EXPLICIT_SPATIAL_INPUT_HEADERS_ONLY',
        completed_UTC=datetime.now(timezone.utc).isoformat(),producer_SHA256=sha(Path(__file__).read_bytes()),schemas=results,
        no_graph_coordinate_cell_rows_or_arrays_written=True,no_new_spatial_or_clinical_effects=True,
        no_false_full_member_or_coordinate_unit_or_raw_graph_acceptance=True))
    save('private_terminal.json',dict(status='COMPLETE_AUTHOR_EXPLICIT_SPATIAL_INPUT_HEADERS_ONLY',completed_UTC=datetime.now(timezone.utc).isoformat(),
        source_SHA256=sha(Path(__file__).read_bytes()),files={p.name:dict(SHA256=sha(p.read_bytes()),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_AUTHOR_EXPLICIT_SPATIAL_INPUT_HEADERS_ONLY',terminal_SHA256=sha((out/'private_terminal.json').read_bytes()))))


if __name__=='__main__':
    try:main()
    except Exception as error:
        out=Path(os.environ['PANCANCER_RAM_OUTPUT']);out.mkdir(mode=0o700,parents=True,exist_ok=True)
        value=dict(status='FAILED_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],no_record_values_in_error=True)
        (out/'private_failed_terminal.json').write_text(json.dumps(value,indent=2)+'\n')
        print(json.dumps(dict(status=value['status'],error_class=value['error_class'])))
        raise SystemExit(1)
