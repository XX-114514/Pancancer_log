"""Complete selected processed cell member in RAM; measurement, not effects.

Read one processed CSV, never raw images/full archive. No individual keys,
records, expression arrays or learned models are written. No protein-positive
threshold, state cluster, spatial association or clinical model is fitted.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import csv,hashlib,importlib.util,io,json,os,re,resource,traceback,zipfile
import numpy as np
import pandas as pd

PROTEINS=['CD45','CD3','CD4','CD8a','CD27','FoxP3','CD127','CD194','CD20','CD38','HLA-DR','CD74','CD68','CD14','CD16','CD11c','CD11b','IDO','Vimentin','SMA','E-cadherin','EpCAM','CAIX','VEGF','PDGFRb','FAP','AXL','Collagen-I','Ki-67','PDL1','OX40','PD1','LAG3','TIM3','ICOS']


def sha(b):return hashlib.sha256(b).hexdigest()


def load(name,path,expected):
    assert sha(path.read_bytes())==expected
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def main():
    helper=Path(os.environ['PANCANCER_IMC_TRANSPORT_HELPER'])
    source=load('published_archive_transport',helper,'9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2')
    retry=load('bounded_archive_transport',Path(os.environ['PANCANCER_IMC_BOUNDED_TRANSPORT']),'709f9b5ced875fea663557bad936f32d004e12d5daf291048b7de476064b0d68')
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert not out.exists();out.mkdir(mode=0o700)
    def save(n,v):(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    save('plan_before_complete_cell_member.json',dict(frozen_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__).read_bytes()),
        selected_member='data/input/sc_data.csv',declared_uncompressed_bytes=172707928,declared_compressed_bytes=51454934,
        scope='All rows of the single processed member in RAM; full CRC/SHA, protein finite/zero/variable measurement, author label names/counts and exact ROI-key coverage. No effects.',
        memory_budget_bytes=2*1024**3,transport_budget='Each HTTP range<=512KiB with default-TLS persistent connections; max3 transient tries, 45s alarm per request, check75min overall before each request; all requested bytes incl failed<=64MiB; no full archive fallback.',
        no_raw_images_sequences_or_learnt_models=True,no_individual_records_or_arrays_written=True,
        no_threshold_or_state_definition_or_patient_prefix_inference=True,no_new_clinical_or_spatial_tests=True))
    api=source.fetch(source.BASE+'api/records/6838169',2*1024**2);d=json.loads(api);f=next(x for x in d['files'] if x['key']=='Melanoma_code_data.zip')
    assert d['metadata']['access_right']=='open' and f['size']==668227028 and f['checksum']=='md5:da74b7228763b76dcd2d9830c1a10cd9'
    rf=retry.RangeFile(source,f['links']['self'],f['size'],checkpoint=lambda v:save('private_transport_progress.json',v))
    with zipfile.ZipFile(rf) as z:
        directory=sha(json.dumps([(e.filename,e.file_size,e.compress_size,e.CRC) for e in z.infolist()],separators=(',',':')).encode())
        assert directory=='26339da50b39ffdf7beb641a67f50523c867f2d3a86b24f42175f3a16dfe8992'
        rb=z.read('data/input/meta/ROI_meta_afterQC.csv')
        assert sha(rb)=='d0fa2438b7b91157a7848654fb1f23e40c9dbcca8020f7f5ff0566ce7acede21'
        registry=list(csv.DictReader(io.StringIO(rb.decode('utf-8-sig'))));roi_keys={r['sampleROI2'] for r in registry}
        info=z.getinfo('data/input/sc_data.csv');assert info.file_size==172707928 and info.compress_size==51454934
        buffer=io.BytesIO();digest=hashlib.sha256();total=0
        with z.open(info) as stream:
            while True:
                chunk=stream.read(4*1024**2)
                if not chunk:break
                total+=len(chunk);assert total<=info.file_size;digest.update(chunk);buffer.write(chunk)
        assert total==info.file_size
        buffer.seek(0);header=buffer.readline();assert sha(header)=='c1dbd75e2834a17585af71efa99b481d1cada905a8ee721c2b78b7e17f13d0ad'
        buffer.seek(0);columns=next(csv.reader([header.decode('utf-8-sig').strip()]));assert len(columns)==44
        used=PROTEINS+['Position','sample_id','celltype4','celltype20']
        assert all(x in columns for x in used)
        buffer.seek(0);df=pd.read_csv(buffer,usecols=used,dtype={**{p:np.float64 for p in PROTEINS},**{p:str for p in used if p not in PROTEINS}},keep_default_na=False,float_precision='round_trip')
        rows=len(df);assert rows>0
        # All selected marker values independently parsed again by stdlib CSV/float.
        scalar_count=0;scalar_zero=np.zeros(len(PROTEINS),np.int64);scalar_negative=np.zeros(len(PROTEINS),np.int64)
        scalar_min=np.full(len(PROTEINS),np.inf);scalar_max=np.full(len(PROTEINS),-np.inf)
        scalar_ROI=Counter();scalar_types={x:Counter() for x in ['celltype4','celltype20']};scalar_format=Counter()
        number=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
        scalar_number=re.compile('^'+number+'$');pair=re.compile(r'^\(?\s*('+number+r')\s*,\s*('+number+r')\s*\)?$')
        values=df[PROTEINS].to_numpy()
        buffer.seek(0)
        for row in csv.DictReader(io.TextIOWrapper(buffer,encoding='utf-8-sig',newline='')):
            v=np.array([float(row[p]) for p in PROTEINS]);assert np.isfinite(v).all()
            scalar_zero+=(v==0);scalar_negative+=(v<0);scalar_min=np.minimum(scalar_min,v);scalar_max=np.maximum(scalar_max,v)
            scalar_ROI[row['sample_id']]+=1
            for c in scalar_types:scalar_types[c][row[c]]+=1
            value=row['Position'].strip()
            scalar_format['numeric_scalar' if scalar_number.fullmatch(value) else 'numeric_pair_string' if pair.fullmatch(value) else 'other_or_empty']+=1
            assert np.array_equal(v,values[scalar_count])
            scalar_count+=1
        assert scalar_count==rows
        assert scalar_ROI==df['sample_id'].value_counts().to_dict()
        for c in scalar_types:assert scalar_types[c]==df[c].value_counts().to_dict()
        values=df[PROTEINS].to_numpy();assert np.isfinite(values).all()
        assert np.array_equal(scalar_zero,(values==0).sum(axis=0)) and np.array_equal(scalar_negative,(values<0).sum(axis=0))
        assert np.array_equal(scalar_min,values.min(axis=0)) and np.array_equal(scalar_max,values.max(axis=0))
        marker_summary=[dict(marker=p,cells=rows,zeros=int(scalar_zero[i]),negative=int(scalar_negative[i]),minimum=float(scalar_min[i]),maximum=float(scalar_max[i]),
            quantiles=np.quantile(values[:,i],[0,.25,.5,.75,1]).tolist(),has_variation=bool(scalar_min[i]<scalar_max[i])) for i,p in enumerate(PROTEINS)]
        # Only recognized biological names or numeric cluster codes may leave RAM.
        safe_label=re.compile(r'^(?:[A-Za-z][A-Za-z0-9 +()_./-]{0,55}|\d{1,3})$')
        labels={c:dict(counts) if all(safe_label.fullmatch(v) for v in counts) else None for c,counts in scalar_types.items()}
        cell_keys=set(scalar_ROI);overlap=cell_keys&roi_keys
        measurement=dict(status='COMPLETE_SELECTED_PROCESSED_CELL_MEMBER_MEASUREMENT_QUALIFICATION',completed_UTC=datetime.now(timezone.utc).isoformat(),
            producer_SHA256=sha(Path(__file__).read_bytes()),complete_member_SHA256=digest.hexdigest(),complete_member_CRC32=info.CRC,complete_member_CRC_checked=True,
            uncompressed_bytes=total,cells=rows,protein_axis=PROTEINS,protein_measurements=marker_summary,author_label_counts=labels,
            source_ROI_key_count=len(roi_keys),cell_table_source_key_count=len(cell_keys),literal_exact_overlap_keys=len(overlap),
            cells_with_literal_ROI_key=sum(scalar_ROI[k] for k in overlap),all_literal_ROI_keys_exact=cell_keys==roi_keys,
            Position_format_counts=dict(scalar_format),Position_not_assumed_2D_centroid_or_um=True,
            complete_independent_CSV_numeric_ROI_and_label_counter_replay=True,
            protein_positive_flags_or_decoy_FDR_thresholds_not_present_in_selected_schema=True,
            no_greater_than_zero_positivity_rule_adopted=True,no_patient_clinical_or_spatial_effects=True,
            no_individual_records_arrays_or_keys_written=True,full_archive_MD5_not_claimed=True)
    with zipfile.ZipFile(rf) as z:
        assert sha(json.dumps([(e.filename,e.file_size,e.compress_size,e.CRC) for e in z.infolist()],separators=(',',':')).encode())==directory
    final=json.loads(source.fetch(source.BASE+'api/records/6838169',2*1024**2));other=next(x for x in final['files'] if x['key']==f['key'])
    assert all(f.get(k)==other.get(k) for k in ['id','size','checksum','links'])
    peak_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*1024
    assert peak_bytes<=2*1024**3
    rf.close()
    measurement.update(published_metadata_and_directory_replayed=True,range_bytes=rf.transferred,requested_bytes_including_failed=rf.requested,transport_failure_classes=rf.failures,
        actual_Linux_peak_RSS_bytes=peak_bytes)
    save('private_complete_cell_measurement_qualification.json',measurement)
    save('private_terminal.json',dict(status=measurement['status'],completed_UTC=datetime.now(timezone.utc).isoformat(),producer_SHA256=sha(Path(__file__).read_bytes()),
        files={p.name:dict(SHA256=sha(p.read_bytes()),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status=measurement['status'],terminal_SHA256=sha((out/'private_terminal.json').read_bytes()))))


if __name__=='__main__':
    try:main()
    except Exception as error:
        out=Path(os.environ['PANCANCER_RAM_OUTPUT']);out.mkdir(mode=0o700,parents=True,exist_ok=True)
        value=dict(status='FAILED_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],no_individual_values_in_error_record=True)
        (out/'private_failed_terminal.json').write_text(json.dumps(value,indent=2)+'\n')
        print(json.dumps(dict(status=value['status'],error_class=value['error_class'])))
        raise SystemExit(1)
