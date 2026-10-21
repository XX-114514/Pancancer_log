"""Qualify explicit imaging patient/ROI keys; ignore learned archetype values.

Complete small CSVs and generic author code are held in RAM. Export only
schema, cross-source agreement counts and masked code; no record keys.
"""
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter,defaultdict
import csv,hashlib,importlib.util,io,json,os,re,traceback,zipfile
import pandas as pd


def sha(b):return hashlib.sha256(b).hexdigest()


def load(name,path,expected):
    assert sha(path.read_bytes())==expected
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def main():
    helper=Path(os.environ['PANCANCER_IMC_TRANSPORT_HELPER'])
    source=load('published_archive_transport',helper,'9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2')
    retry=load('bounded_archive_transport',helper.parent.parent/'pancancer-CM02-melanoma-ROI-registry-20261021-v2'/'qualify_melanoma_IMC_ROI_registry_v2.py','190b81b9c771d06e1b3b63506a77c13ef29202f96553b456973796cec5062910')
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert not out.exists();out.mkdir(mode=0o700)
    def save(n,v):(out/n).write_text(json.dumps(v,indent=2)+'\n')
    save('plan_before_imaging_key_rows.json',dict(frozen_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__).read_bytes()),
        selected_complete_metadata=['data/input/meta/ROI_meta_afterQC.csv','data/input/hc_result/hc_result_sample_anno_invasive.csv'],
        only_explicit_patient_ROI_region_response_fields_used=True,archetype_values_ignored=True,
        no_prefix_identity_inference=True,no_record_keys_or_rows_or_arrays_written=True,no_new_effects=True,
        schema_header_safety='Only declared generic clinical fields exported; noncontract coded column headers suppressed.',
        complete_bridge_subset_allowed='An exact subset of registry ROIs is a valid partial bridge; uncovered regions remain unmapped.'))
    d=json.loads(source.fetch(source.BASE+'api/records/6838169',2*1024**2));f=next(x for x in d['files'] if x['key']=='Melanoma_code_data.zip')
    assert f['size']==668227028 and f['checksum']=='md5:da74b7228763b76dcd2d9830c1a10cd9'
    rf=retry.BoundedRetryRangeFile(source,f['links']['self'],f['size']);tables={};proofs=[];code_paths=[];contexts=[]
    with zipfile.ZipFile(rf) as z:
        for path in ['data/input/meta/ROI_meta_afterQC.csv','data/input/hc_result/hc_result_sample_anno_invasive.csv']:
            info=z.getinfo(path);assert info.file_size<=65536;b=z.read(info)
            rows=list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))));cols=list(rows[0])
            # Other columns may carry learned archetype labels: do not inspect or export their values.
            df=pd.read_csv(io.BytesIO(b),dtype=str,keep_default_na=False)
            used=[c for c in cols if c in ['', 'sampleROI2','sampleROI','sample_id','ROI','roi','core','sample','patient','patient_id','patientID','region','response']]
            for c in used:assert df.iloc[:,cols.index(c)].tolist()==[r[c] for r in rows]
            selected=[{c:r[c] for c in used} for r in rows]
            tables[path]=selected
            proofs.append(dict(generic_path=path,complete_SHA256=sha(b),CRC_checked=True,rows=len(rows),safe_schema_columns=[c for c in cols if c in set(used+['QUPATH','TLS','response','type','TNMstage','Site','QC'])],
                unselected_noncontract_column_count=len([c for c in cols if c not in set(used+['QUPATH','TLS','response','type','TNMstage','Site','QC'])]),used_columns=used,
                exact_selected_csv_Pandas_value_replay=True,no_archetype_values_selected=True))
        # Generic file paths are source code, not any data record. Reject coded file names.
        for path in ['Rscript/Spatial_analysis.Rmd','Rscript/functions.R','Rscript/generate_Figs.Rmd']:
            b=z.read(path);text=b.decode('utf-8-sig');lines=text.splitlines()
            for m in re.finditer(r'''["']((?:\.\./|\./)*(?:data|Rscript)/[A-Za-z_./ -]+\.(?:csv|R))["']''',text):
                normalized=m.group(1)
                while normalized.startswith('../') or normalized.startswith('./'):normalized=normalized.split('/',1)[1]
                if normalized in z.namelist():code_paths.append(dict(source=path,generic_path=normalized))
            targets=set()
            for i,line in enumerate(lines):
                if re.search(r'renamePID\s*(?:<-|=)\s*function|n_dat\s*<-\s*read|ROI.meta\s*<-\s*read|anno3\s*<-\s*read|colnames\(sc_data\)|colnames\(n_dat\)',line):
                    targets.update(range(max(0,i-1),min(len(lines),i+9)))
            safe={'patient','patient_id','sample','sample_id','sampleROI2','region','response','ID','CellID','cell_id','x','y','X','Y','core','ROI','roi','_','-'}
            quote=re.compile(r'''(["'])(?:\\.|(?!\1).)*?\1''')
            for i in sorted(targets):
                line=quote.sub(lambda m:m.group() if m.group()[1:-1] in safe else '"[masked_literal]"',lines[i])
                contexts.append(dict(source=path,line=i+1,masked_generic_code=line[:600]))
    roi=tables['data/input/meta/ROI_meta_afterQC.csv'];bridge=tables['data/input/hc_result/hc_result_sample_anno_invasive.csv']
    roi_index={r['sampleROI2']:r for r in roi};assert len(roi_index)==len(roi)
    patient_fields=[c for c in bridge[0] if c.lower().replace('_','') in ['patient','patientid']]
    joins=[]
    for field in bridge[0]:
        keys=[r[field] for r in bridge]
        if len(set(keys))!=len(keys) or not set(keys)<=set(roi_index):continue
        for patient in patient_fields:
            groups=defaultdict(list)
            for row in bridge:groups[row[patient]].append(roi_index[row[field]])
            responses={p:{r['response'] for r in rs} for p,rs in groups.items()}
            consistent=all(len(v)==1 for v in responses.values())
            counts=Counter(next(iter(v)) for v in responses.values()) if consistent else Counter()
            assert len(groups)==pd.DataFrame(bridge)[patient].nunique()
            joins.append(dict(bridge_ROI_field=field,explicit_patient_field=patient,exact_all_bridge_ROIs_in_registry=True,
                full_registry_covered=set(keys)==set(roi_index),
                bridge_ROIs=len(keys),regions_of_matched_ROIs=dict(Counter(roi_index[k]['region'] for k in keys)),
                direct_response_field_disagreements=sum(row.get('response')!=roi_index[row[field]]['response'] for row in bridge),
                source_patient_units=len(groups),all_units_nonempty=all(groups),all_responses_consistent=consistent,
                source_units_by_response=dict(counts) if set(counts)<=set(['R','NR']) else None,
                matches_published26_and14R12NR=consistent and len(groups)==26 and counts=={'R':14,'NR':12},
                no_prefix_rule_or_outcome_inference_used=True))
    last=json.loads(source.fetch(source.BASE+'api/records/6838169',2*1024**2));other=next(x for x in last['files'] if x['key']==f['key'])
    assert all(f.get(k)==other.get(k) for k in ['id','size','checksum','links'])
    save('private_complete_imaging_key_qualification.json',dict(status='COMPLETE_EXPLICIT_IMAGING_KEY_QUALIFICATION',completed_UTC=datetime.now(timezone.utc).isoformat(),
        source_SHA256=sha(Path(__file__).read_bytes()),metadata_proofs=proofs,exact_explicit_patient_ROI_joins=joins,
        generic_existing_code_paths=code_paths,masked_generic_code_contexts=contexts,
        no_patient_prefix_inference=True,no_cell_rows_or_spatial_or_clinical_effects_read=True,
        range_bytes=rf.transferred,requested_bytes_including_failed=rf.requested,transport_failure_classes=rf.failures))
    save('private_terminal.json',dict(status='COMPLETE_EXPLICIT_IMAGING_KEY_QUALIFICATION',completed_UTC=datetime.now(timezone.utc).isoformat(),
        source_SHA256=sha(Path(__file__).read_bytes()),files={p.name:dict(SHA256=sha(p.read_bytes()),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status='COMPLETE_EXPLICIT_IMAGING_KEY_QUALIFICATION',terminal_SHA256=sha((out/'private_terminal.json').read_bytes()))))


if __name__=='__main__':
    try:main()
    except Exception as error:
        out=Path(os.environ['PANCANCER_RAM_OUTPUT']);out.mkdir(mode=0o700,parents=True,exist_ok=True)
        value=dict(status='FAILED_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],no_individual_values_in_error_record=True)
        (out/'private_failed_terminal.json').write_text(json.dumps(value,indent=2)+'\n')
        print(json.dumps(dict(status=value['status'],error_class=value['error_class'])))
        raise SystemExit(1)
