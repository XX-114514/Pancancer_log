"""Outcome-blind patient/type/marker measurement coverage, RAM only.

Requires a completed, SHA-checked cell-measurement qualification first. It is
not an association test, CM projection or new cell-state annotation. Multiple
ROIs are pooled within the explicit source patient. No IDs/rows/arrays saved.
"""
from collections import Counter,defaultdict
from datetime import datetime,timezone
from itertools import combinations
from pathlib import Path
import csv,hashlib,importlib.util,io,json,os,re,resource,sys,traceback,zipfile
import numpy as np
import pandas as pd

MARKERS=['CD3','CD4','CD8a','FoxP3','Ki-67','LAG3','PD1','TIM3','OX40','CD127']
PARENT_CELL_MINIMUM=20


def sha(b):return hashlib.sha256(b).hexdigest()


def load(name,path,expected):
    assert sha(path.read_bytes())==expected
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def explicit_margin_map(registry,bridge):
    """No prefix rules; ignore outcome/archetype/coded wide-column values."""
    index={r['sampleROI2']:r['region'] for r in registry}
    assert len(index)==len(registry) and all(index)
    mapping={}
    for row in bridge:
        key=row['ROI'];patient=row['patient']
        assert key not in mapping and key in index and patient
        assert index[key]=='invasive'
        mapping[key]=patient
    assert mapping and set(mapping)=={k for k,v in index.items() if v=='invasive'}
    return mapping,set(index)


def summarize(frame,roi_to_patient,registry_keys):
    """All-author-type opportunities; zeros are counts, absent means stay NaN."""
    assert set(frame['sample_id'])==registry_keys
    assert np.isfinite(frame[MARKERS].to_numpy()).all()
    types=sorted(set(frame['celltype20']));patients=sorted(set(roi_to_patient.values()))
    assert 1<=len(types)<=64 and all(isinstance(x,str) and x for x in types)
    selected=frame['sample_id'].isin(roi_to_patient)
    dat=frame.loc[selected,['sample_id','celltype20']+MARKERS].copy()
    dat['source_patient']=dat['sample_id'].map(roi_to_patient)
    assert dat['source_patient'].notna().all()
    cells=dat.groupby(['source_patient','celltype20'],sort=True).size().unstack(fill_value=0).reindex(index=patients,columns=types,fill_value=0)
    means=dat.groupby(['source_patient','celltype20'],sort=True)[MARKERS].mean()
    roi_means=dat.groupby(['source_patient','celltype20','sample_id'],sort=True)[MARKERS].mean()
    balanced=roi_means.groupby(level=[0,1]).mean()
    assert int(cells.to_numpy().sum())==len(dat)
    cell_profiles=[];marker_profiles=[]
    safe=re.compile(r'^(?:[A-Za-z][A-Za-z0-9 +()_./-]{0,55}|\d{1,3})$')
    labels={x:x if safe.fullmatch(x) else 'unexported_author_label_'+str(i+1) for i,x in enumerate(types)}
    for typ in types:
        count=cells[typ].to_numpy();eligible=count>=PARENT_CELL_MINIMUM
        cell_profiles.append(dict(author_type=labels[typ],source_patient_units=len(patients),
            source_cells=int(count.sum()),patients_with_any=int((count>0).sum()),patients_ge20=int(eligible.sum()),
            patients_ge50=int((count>=50).sum()),patients_without_type=int((count==0).sum()),
            all_source_patients_preserved=True,author_identity_not_new_formal_annotation=True))
        for marker in MARKERS:
            for representation,table in [('cell_pooled',means),('ROI_equal_weighted',balanced)]:
                series=table[marker].unstack().reindex(index=patients,columns=types)[typ]
                value=series.to_numpy(dtype=float);assert np.isnan(value[count==0]).all()
                take=value[eligible];assert np.isfinite(take).all()
                marker_profiles.append(dict(author_type=labels[typ],marker=marker,representation=representation,
                    patients_parent_ge20=len(take),distinct_patient_means=len(np.unique(take)),
                    variable_between_patients=bool(len(take)>=2 and np.ptp(take)>0),
                    mean_IQR_quartiles=np.quantile(take,[.25,.5,.75]).tolist() if len(take)>=10 else None,
                    small_group_values_suppressed=len(take)<10,not_positive_cell_fraction=True))
    pair_profiles=[]
    for a,b in combinations(types,2):
        n=int(((cells[a]>=PARENT_CELL_MINIMUM)&(cells[b]>=PARENT_CELL_MINIMUM)).sum())
        pair_profiles.append(dict(author_type_A=labels[a],author_type_B=labels[b],
            source_patient_units=len(patients),both_types_ge20_patients=n,at_least20_joint_patients=n>=20,
            coverage_only_not_full_model_estimability=True))
    return dict(author_type_coverage=cell_profiles,all_author_type_pair_coverage=pair_profiles,
        all_marker_mean_coverage=marker_profiles,source_patient_units=len(patients),
        selected_margin_cells=len(dat),known_nonmargin_cells=int((~selected).sum()),
        source_ROIs_with_explicit_patient=len(roi_to_patient),
        zero_type_counts_and_absent_marker_means_kept_distinct=True,
        outcome_and_archetype_values_not_used=True,no_new_clinical_or_spatial_effects=True),cells,means


def scalar_replay(blob,roi_to_patient,cells,means):
    counter=Counter();sums=defaultdict(lambda:np.zeros(len(MARKERS),dtype=float));n=0
    for row in csv.DictReader(io.StringIO(blob.decode('utf-8-sig'))):
        if row['sample_id'] not in roi_to_patient:continue
        key=(roi_to_patient[row['sample_id']],row['celltype20'])
        value=np.array([float(row[m]) for m in MARKERS]);assert np.isfinite(value).all()
        counter[key]+=1;sums[key]+=value;n+=1
    for patient in cells.index:
        for typ in cells.columns:assert int(cells.loc[patient,typ])==counter[(patient,typ)]
    for key,count in counter.items():
        assert np.allclose(sums[key]/count,means.loc[key,MARKERS].to_numpy(dtype=float),rtol=1e-10,atol=1e-10)
    assert n==int(cells.to_numpy().sum())
    return dict(all_CSV_dictionary_patient_type_counts_and_marker_means_replayed=True,
        independent_reduction_not_independent_input_bytes=True)


def main():
    # A partial file, progress record or process lock cannot satisfy this gate.
    accepted=Path(os.environ['PANCANCER_IMC_ACCEPTED_MEASUREMENT'])
    terminal=accepted.parent/'private_terminal.json';t=json.loads(terminal.read_text());prior=json.loads(accepted.read_text())
    assert t['status']==prior['status']=='COMPLETE_SELECTED_PROCESSED_CELL_MEMBER_MEASUREMENT_QUALIFICATION'
    assert t['producer_SHA256']=='eb142bcdeb172f068f9d182a6ba89866d3d07e852c75f0b976ec4c5b43cd436d'
    assert t['files'][accepted.name]['SHA256']==sha(accepted.read_bytes())
    assert prior['complete_member_CRC_checked'] and prior['complete_independent_CSV_numeric_ROI_and_label_counter_replay']
    assert prior['published_metadata_and_directory_replayed'] and prior['all_literal_ROI_keys_exact']
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert not out.exists();out.mkdir(mode=0o700)
    def save(n,v):(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    save('plan_before_patient_coverage.json',dict(frozen_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__).read_bytes()),
        accepted_measurement_SHA256=sha(accepted.read_bytes()),scope='All author celltype20 labels, all distinct pairs and ten fixed marker means; outcome-blind eligibility only.',
        markers=MARKERS,parent_cell_minimum=PARENT_CELL_MINIMUM,source_patient_unit='Explicit patient field and invasive-margin literal ROI map only.',
        no_ROI_pseudoreplication_or_prefix_inference=True,no_protein_positive_threshold_or_clinical_tests=True,
        no_individual_keys_rows_expression_or_learned_arrays_written=True))
    helper=Path(os.environ['PANCANCER_IMC_TRANSPORT_HELPER'])
    source=load('published_archive_transport',helper,'9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2')
    transport=load('bounded_archive_transport',Path(os.environ['PANCANCER_IMC_BOUNDED_TRANSPORT']),'709f9b5ced875fea663557bad936f32d004e12d5daf291048b7de476064b0d68')
    d=json.loads(source.fetch(source.BASE+'api/records/6838169',2*1024**2));f=next(x for x in d['files'] if x['key']=='Melanoma_code_data.zip')
    assert f['size']==668227028 and f['checksum']=='md5:da74b7228763b76dcd2d9830c1a10cd9'
    rf=transport.RangeFile(source,f['links']['self'],f['size'],checkpoint=lambda v:save('private_transport_progress.json',v))
    with zipfile.ZipFile(rf) as z:
        a=z.read('data/input/meta/ROI_meta_afterQC.csv');b=z.read('data/input/hc_result/hc_result_sample_anno_invasive.csv')
        assert sha(a)=='d0fa2438b7b91157a7848654fb1f23e40c9dbcca8020f7f5ff0566ce7acede21'
        assert sha(b)=='59c9f04e4e12cc7484f1cc262b96b77c7a5eae4efff7d2510fbc23921b6f55dc'
        registry=list(csv.DictReader(io.StringIO(a.decode('utf-8-sig'))))
        # The wide table has coded headers; export neither headers nor unused values.
        bridge=[dict(ROI=r['ROI'],patient=r['patient']) for r in csv.DictReader(io.StringIO(b.decode('utf-8-sig')))]
        mapping,keys=explicit_margin_map(registry,bridge)
        buffer=io.BytesIO();total=0
        info=z.getinfo('data/input/sc_data.csv');assert info.file_size==prior['uncompressed_bytes']
        with z.open(info) as stream:
            while True:
                chunk=stream.read(4*1024**2)
                if not chunk:break
                total+=len(chunk);assert total<=info.file_size;buffer.write(chunk)
        assert total==info.file_size;blob=buffer.getvalue();assert sha(blob)==prior['complete_member_SHA256']
        frame=pd.read_csv(io.BytesIO(blob),usecols=['sample_id','celltype20']+MARKERS,
            dtype={**{m:np.float64 for m in MARKERS},'sample_id':str,'celltype20':str},keep_default_na=False,float_precision='round_trip')
        assert len(frame)==prior['cells']
        assert frame['celltype20'].value_counts().to_dict()==prior['author_label_counts']['celltype20']
        result,cells,means=summarize(frame,mapping,keys);result.update(scalar_replay(blob,mapping,cells,means))
    final=json.loads(source.fetch(source.BASE+'api/records/6838169',2*1024**2));other=next(x for x in final['files'] if x['key']==f['key'])
    assert all(f.get(k)==other.get(k) for k in ['id','size','checksum','links'])
    assert int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*1024<=2*1024**3
    rf.close();result.update(status='COMPLETE_OUTCOME_BLIND_PATIENT_TYPE_MARKER_COVERAGE',completed_UTC=datetime.now(timezone.utc).isoformat(),
        producer_SHA256=sha(Path(__file__).read_bytes()),accepted_measurement_SHA256=sha(accepted.read_bytes()),complete_cell_member_SHA256=sha(blob),
        complete_source_CRC_SHA_and_published_metadata_checked=True,no_individual_values_or_arrays_written=True)
    save('private_complete_patient_type_marker_coverage.json',result)
    save('private_terminal.json',dict(status=result['status'],completed_UTC=datetime.now(timezone.utc).isoformat(),producer_SHA256=sha(Path(__file__).read_bytes()),
        files={p.name:dict(SHA256=sha(p.read_bytes()),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status=result['status'],terminal_SHA256=sha((out/'private_terminal.json').read_bytes()))))


if __name__=='__main__':
    try:main()
    except Exception as error:
        out=Path(os.environ['PANCANCER_RAM_OUTPUT']);out.mkdir(mode=0o700,parents=True,exist_ok=True)
        value=dict(status='FAILED_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],no_record_values_in_error=True)
        (out/'private_failed_terminal.json').write_text(json.dumps(value,indent=2)+'\n')
        print(json.dumps(dict(status=value['status'],error_class=value['error_class'])))
        raise SystemExit(1)
