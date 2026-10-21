"""Same-source local patient/type/marker coverage with exact blocked neighbors.

Requires a completed, SHA-checked cell-measurement qualification first. It is
not an association test, CM projection or new cell-state annotation. Multiple
ROIs are pooled within the explicit source patient. No IDs/rows/arrays saved.
"""
from collections import Counter,defaultdict
from datetime import datetime,timezone
from itertools import combinations
from pathlib import Path
import csv,hashlib,importlib.util,io,json,os,re,resource,sys,time,traceback,zipfile
import numpy as np
import pandas as pd

MARKERS=['CD3','CD4','CD8a','FoxP3','Ki-67','LAG3','PD1','TIM3','OX40','CD127']
PARENT_CELL_MINIMUM=20


def sha(b):return hashlib.sha256(b).hexdigest()


def load(name,path,expected):
    assert sha(path.read_bytes())==expected
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


_local_neighbors=load('exact_local_neighbors',Path(__file__).with_name('melanoma_IMC_numpy_neighbors_v1.py'),'acbde69dae387ee79d6d8397899681bedf00553dd3da08878abd16b392858471')
cKDTree=_local_neighbors.ExactEuclideanNeighbors
peak_rss_bytes=_local_neighbors.peak_rss_bytes


def directory_signature(entries):
    return sha(json.dumps([(e.filename,e.file_size,e.compress_size,e.CRC) for e in entries],separators=(',',':')).encode())


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
    safe=re.compile(r'^[\w][\w +():./^=\-]{0,80}$')
    labels={x:x if safe.fullmatch(x) and x not in set(patients)|registry_keys else 'unexported_author_label_'+str(i+1) for i,x in enumerate(types)}
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


def scalar_replay(blob,roi_to_patient,cells,means,raw_to_canonical=None):
    counter=Counter();sums=defaultdict(lambda:np.zeros(len(MARKERS),dtype=float));n=0
    for row in csv.DictReader(io.TextIOWrapper(io.BytesIO(blob),encoding='utf-8-sig',newline='')):
        roi=raw_to_canonical[row['sample_id']] if raw_to_canonical is not None else row['sample_id']
        if roi not in roi_to_patient:continue
        key=(roi_to_patient[roi],row['celltype20'])
        value=np.array([float(row[m]) for m in MARKERS]);assert np.isfinite(value).all()
        counter[key]+=1;sums[key]+=value;n+=1
    for patient in cells.index:
        for typ in cells.columns:assert int(cells.loc[patient,typ])==counter[(patient,typ)]
    for key,count in counter.items():
        assert np.allclose(sums[key]/count,means.loc[key,MARKERS].to_numpy(dtype=float),rtol=1e-10,atol=1e-10)
    assert n==int(cells.to_numpy().sum())
    return dict(all_CSV_dictionary_patient_type_counts_and_marker_means_replayed=True,
        independent_reduction_not_independent_input_bytes=True)


def author_ROI_infix(text):
    assert "if (stringr::str_detect(each, 'F'))" in text
    assert "splitID = str_split(each,'_')[[1]]" in text
    line=next(s for s in text.splitlines() if 'a = paste0(splitID[1]' in s)
    match=re.fullmatch(r'''\s*a\s*=\s*paste0\(\s*splitID\[1\]\s*,\s*(["'])([^"'\\\r\n]{0,64})\1\s*,\s*splitID\[2\]\s*\)\s*''',line)
    assert match is not None
    return match.group(2)


def canonical_ROIs(raw_keys,registry_keys,infix):
    mapping={}
    for raw in raw_keys:
        if 'F' in raw:canonical=raw
        else:
            parts=raw.split('_');assert len(parts)>=2
            canonical=parts[0]+infix+parts[1]
        mapping[raw]=canonical
    assert len(set(mapping.values()))==len(mapping)
    assert set(mapping.values())==registry_keys
    return mapping


def coordinate_array(series):
    number=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
    clean=series.str.strip()
    assert (clean.str.startswith('(')==clean.str.endswith(')')).all()
    parts=clean.str.extract(r'^\(?\s*('+number+r')\s*,\s*('+number+r')\s*\)?$')
    result=parts.to_numpy(dtype=float);assert np.isfinite(result).all()
    return result


def coordinate_summary(frame,coords):
    assert coords.shape==(len(frame),2) and np.isfinite(coords).all()
    duplicates=0
    for indices in frame.groupby('sample_id',sort=False).indices.values():
        duplicates+=len(indices)-len(np.unique(coords[indices],axis=0))
    return dict(native_coordinate_pairs=len(coords),finite_all=True,within_ROI_duplicate_coordinate_excess=int(duplicates),
        different_ROI_coordinate_frames_never_stitched=True,physical_um_scale_not_accepted=True)


def native_knn_edges(coords,k=10):
    """Union of kNN with all exact boundary ties; uniform-scale invariant."""
    n=len(coords)
    if n<2:return np.empty((0,2),dtype=np.int64)
    _,multiplicity=np.unique(coords,axis=0,return_counts=True)
    assert multiplicity.max()<=50,'degenerate_coincident_geometry_requires_review'
    k=min(k,n-1);tree=cKDTree(coords)
    distances,_=tree.query(coords,k=k+1,workers=1)
    radius=np.nextafter(distances[:,-1],np.inf)
    neighborhoods=tree.query_ball_point(coords,radius,workers=1)
    edges=[]
    for i,neighbors in enumerate(neighborhoods):
        edges.extend((min(i,j),max(i,j)) for j in neighbors if i!=j)
    return np.unique(np.asarray(edges,dtype=np.int64).reshape(-1,2),axis=0)


def graph_coverage(frame,coords,mapping,cells,k=10):
    types=list(cells.columns);patients=list(cells.index);lookup={v:i for i,v in enumerate(types)}
    patient_index={v:i for i,v in enumerate(patients)};counts=np.zeros((len(patients),len(types),len(types)),dtype=np.int64)
    total_edges=0
    for roi,indices in frame.groupby('sample_id',sort=False).indices.items():
        if roi not in mapping:continue
        edges=native_knn_edges(coords[indices],k);total_edges+=len(edges)
        codes=np.array([lookup[x] for x in frame.iloc[indices]['celltype20']],dtype=np.int64)
        if len(edges):
            a=np.minimum(codes[edges[:,0]],codes[edges[:,1]])
            b=np.maximum(codes[edges[:,0]],codes[edges[:,1]])
            np.add.at(counts[patient_index[mapping[roi]]],(a,b),1)
    safe=re.compile(r'^[\w][\w +():./^=\-]{0,80}$');sensitive=set(patients)|set(frame['sample_id'])
    labels={t:t if safe.fullmatch(t) and t not in sensitive else 'unexported_author_label_'+str(i+1) for i,t in enumerate(types)}
    rows=[]
    for i,a in enumerate(types):
        for j in range(i+1,len(types)):
            b=types[j];eligible=(cells[a]>=20)&(cells[b]>=20);edges=counts[:,i,j]
            rows.append(dict(author_type_A=labels[a],author_type_B=labels[b],
                patients_any_native_neighbor_edge=int((edges>0).sum()),
                patients_both_types_ge20_and_atleast20_edges=int((eligible.to_numpy()&(edges>=20)).sum())))
    return dict(k=k,ties_included=True,source_margin_undirected_edges=int(total_edges),pair_opportunities=rows,
        within_ROI_only=True,no_micron_distance_claim=True,no_state_colocalization_or_outcome_test=True)


def source_summary_only(value,sensitive):
    """Prevent accidental return of records/arrays from the next own module."""
    if isinstance(value,dict):
        forbidden={'patient','patient_id','sample_id','cell_id','barcode','donor','matrix','weights','patient_scores'}
        assert not (set(value)&forbidden)
        for k,v in value.items():
            assert str(k) not in sensitive;source_summary_only(v,sensitive)
    elif isinstance(value,list):
        assert len(value)<=1000
        if value and all(isinstance(x,(float,int)) for x in value):assert len(value)<=5
        for v in value:source_summary_only(v,sensitive)
    elif isinstance(value,str):assert value not in sensitive
    else:assert value is None or isinstance(value,(int,float,bool))


def hold_for_frozen_step(out,frame,coords,mapping,aliases,save,author_nodes):
    root=Path(__file__).parent;started=time.monotonic();done=set()
    bundle=dict(frame=frame,native_coordinates=coords,ROI_to_source_patient=mapping,raw_ROI_aliases=aliases,
        native_knn_edges=native_knn_edges,author_node_table=author_nodes,source_code_SHA256=sha(Path(__file__).read_bytes()))
    def session(state,**extra):
        save('private_RAM_session_status.json',dict(status=state,observed_UTC=datetime.now(timezone.utc).isoformat(),
            no_individual_arrays_or_rows_written=True,**extra))
    session('READY_WITH_RAM_FOR_SEPARATELY_FROZEN_OWN_ANALYSIS',maximum_hold_seconds=75*60)
    sensitive=set(mapping)|set(mapping.values())|set(aliases)
    while time.monotonic()-started<75*60:
        requests=sorted(root.glob('frozen_step_request_*.json'))
        for request in requests:
            if request.name in done:continue
            try:
                protocol=json.loads(request.read_text());assert set(protocol)=={'step','code_file','code_SHA256','parameters'}
                step=protocol['step'];assert re.fullmatch('[A-Za-z0-9_]{1,60}',step)
                assert not (out/('private_frozen_step_'+step+'_source_summary.json')).exists()
                filename=protocol['code_file'];assert Path(filename).name==filename and filename.endswith('.py')
                module_path=root/filename;assert module_path.resolve().parent==root.resolve() and not module_path.is_symlink()
                module=load('frozen_'+step,module_path,protocol['code_SHA256'])
                session('RUNNING_SEPARATELY_FROZEN_STEP',step=step,protocol_SHA256=sha(request.read_bytes()),code_SHA256=protocol['code_SHA256'])
                independent_bundle=dict(bundle,frame=frame.copy(deep=True),native_coordinates=coords.copy(),
                    ROI_to_source_patient=dict(mapping),raw_ROI_aliases=dict(aliases))
                result=module.analyze(independent_bundle,protocol['parameters'])
                source_summary_only(result,sensitive)
                assert peak_rss_bytes()<=2*1024**3
                save('private_frozen_step_'+step+'_source_summary.json',result)
            except Exception as error:
                save('private_frozen_step_'+request.stem+'_not_accepted.json',dict(status='FROZEN_STEP_NOT_ACCEPTED',
                    error_class=type(error).__name__,frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],
                    completed_coverage_phase_remains_separate=True,no_record_values_in_error=True))
            finally:
                if 'independent_bundle' in locals():del independent_bundle
            done.add(request.name)
            session('READY_WITH_RAM_FOR_SEPARATELY_FROZEN_OWN_ANALYSIS',completed_steps=sorted(done))
        if len(done)>=3:break
        time.sleep(5)
    session('RAM_SESSION_CLOSED',completed_steps=sorted(done))


def main():
    # A partial file, progress record or process lock cannot satisfy this gate.
    accepted=Path(os.environ['PANCANCER_IMC_ACCEPTED_MEASUREMENT'])
    terminal=accepted.parent/'private_terminal.json';t=json.loads(terminal.read_text());prior=json.loads(accepted.read_text())
    assert t['status']==prior['status']=='COMPLETE_SELECTED_PROCESSED_CELL_MEMBER_MEASUREMENT_QUALIFICATION'
    assert t['producer_SHA256']=='eb142bcdeb172f068f9d182a6ba89866d3d07e852c75f0b976ec4c5b43cd436d'
    assert t['files'][accepted.name]['SHA256']==sha(accepted.read_bytes())
    assert prior['complete_member_CRC_checked'] and prior['complete_independent_CSV_numeric_ROI_and_label_counter_replay']
    assert prior['published_metadata_and_directory_replayed']
    assert prior['source_ROI_key_count']==prior['cell_table_source_key_count']
    out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert not out.exists();out.mkdir(mode=0o700)
    def save(n,v):(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
    save('plan_before_patient_coverage.json',dict(frozen_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(Path(__file__).read_bytes()),
        accepted_measurement_SHA256=sha(accepted.read_bytes()),scope='All author types/pairs/ten marker means, source-declared ROI aliases, native XY and fixed k10 graph coverage. No marker-spatial or outcome association.',
        markers=MARKERS,parent_cell_minimum=PARENT_CELL_MINIMUM,source_patient_unit='Explicit patient field and invasive-margin literal ROI map only.',
        no_ROI_pseudoreplication_or_prefix_inference=True,
        author_defined_ROI_format_derived_from_fixed_source_only=True,native_kNN_k=10,no_physical_um_scale_assumed=True,
        optional_bounded_RAM_hold_for_separately_frozen_own_code=True,
        cached_source_contract_bound_by_prior_complete_member_SHA=True,fresh_metadata_API_not_selected=True,
        local_exact_neighbors_helper_SHA256='acbde69dae387ee79d6d8397899681bedf00553dd3da08878abd16b392858471',local_range_max_chunk_bytes=262144,coalesced_RAM_cache_blocks=8,fresh_wire_directory_replay=True,no_protein_positive_threshold_or_clinical_tests=True,
        no_individual_keys_rows_expression_or_learned_arrays_written=True))
    helper=Path(os.environ['PANCANCER_IMC_TRANSPORT_HELPER'])
    source=load('published_archive_transport',helper,'9fb200d98e321584fde8bca6cb92f72151ac33ec00d64e04d93e7e2b33e870e2')
    transport=load('bounded_archive_transport',Path(os.environ['PANCANCER_IMC_BOUNDED_TRANSPORT']),'b099f8ffe429a856a70a5b9006fe871e944dc539943fb28bbdda354dd9e222c8')
    f=dict(key='Melanoma_code_data.zip',size=668227028,checksum='md5:da74b7228763b76dcd2d9830c1a10cd9',links=dict(self='https://zenodo.org/api/records/6838169/files/Melanoma_code_data.zip/content'))
    assert f['size']==668227028 and f['checksum']=='md5:da74b7228763b76dcd2d9830c1a10cd9'
    rf=transport.RangeFile(source,f['links']['self'],f['size'],checkpoint=lambda v:save('private_transport_progress.json',v))
    rf.max_chunk=256*1024
    with zipfile.ZipFile(rf) as z:
        assert directory_signature(z.infolist())=='26339da50b39ffdf7beb641a67f50523c867f2d3a86b24f42175f3a16dfe8992'
        a=z.read('data/input/meta/ROI_meta_afterQC.csv');b=z.read('data/input/hc_result/hc_result_sample_anno_invasive.csv')
        assert sha(a)=='d0fa2438b7b91157a7848654fb1f23e40c9dbcca8020f7f5ff0566ce7acede21'
        assert sha(b)=='59c9f04e4e12cc7484f1cc262b96b77c7a5eae4efff7d2510fbc23921b6f55dc'
        registry=list(csv.DictReader(io.StringIO(a.decode('utf-8-sig'))))
        # The wide table has coded headers; export neither headers nor unused values.
        bridge=[dict(ROI=r['ROI'],patient=r['patient']) for r in csv.DictReader(io.StringIO(b.decode('utf-8-sig')))]
        mapping,keys=explicit_margin_map(registry,bridge)
        author_code=z.read('Rscript/functions.R');assert sha(author_code)=='bfb6f99c80947643655aa5d4d8ef310e5f0a1c24aac86326482cb32c123565ee'
        infix=author_ROI_infix(author_code.decode('utf-8-sig'))
        # Read the small complete author node table FIRST and retain RAM.
        spatial_code=z.read('Rscript/Spatial_analysis.Rmd')
        assert sha(spatial_code)=='d8faddad268f2900a83f8cfa3f8352dd697128d408a528181fb997c5773dba56'
        line=next(v for v in spatial_code.decode('utf-8-sig').splitlines() if re.match(r'\s*sc_data\s*<-\s*read.csv',v))
        literal=re.search(r"read.csv\(\s*([\"\'])([^\"\'\r\n]+)\1",line).group(2)
        while literal.startswith('../') or literal.startswith('./'):literal=literal.split('/',1)[1]
        assert z.getinfo(literal).file_size<=40*1024**2
        node_blob=z.read(literal)
        assert sha(node_blob)=='1b3f8e10a296794d97dc467681dde9f22cfa6935a7137e4c4230b397dc554f6b'
        nodes=pd.read_csv(io.BytesIO(node_blob),dtype=str,keep_default_na=False).rename(columns={'Unnamed: 0':''})
        scalar_nodes=list(csv.DictReader(io.StringIO(node_blob.decode('utf-8-sig'))))
        assert len(scalar_nodes)==len(nodes) and list(scalar_nodes[0])==list(nodes)
        assert all([v[col] for v in scalar_nodes]==nodes[col].tolist() for col in nodes)
        del node_blob,scalar_nodes
        save('private_node_preloaded_qualification.json',dict(complete_node_member_CRC_checked=True,
            complete_node_member_SHA256='1b3f8e10a296794d97dc467681dde9f22cfa6935a7137e4c4230b397dc554f6b',
            independent_CSV_Pandas_node_replay=True,no_identity_accepted_before_main=True,all_individual_arrays_RAM_only=True))
        buffer=io.BytesIO();total=0
        info=z.getinfo('data/input/sc_data.csv');assert info.file_size==prior['uncompressed_bytes']
        with z.open(info) as stream:
            while True:
                chunk=stream.read(4*1024**2)
                if not chunk:break
                total+=len(chunk);assert total<=info.file_size;buffer.write(chunk)
        assert total==info.file_size;blob=buffer.getvalue();assert sha(blob)==prior['complete_member_SHA256']
        frame=pd.read_csv(io.BytesIO(blob),usecols=['sample_id','celltype20','celltype4','Position','X','Unnamed: 0','Area']+MARKERS,
            dtype={**{m:np.float64 for m in MARKERS},'sample_id':str,'celltype20':str,'celltype4':str,'Position':str,'X':str,'Unnamed: 0':str,'Area':str},keep_default_na=False,float_precision='round_trip')
        assert len(frame)==prior['cells']
        if prior['author_label_counts']['celltype20'] is not None:
            assert frame['celltype20'].value_counts().to_dict()==prior['author_label_counts']['celltype20']
        assert frame.groupby('celltype20')['celltype4'].nunique().max()==1
        raw_to_canonical=canonical_ROIs(set(frame['sample_id']),keys,infix)
        frame['sample_id']=frame['sample_id'].map(raw_to_canonical)
        counter_bridge=load('explicit_source_counter_bridge',Path(__file__).with_name('melanoma_IMC_explicit_counter_bridge_v1.py'),'9eaee9aae48eb448a7164464d8e0f2a63e29afff9872aa706342e8da687d67da')
        _,_,counter_proof=counter_bridge.bind(frame,nodes,raw_to_canonical)
        save('private_author_node_counter_qualification.json',dict(counter_proof,
            status='COMPLETE_EXPLICIT_COUNTER_BASE_AND_LABEL_BIJECTION',
            complete_node_member_CRC_checked=True,independent_CSV_Pandas_node_replay=True,
            complete_node_member_SHA256='1b3f8e10a296794d97dc467681dde9f22cfa6935a7137e4c4230b397dc554f6b'))
        result,cells,means=summarize(frame,mapping,keys);result.update(scalar_replay(blob,mapping,cells,means,raw_to_canonical))
        save('private_basic_patient_type_marker_coverage.json',dict(result,stage='COMPLETE_BASIC_COVERAGE_BEFORE_GEOMETRY',
            producer_SHA256=sha(Path(__file__).read_bytes()),complete_cell_member_SHA256=sha(blob)))
        coords=coordinate_array(frame['Position'])
        result['native_coordinate_qualification']=coordinate_summary(frame,coords)
        try:result['native_k10_graph_coverage']=graph_coverage(frame,coords,mapping,cells)
        except Exception as error:
            result['native_k10_graph_coverage']=dict(status='NATIVE_GRAPH_NOT_ACCEPTED',error_class=type(error).__name__,
                basic_patient_measurement_kept_separate=True,no_marker_or_outcome_test_computed=True)
        result['ROI_alias_qualification']=dict(author_source_SHA256=sha(author_code),infix_SHA256=sha(infix.encode()),
            raw_ROIs=len(raw_to_canonical),changed_ROIs=sum(k!=v for k,v in raw_to_canonical.items()),
            all_canonical_ROIs_exactly_equal_registry=True,no_many_to_one_merge=True)
        source_blob_SHA256=sha(blob);del blob,buffer
    rf.clear_cache()
    with zipfile.ZipFile(rf) as replay:
        assert directory_signature(replay.infolist())=='26339da50b39ffdf7beb641a67f50523c867f2d3a86b24f42175f3a16dfe8992'
    assert peak_rss_bytes()<=2*1024**3
    rf.close();result.update(status='COMPLETE_OUTCOME_BLIND_PATIENT_TYPE_MARKER_COVERAGE',completed_UTC=datetime.now(timezone.utc).isoformat(),
        producer_SHA256=sha(Path(__file__).read_bytes()),accepted_measurement_SHA256=sha(accepted.read_bytes()),complete_cell_member_SHA256=source_blob_SHA256,
        complete_source_CRC_SHA_and_cached_source_contract_checked=True,fresh_metadata_API_replayed=False,
        neighbor_backend='exact_blocked_full_Euclidean_NumPy_not_approximate',peak_RSS_bytes=peak_rss_bytes(),no_individual_values_or_arrays_written=True)
    save('private_complete_patient_type_marker_coverage.json',result)
    save('private_terminal.json',dict(status=result['status'],completed_UTC=datetime.now(timezone.utc).isoformat(),producer_SHA256=sha(Path(__file__).read_bytes()),
        files={p.name:dict(SHA256=sha(p.read_bytes()),bytes=p.stat().st_size) for p in out.iterdir()}))
    print(json.dumps(dict(status=result['status'],terminal_SHA256=sha((out/'private_terminal.json').read_bytes()))),flush=True)
    if os.environ.get('PANCANCER_RETAIN_RAM_FOR_FROZEN_STEP')=='1':
        hold_for_frozen_step(out,frame,coords,mapping,raw_to_canonical,save,nodes)


if __name__=='__main__':
    try:main()
    except Exception as error:
        out=Path(os.environ['PANCANCER_RAM_OUTPUT']);out.mkdir(mode=0o700,parents=True,exist_ok=True)
        value=dict(status='FAILED_NOT_ACCEPTED',observed_UTC=datetime.now(timezone.utc).isoformat(),error_class=type(error).__name__,
            frames=[dict(file=Path(x.filename).name,line=x.lineno,function=x.name) for x in traceback.extract_tb(error.__traceback__)],no_record_values_in_error=True)
        (out/'private_failed_terminal.json').write_text(json.dumps(value,indent=2)+'\n')
        print(json.dumps(dict(status=value['status'],error_class=value['error_class'])))
        raise SystemExit(1)
