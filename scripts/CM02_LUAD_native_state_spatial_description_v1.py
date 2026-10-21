"""RAM-only fixed RNA-pattern and native-neighborhood descriptive measurements.

No clinical/identity rows, clonotypes, P values, learned weights or saved arrays.
Author cell bins and deposited positive integers are not verified original UMIs.
"""
from collections import Counter
from datetime import datetime, timezone
import hashlib
import io
import json
import resource
import numpy as np

STATES_R = ['R_TNFRSF9', 'R_LAG3', 'R_cycling']
STATES_E = ['E_dysfunction', 'E_cycling', 'E_memory']


def text(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def summarize(values):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if not len(values):
        return dict(n=0, median=None, IQR=None, positive=0, negative=0, leave_one_image_out_median_range=None)
    loo = [float(np.median(np.delete(values, i))) for i in range(len(values))] if len(values)>1 else []
    return dict(n=len(values), median=float(np.median(values)),
                IQR=np.quantile(values, [.25, .75]).tolist(),
                positive=int((values>0).sum()), negative=int((values<0).sum()),
                leave_one_image_out_median_range=[min(loo),max(loo)] if loo else None)


def phenotype(counts, evidence):
    good = np.ones(len(counts), bool)
    f = evidence.vector_evidence(counts, good, 'dual_regulatory_evidence')
    np.testing.assert_array_equal(f, evidence.scalar_evidence(counts, good, 'dual_regulatory_evidence'))
    other = f[:,0] & ~f[:,1:4].any(axis=1)
    x = {g:counts[:,i]>0 for i,g in enumerate(evidence.GENES)}
    cycle = x['MKI67'] & (x['TOP2A'] | x['UBE2C'])
    return f, other, np.column_stack([other & x['TNFRSF9'], other & x['LAG3'], other & cycle])


def gene_capture(h, evidence):
    var = h['var']; index = text(var.attrs.get('_index','_index'))
    names = [text(x).strip().upper() for x in var[index][:]]
    counts = Counter(names)
    assert all(counts[g]==1 for g in evidence.GENES), 'FIXED_GENES_NOT_UNIQUELY_MEASURED'
    n, nf = tuple(h['X'].attrs['shape'])
    assert nf==len(names)
    lookup = {g:i for i,g in enumerate(names)}
    mapping = np.full(nf, -1, int)
    for i,g in enumerate(evidence.GENES):mapping[lookup[g]]=i
    matrix = h['X']; values = np.zeros((n,len(evidence.GENES)),np.int32)
    total = np.zeros(n); detected = np.zeros(n,np.int32)
    assert text(matrix.attrs['encoding-type'])=='csr_matrix'
    ptr = np.asarray(matrix['indptr'],np.int64)
    assert len(ptr)==n+1 and ptr[0]==0 and (np.diff(ptr)>=0).all()
    assert ptr[-1]==len(matrix['data'])==len(matrix['indices'])
    for start in range(0,n,4096):
        stop=min(start+4096,n);lo,hi=int(ptr[start]),int(ptr[stop]);p=ptr[start:stop+1]-lo
        data=np.asarray(matrix['data'][lo:hi]);indices=np.asarray(matrix['indices'][lo:hi],int)
        assert np.isfinite(data).all() and (data>0).all() and np.array_equal(data,np.floor(data))
        assert len(data)==0 or data.max()<2**31
        assert (indices>=0).all() and (indices<nf).all()
        row=np.repeat(np.arange(stop-start),np.diff(p))
        if len(indices)>1:assert ((indices[1:]>indices[:-1]) | (row[1:]!=row[:-1])).all()
        totals=np.r_[0.,np.cumsum(data,dtype=float)]
        total[start:stop]=totals[p[1:]]-totals[p[:-1]];detected[start:stop]=np.diff(p)
        take=mapping[indices]>=0;values[start+row[take],mapping[indices[take]]]=data[take]
        # Independent direct indexed extraction checks every requested column in each block.
        for j in range(len(evidence.GENES)):
            hit=indices==lookup[evidence.GENES[j]]
            expected=np.bincount(row[hit],weights=data[hit],minlength=stop-start)
            np.testing.assert_array_equal(values[start:stop,j],expected)
    coords=np.asarray(h['obsm']['spatial'])
    assert len(coords)==n and coords.ndim==2 and coords.shape[1]>=2 and np.isfinite(coords).all()
    return values,coords[:,:2].astype(float),total,detected


def graph(xy,k):
    from scipy.spatial import cKDTree
    assert len(xy)>k
    tree=cKDTree(xy);distance,_=tree.query(xy,k=k+1,workers=1)
    median=float(np.median(distance[:,1]));assert median>0
    radii=np.minimum(distance[:,k],5*median)*(1+1e-12)
    lengths=tree.query_ball_point(xy,radii,return_length=True,workers=1)
    assert lengths.max()<=4096 and lengths.sum()<=128*len(xy), 'NATIVE_GRAPH_RESOURCE_CAP'
    sources=[];targets=[]
    for lo in range(0,len(xy),10000):
        near=tree.query_ball_point(xy[lo:lo+10000],radii[lo:lo+10000],return_sorted=True,workers=1)
        for j,items in enumerate(near):
            anchor=lo+j;items=np.asarray(items,int);items=items[items!=anchor]
            sources.append(np.full(len(items),anchor,int));targets.append(items)
    return np.concatenate(sources),np.concatenate(targets),median


def edge_effect(parent, state, receiver, receiver_state, source, target):
    group=parent & state;comparison=parent & ~state
    valid_parent=np.bincount(source,weights=receiver[target],minlength=len(parent))>0
    gate=dict(parent=int(parent.sum()),receiver_parent=int(receiver.sum()),
              positive_anchors=int(group.sum()),negative_anchors=int(comparison.sum()),
              receiver_positive=int(receiver_state.sum()),receiver_negative=int((receiver & ~receiver_state).sum()),
              positive_anchors_with_receiver_neighbors=int((group & valid_parent).sum()),
              negative_anchors_with_receiver_neighbors=int((comparison & valid_parent).sum()))
    if min(gate['parent'],gate['receiver_parent'])<20:return None,'PARENT_BELOW20'
    if min(gate[k] for k in ['positive_anchors','negative_anchors','receiver_positive','receiver_negative',
                           'positive_anchors_with_receiver_neighbors','negative_anchors_with_receiver_neighbors'])<10:
        return None,'STATE_OR_CONTRAST_CAPTURE_BELOW10'
    categories=[group,comparison];targets=[receiver_state,receiver & ~receiver_state]
    observed=np.asarray([[np.sum(a[source] & b[target]) for b in targets] for a in categories],np.int64)
    # A second implementation collapses receiver mark counts at each anchor first.
    by_anchor=[np.bincount(source,weights=b[target],minlength=len(parent)).astype(np.int64) for b in targets]
    independent=np.asarray([[v[a].sum() for v in by_anchor] for a in categories],np.int64)
    np.testing.assert_array_equal(observed,independent)
    a,b,c,d=observed.ravel()
    return float(np.log((a+.5)*(d+.5)/((b+.5)*(c+.5)))),'DESCRIPTIVE_ELIGIBLE_NO_CELL_LEVEL_P'


def describe(values,xy,evidence):
    flags,other,other_states=phenotype(values,evidence)
    _,inverse,frequencies=np.unique(xy,axis=0,return_inverse=True,return_counts=True)
    unique=frequencies[inverse]==1
    xy=xy[unique];flags=flags[unique];other=other[unique];other_states=other_states[unique]
    effects=[]
    for k in [5,10,20]:
        if len(xy)<=k:
            for r in STATES_R:
                for e in STATES_E:effects.append(dict(k=k,R=r,E=e,R_value=None,O_value=None,R_status='GEOMETRY_TOO_SMALL',O_status='GEOMETRY_TOO_SMALL'))
            continue
        source,target,median=graph(xy,k)
        for i,r in enumerate(STATES_R):
            for j,e in enumerate(STATES_E):
                v,why=edge_effect(flags[:,1],flags[:,5+i],flags[:,2],flags[:,8+j],source,target)
                o,owhy=edge_effect(other,other_states[:,i],flags[:,2],flags[:,8+j],source,target)
                effects.append(dict(k=k,R=r,E=e,R_value=v,O_value=o,R_status=why,O_status=owhy))
    return effects,dict(author_rows=len(values),unique_coordinate_rows=int(unique.sum()),
        duplicate_coordinate_rows_excluded=int((~unique).sum()),
        total_pattern_counts={name:int(flags[:,i].sum()) for i,name in enumerate(evidence.METRICS)})


def synthetic_native_checks(evidence):
    import h5py
    from scipy.sparse import csr_matrix
    rng=np.random.default_rng(72833);values=rng.integers(0,3,size=(64,len(evidence.GENES)),dtype=np.int32)
    xy=np.column_stack(np.unravel_index(np.arange(64),(8,8))).astype(float)
    matrix=csr_matrix(values[:,::-1]);buffer=io.BytesIO()
    with h5py.File(buffer,'w') as h:
        var=h.create_group('var');var.attrs['_index']='_index'
        var.create_dataset('_index',data=np.array(evidence.GENES[::-1],dtype=h5py.string_dtype()))
        x=h.create_group('X');x.attrs['shape']=matrix.shape;x.attrs['encoding-type']='csr_matrix'
        for name in ['data','indices','indptr']:x.create_dataset(name,data=getattr(matrix,name))
        h.create_group('obsm').create_dataset('spatial',data=xy)
        got,coords,total,detected=gene_capture(h,evidence)
        np.testing.assert_array_equal(got,values);np.testing.assert_array_equal(coords,xy)
        np.testing.assert_array_equal(total,values.sum(axis=1));np.testing.assert_array_equal(detected,(values>0).sum(axis=1))
    for k in [5,10,20]:
        source,target,median=graph(xy,k)
        distances=np.sqrt(((xy[:,None]-xy[None,:])**2).sum(axis=2))
        radius=np.minimum(np.sort(distances,axis=1)[:,k],5*np.median(np.sort(distances,axis=1)[:,1]))*(1+1e-12)
        expected=(distances<=radius[:,None]) & ~np.eye(len(xy),dtype=bool)
        observed=np.zeros_like(expected);observed[source,target]=True
        np.testing.assert_array_equal(observed,expected)
    effects,coverage=describe(values,xy,evidence);assert len(effects)==27 and coverage['author_rows']==64
    buffer.close()


def remote_main(scope,evidence):
    import h5py
    import sys
    minimum=int(next(s.split()[1] for s in open('/proc/meminfo') if s.startswith('MemAvailable:')))*1024
    assert minimum>=9*1024**3,'SOURCE_MEMORY_HEADROOM_NOT_ACCEPTED'
    synthetic_native_checks(evidence)
    print(json.dumps(dict(kind='ready',status='SYNTHETIC_H5AD_ALL_TARGETS_AND_TIE_SAFE_GEOMETRY_ACCEPTED',
        observed_UTC=datetime.now(timezone.utc).isoformat())),flush=True)
    all_effects=[];coverage=[];statuses=[]
    for ordinal,expected in enumerate(scope['sources']):
        header=json.loads(sys.stdin.buffer.readline());assert header==expected
        size=int(expected['bytes']);buffer=io.BytesIO();digest=hashlib.sha256()
        while buffer.tell()<size:
            chunk=sys.stdin.buffer.read(min(4*1024**2,size-buffer.tell()))
            assert chunk,'SOURCE_PIPE_ENDED_BEFORE_DECLARED_SIZE'
            buffer.write(chunk);digest.update(chunk)
        assert digest.hexdigest()==expected['SHA256'],'SOURCE_FILE_NOT_IDENTICAL_TO_PREVIOUS_FULL_QUALIFICATION'
        buffer.seek(0)
        try:
            with h5py.File(buffer,'r') as h:values,xy,total,detected=gene_capture(h,evidence)
            buffer.close()
            effects,counts=describe(values,xy,evidence)
            all_effects.append(effects);coverage.append(counts);statuses.append('COMPLETE_ALL_FIXED_DESCRIPTIVE_SLOTS')
        except (AssertionError,ValueError,KeyError) as error:
            buffer.close();statuses.append('UNQUALIFIED_'+type(error).__name__)
        assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024<=6*1024**3,'SOURCE_RSS_BUDGET_EXCEEDED'
        print(json.dumps(dict(kind='progress',finished_images=ordinal+1,total_images=len(scope['sources']),
            status=statuses[-1],observed_UTC=datetime.now(timezone.utc).isoformat())),flush=True)
    summary=[]
    for k in [5,10,20]:
        for r in STATES_R:
            for e in STATES_E:
                records=[next(v for v in image if (v['k'],v['R'],v['E'])==(k,r,e)) for image in all_effects]
                rv=[x['R_value'] for x in records if x['R_value'] is not None]
                ov=[x['O_value'] for x in records if x['O_value'] is not None]
                delta=[x['R_value']-x['O_value'] for x in records if x['R_value'] is not None and x['O_value'] is not None]
                summary.append(dict(k=k,R_state=r,E_state=e,R_summary=summarize(rv),other_T=summarize(ov),matched_R_minus_other=summarize(delta),
                    R_status_counts=dict(Counter(x['R_status'] for x in records)),other_status_counts=dict(Counter(x['O_status'] for x in records))))
    result=dict(status='COMPLETE_ALL11_CURRENT_NATIVE_SCOPE_DESCRIPTIVE_NO_INFERENCE',
        source='GSE328481',source_status_counts=dict(Counter(statuses)),
        native_author_rows=sum(x['author_rows'] for x in coverage),
        duplicate_coordinate_rows_excluded=sum(x['duplicate_coordinate_rows_excluded'] for x in coverage),
        source_total_pattern_counts={name:sum(x['total_pattern_counts'][name] for x in coverage) for name in evidence.METRICS},
        all27_fixed_neighborhood_pair_slots=summary,
        all_declared_sources_accepted=len(all_effects)==len(scope['sources']),
        all_accepted_cell_scalar_patterns_and_full_target_extraction_verified=True,
        no_UMI_scale_or_physical_micron_or_patient_independence_or_treatment_phase_claim=True,
        peak_RSS_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*1024,
        individual_arrays_and_source_payload_RAM_only=True,completed_UTC=datetime.now(timezone.utc).isoformat())
    print(json.dumps(dict(kind='source_aggregate',value=result),allow_nan=False),flush=True)
