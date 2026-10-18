"""Pure all-layer HDF5 numerical/structural checks shared across locations.

Sharing this kernel is not independent algorithmic validation.
"""
import h5py
import numpy as np

def string(x):
    return x.decode('utf-8') if isinstance(x,bytes) else str(x)


def matrix_structure(group):
    if isinstance(group,h5py.Dataset):return dict(storage='dense',shape=list(group.shape),dtype=str(group.dtype))
    encoding=string(group.attrs.get('encoding-type','unknown'))
    assert encoding in ['csr_matrix','csc_matrix'], 'Unsupported matrix encoding: '+encoding
    shape=list(map(int,group.attrs['shape']));assert len(shape)==2
    ptr=np.asarray(group['indptr'],np.int64);axis=0 if encoding=='csr_matrix' else 1
    assert len(ptr)==shape[axis]+1 and ptr[0]==0 and np.all(np.diff(ptr)>=0) and ptr[-1]==len(group['data'])==len(group['indices'])
    limit=shape[1-axis]
    for start in range(0,len(group['indices']),1000000):
        a=np.asarray(group['indices'][start:start+1000000]);assert (a>=0).all() and (a<limit).all()
    return dict(storage=encoding,shape=shape,dtype=str(group['data'].dtype),nonzero_stored=int(ptr[-1]),empty_axis_entries=int((np.diff(ptr)==0).sum()))


def numerical(group):
    d=group if isinstance(group,h5py.Dataset) else group['data']
    minimum=np.inf;maximum=-np.inf;negative=zeros=noninteger=total=0;sum_values=0.
    if d.ndim==1:
        chunks=(np.asarray(d[start:start+1000000]) for start in range(0,len(d),1000000))
    else:chunks=(np.asarray(d[start:start+128,:]) for start in range(0,d.shape[0],128))
    for a in chunks:
        assert np.isfinite(a).all(),'Nonfinite matrix value'
        total+=a.size
        if a.size:minimum=min(minimum,float(a.min()));maximum=max(maximum,float(a.max()))
        negative+=int((a<0).sum());zeros+=int((a==0).sum());noninteger+=int((a!=np.floor(a)).sum());sum_values+=float(a.astype(float).sum())
    return dict(all_finite=True,stored_values=total,min=None if total==0 else minimum,max=None if total==0 else maximum,negative_stored_values=negative,zero_stored_values=zeros,noninteger_stored_values=noninteger,sum_stored_values=sum_values,
                integer_values_alone_do_not_prove_raw=True)


