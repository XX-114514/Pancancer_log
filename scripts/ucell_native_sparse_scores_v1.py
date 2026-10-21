"""Independent positive-signature UCell math, capped average within-cell ranks.

Formula follows UCell >=2.7.6: 1-(sum(r)-m(m+1)/2)/(m*R-m(m+1)/2).
The caller provides measured, stable human features only; missing signature
features are excluded before calling, never imputed as zero measurement.
Implicit cell zeros of measured genes retain their full-universe average rank.
No gene names, learned weights, identities or arrays are saved by this module.
"""
import numpy as np


def score_positive_sparse_cell(positive_values,signature_values,universe_size,max_rank=1500):
    """Scores using the positive-count histogram, exact ties and implicit zeros.

positive_values includes every measured-human feature with count >0 once.
Each signature vector contains its measured feature counts, including cell
zeros. A canonical CSR row and unique gene mapping are caller invariants.
"""
    values=np.asarray(positive_values,dtype=float)
    assert values.ndim==1 and np.isfinite(values).all() and (values>0).all()
    assert (values==np.floor(values)).all() and values.max(initial=0)<2**53
    assert isinstance(universe_size,(int,np.integer)) and len(values)<=universe_size
    assert isinstance(max_rank,(int,np.integer)) and max_rank>0
    cutoff=min(int(max_rank),int(universe_size))
    unique,frequency=np.unique(values,return_counts=True)
    greater=np.cumsum(frequency[::-1])[::-1]-frequency
    average=greater+(frequency+1)/2
    zero_count=universe_size-len(values)
    zero_rank=len(values)+(zero_count+1)/2
    scores=[]
    for signature in signature_values:
        query=np.asarray(signature,dtype=float)
        assert query.ndim==1 and len(query)>0 and len(query)<=cutoff
        assert np.isfinite(query).all() and (query>=0).all() and (query==np.floor(query)).all()
        assert zero_count>0 or not (query==0).any()
        ranks=np.full(len(query),float(zero_rank))
        positive=query>0
        if positive.any():
            index=np.searchsorted(unique,query[positive])
            assert (index<len(unique)).all() and np.array_equal(unique[index],query[positive])
            ranks[positive]=average[index]
        ranks=np.minimum(ranks,cutoff)
        m=len(query);minimum=m*(m+1)/2;denominator=m*cutoff-minimum
        assert denominator>0
        result=1-(float(ranks.sum())-minimum)/denominator
        assert -1e-12<=result<=1+1e-12
        scores.append(float(np.clip(result,0,1)))
    return np.array(scores)


def score_selected_CSR_rows(indptr,indices,data,selected_rows,human,measured,program_columns,max_rank=1500):
    """Canonicalize selected raw rows and score only the declared native universe.

Each programme contains unique measured-human feature columns. Duplicate raw
CSR entries are summed before ranks/detection counts. No rank matrix is built.
First/last selected rows get an independent Python dictionary count replay.
"""
    def integer_indices(value):
        raw=np.asarray(value)
        assert raw.size==0 or np.issubdtype(raw.dtype,np.integer)
        return raw.astype(np.int64)
    ptr=integer_indices(indptr);columns=integer_indices(indices);values=np.asarray(data)
    chosen=integer_indices(selected_rows);human=np.asarray(human,dtype=bool);measured=np.asarray(measured,dtype=bool)
    assert ptr.ndim==1 and len(ptr)>0 and ptr[0]==0 and np.all(np.diff(ptr)>=0)
    assert ptr[-1]==len(columns)==len(values) and len(human)==len(measured)
    assert chosen.ndim==1 and np.all(np.diff(chosen)>0) and (chosen>=0).all() and (chosen<len(ptr)-1).all()
    assert not np.any(human & ~measured)
    programmes=[integer_indices(x) for x in program_columns]
    for p in programmes:
        assert p.ndim==1 and len(p)>0 and len(p)==len(np.unique(p))
        assert (p>=0).all() and (p<len(human)).all() and human[p].all()
    scores=np.zeros((len(chosen),len(programmes)));totals=np.zeros(len(chosen),dtype=np.int64);detected=totals.copy()
    for position,row in enumerate(chosen):
        a,b=int(ptr[row]),int(ptr[row+1]);idx=columns[a:b];val=values[a:b]
        assert (idx>=0).all() and (idx<len(human)).all()
        assert np.isfinite(val).all() and (val>=0).all() and (val==np.floor(val)).all()
        assert float(np.sum(val,dtype=np.float64))<2**53
        assert measured[idx[val>0]].all()
        take=human[idx] & (val>0);native_idx=idx[take];native_val=val[take].astype(np.int64)
        unique,inverse=np.unique(native_idx,return_inverse=True)
        canonical=np.zeros(len(unique),dtype=np.int64);np.add.at(canonical,inverse,native_val)
        totals[position]=canonical.sum();detected[position]=len(unique)
        query=[]
        for p in programmes:
            found=np.searchsorted(unique,p);ok=found<len(unique)
            good=np.zeros(len(p),dtype=bool)
            good[ok]=unique[found[ok]]==p[ok]
            counts=np.zeros(len(p),dtype=np.int64);counts[good]=canonical[found[good]];query.append(counts)
        scores[position]=score_positive_sparse_cell(canonical,query,int(human.sum()),max_rank)
        if position in [0,len(chosen)-1]:
            ordinary={}
            for feature,count in zip(idx,val):
                if human[feature] and count>0:ordinary[int(feature)]=ordinary.get(int(feature),0)+int(count)
            assert sum(ordinary.values())==totals[position] and len(ordinary)==detected[position]
            for p,q in zip(programmes,query):np.testing.assert_array_equal(q,[ordinary.get(int(k),0) for k in p])
    return scores,totals,detected
