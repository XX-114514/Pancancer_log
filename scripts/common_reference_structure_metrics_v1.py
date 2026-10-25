"""Fixed native-state structure comparisons with unit exclusion and common queries."""
import numpy as np
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors
from annotation_graph_diagnostics_v1 import unit_out_graph_state,direct_unit_out_graph_state
from paired_granularity_review_v1 import paired_scores


def exact_neighbor_graph(z,neighbors=30):
    z=np.asarray(z)
    if z.ndim!=2 or len(z)<2 or not np.isfinite(z).all():raise ValueError('Invalid finite structure coordinates')
    k=min(neighbors,len(z)-1)
    _,near=NearestNeighbors(n_neighbors=k+1,n_jobs=4).fit(z).kneighbors(z)
    near=np.asarray([r[r!=i][:k] for i,r in enumerate(near)])
    if near.shape!=(len(z),k):raise ValueError('Exact self exclusion failed')
    return sp.csr_matrix((np.ones(len(z)*k),near.ravel(),np.arange(0,len(z)*k+1,k)),shape=(len(z),len(z)))


def paired_native_states(coordinates,units,depth,values,names,minimum_neighbors=5,minimum_train=30,audit=False,graphs=None):
    if values.shape!=(len(units),len(names)):raise ValueError('State matrix axes differ')
    results=[];pairs=[];cache={};max_error=0.
    for variant,z in coordinates.items():
        graph=exact_neighbor_graph(z) if graphs is None else graphs[variant];cache[variant]={}
        for j,name in enumerate(names):
            item=unit_out_graph_state(graph,units,depth,values[:,j],minimum_neighbors,minimum_train,True);cache[variant][name]=item
            results.append({'variant':variant,'state':name,**{k:v for k,v in item.items() if k not in ['prediction','baseline']}})
            if audit:
                pred,base=direct_unit_out_graph_state(graph,units,depth,values[:,j],minimum_neighbors,minimum_train)
                np.testing.assert_allclose(pred,item['prediction'],atol=1e-9,rtol=1e-7,equal_nan=True);np.testing.assert_allclose(base,item['baseline'],atol=1e-9,rtol=1e-7,equal_nan=True)
                finite=np.isfinite(pred)&np.isfinite(item['prediction'])
                if finite.any():max_error=max(max_error,float(np.max(np.abs(pred[finite]-item['prediction'][finite]))))
    variants=list(coordinates)
    for i,a in enumerate(variants):
        for b in variants[i+1:]:
            for j,name in enumerate(names):pairs.append({'first_variant':a,'second_variant':b,'state':name,**paired_scores(values[:,j],units,cache[a][name],cache[b][name])})
    return results,pairs,max_error
