"""Use unchanged CellTypist functions with lossless Unicode gene metadata.

NumPy's object-array membership path repeatedly compares each model feature.
Only the representation of gene-name arrays is changed. The official
Classifier constructor, scaling, feature ordering and celltype method remain
unchanged; expression and model numerical parameters are not replaced.
"""
import copy
import numpy as np
import celltypist
from celltypist.classifier import Classifier

def lossless_unicode(values):
    original=np.asarray(values)
    if original.ndim!=1 or not all(isinstance(v,str) for v in original.tolist()):
        raise ValueError('Expected a one-dimensional exact string gene axis')
    converted=np.asarray(original.tolist(),dtype=np.str_)
    np.testing.assert_array_equal(converted,original)
    if converted.tolist()!=original.tolist():raise ValueError('Gene names changed')
    return converted

def official_annotate_unicode_axis(adata,model,mode='best match',p_thres=.5):
    if celltypist.__version__!='1.7.1':raise ValueError('Revalidate official API version before changing runtime')
    # annotate(... majority_voting=False) uses exactly these two calls.
    # Clone only the model to preserve its public cache and caller's object.
    classifier=Classifier(filename=adata,model=copy.deepcopy(model),transpose=False,gene_file=None,cell_file=None)
    classifier.indata_genes=lossless_unicode(classifier.indata_genes)
    classifier.model.classifier.features=lossless_unicode(classifier.model.classifier.features)
    return classifier.celltype(mode=mode,p_thres=p_thres)

def independent_exact_membership_check(query,features):
    a=lossless_unicode(query);b=lossless_unicode(features);names=set(b.tolist())
    expected=np.fromiter((g in names for g in a.tolist()),dtype=bool,count=len(a))
    np.testing.assert_array_equal(np.isin(a,b),expected)
    return expected
