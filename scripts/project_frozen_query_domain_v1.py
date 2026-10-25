"""Whole-domain query adaptation with explicit archive and model-signal masks."""
import numpy as np
from frozen_symphony_reference_v1 import map_query,reference_fingerprint


def project_domain(reference,x,genes,source,domain,fit_pool,positive_RNA,mapper=map_query,symphony=None):
    n=x.shape[0];domain=np.asarray(domain,bool);fit_pool=np.asarray(fit_pool,bool);positive=np.asarray(positive_RNA,bool)
    if any(a.shape!=(n,) for a in [domain,fit_pool,positive]) or list(genes)!=reference.var_names.tolist():raise ValueError('Whole-source model/row axes differ')
    raw=np.full((n,reference.varm['PCs'].shape[1]),np.nan,np.float32);corrected=raw.copy();signal=np.asarray((x>0).getnnz(axis=1)>0).ravel()
    raw_rows=np.flatnonzero(domain&positive)
    for low in range(0,len(raw_rows),8192):
        q=raw_rows[low:low+8192];raw[q]=np.asarray(x[q]@reference.varm['PCs'])-reference.var['mean'].to_numpy()@reference.varm['PCs']
    numerical=domain&positive&np.isfinite(raw).all(axis=1)&(np.linalg.norm(np.nan_to_num(raw),axis=1)>1e-10)&(np.max(np.nan_to_num(raw,nan=-np.inf),axis=1)>0)
    adaptation_eligible=domain&fit_pool&positive&signal
    q=np.flatnonzero(adaptation_eligible);included=np.zeros(n,bool);before=reference_fingerprint(reference)
    # One call for the complete declared source/model domain. Only arithmetic
    # projection is chunked above; query-batch coefficients are not refit in chunks.
    if len(q):
        projected,mapped,supported=mapper(reference,x[q],genes,np.repeat(str(source),len(q)),symphony)
        np.testing.assert_allclose(projected,raw[q],rtol=2e-5,atol=2e-4)
        np.testing.assert_array_equal(supported,numerical[q]);corrected[q]=mapped;included[q[supported]]=True
    if reference_fingerprint(reference)!=before:raise ValueError('Whole-source adaptation mutated frozen reference')
    if np.any(included&~fit_pool) or np.any(included&~signal) or np.any(np.isfinite(corrected).any(axis=1)&~included):raise ValueError('Archive exclusion or zero model signal entered correction')
    if included.any() and not np.isfinite(corrected[included]).all():raise ValueError('Nonfinite supported whole-source correction')
    return {'raw_PCA':raw,'candidate_Harmony':corrected,'model_gene_signal':signal,'projection_numerical_supported':numerical,'adaptation_eligible':adaptation_eligible,'adaptation_included':included}


def project_raw_context(x,mean,loadings,positive_RNA):
    """Uncorrected context projection; no query batch is learned or removed."""
    positive=np.asarray(positive_RNA,bool);n=x.shape[0]
    if positive.shape!=(n,) or x.shape[1]!=len(mean) or loadings.shape[0]!=len(mean):raise ValueError('Raw-context projection axes differ')
    raw=np.full((n,loadings.shape[1]),np.nan,np.float32);q=np.flatnonzero(positive)
    for low in range(0,len(q),8192):
        take=q[low:low+8192];raw[take]=np.asarray(x[take]@loadings)-mean@loadings
    signal=np.asarray((x>0).getnnz(axis=1)>0).ravel()
    return {'raw_PCA':raw,'candidate_Harmony':np.full_like(raw,np.nan),'model_gene_signal':signal,'projection_numerical_supported':np.zeros(n,bool),'adaptation_eligible':np.zeros(n,bool),'adaptation_included':np.zeros(n,bool)}
