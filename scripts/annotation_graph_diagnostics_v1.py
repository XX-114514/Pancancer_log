"""Conditional native RNA/identity diagnostics on actual integration graphs.

Graphs are transductive. Expression units are not assumed to be patients.
Candidate identity strata are not independent biological truth.
"""
import numpy as np
import scipy.sparse as sp


def checked_graph(graph):
    w = graph.tocsr().astype(np.float64, copy=True)
    if w.shape[0] != w.shape[1] or not np.isfinite(w.data).all() or (w.data < 0).any():
        raise ValueError('Invalid graph weights or shape')
    w.sum_duplicates();w.setdiag(0);w.eliminate_zeros();w.sort_indices()
    return w


def unit_out_graph_state(graph, units, depth, values, minimum_neighbors=5,
                         minimum_train=30, return_predictions=False):
    """Weighted neighbor residual prediction, excluding the entire query unit."""
    w = checked_graph(graph)
    units, depth, values = np.asarray(units), np.asarray(depth, float), np.asarray(values, float)
    n = len(units)
    if w.shape != (n,n) or len(depth) != n or len(values) != n:
        raise ValueError('Graph-state axes differ')
    good = np.isfinite(depth) & np.isfinite(values)
    # Remove unavailable target cells and same-unit edges before all sums.
    row = np.repeat(np.arange(n), np.diff(w.indptr));col = w.indices
    retain = good[row] & good[col] & (units[row] != units[col])
    w.data[~retain] = 0;w.eliminate_zeros()
    degree = np.diff(w.indptr);weight = np.asarray(w.sum(axis=1)).ravel()
    neighbor_y = np.asarray(w @ np.where(good, values, 0)).ravel()
    neighbor_x = np.asarray(w @ np.where(good, depth, 0)).ravel()
    pred, base = np.full(n,np.nan), np.full(n,np.nan)
    fitted = 0
    for unit in np.unique(units[good]):
        train = good & (units != unit)
        query = good & (units == unit) & (degree >= minimum_neighbors) & (weight > 0)
        if train.sum() < minimum_train or not query.any():continue
        x = np.column_stack([np.ones(train.sum()),depth[train]])
        beta = np.linalg.lstsq(x,values[train],rcond=None)[0]
        base[query] = beta[0] + beta[1] * depth[query]
        pred[query] = base[query] + neighbor_y[query]/weight[query] - beta[0] - beta[1]*neighbor_x[query]/weight[query]
        fitted += 1
    ok = np.isfinite(pred)
    out = {'status':'insufficient_comparable_cells','n_input_cells':n,'n_evaluable':int(ok.sum()),
           'evaluable_fraction':float(ok.mean()) if n else 0.,
           'n_expression_units_evaluable_not_patients':fitted,
           'median_other_unit_neighbor_count':float(np.median(degree[ok])) if ok.any() else None,
           'same_unit_and_self_excluded':True,'graph_is_transductive':True}
    if ok.sum() >= 100:
        baseline_sse = float(np.sum((values[ok]-base[ok])**2))
        sse = float(np.sum((values[ok]-pred[ok])**2))
        out.update(status='estimated' if baseline_sse >= 1e-10 else 'no_variable_depth_residual',
                   depth_baseline_SSE=baseline_sse,neighbor_plus_depth_SSE=sse,
                   RNA_prediction_improvement_over_depth=1-sse/baseline_sse if baseline_sse>=1e-10 else None)
    if return_predictions:out.update(prediction=pred,baseline=base)
    return out


def direct_unit_out_graph_state(graph, units, depth, values, minimum_neighbors=5, minimum_train=30):
    """Independent centered OLS and row-wise dot products for audit."""
    w = graph.tocsr();units=np.asarray(units);depth=np.asarray(depth,float);values=np.asarray(values,float)
    good=np.isfinite(depth)&np.isfinite(values);pred=np.full(len(units),np.nan);base=pred.copy()
    for unit in np.unique(units[good]):
        train=good&(units!=unit)
        if train.sum()<minimum_train:continue
        x,y=depth[train],values[train];xc=x-x.mean();denom=float(xc@xc)
        slope=float(xc@(y-y.mean()))/denom if denom>1e-20 else 0.
        intercept=float(y.mean()-slope*x.mean())
        for i in np.flatnonzero(good&(units==unit)):
            low,high=w.indptr[i:i+2];j=w.indices[low:high];weight=w.data[low:high]
            take=(j!=i)&(units[j]!=unit)&good[j]&(weight>0)
            if take.sum()<minimum_neighbors:continue
            j,weight=j[take],weight[take]
            base[i]=intercept+slope*depth[i]
            pred[i]=base[i]+np.dot(weight,values[j]-intercept-slope*depth[j])/weight.sum()
    return pred,base


def identity_neighbor_profiles(graph, source, candidate, minimum_source_cells=20):
    w=checked_graph(graph);source=np.asarray(source);candidate=np.asarray(candidate)
    if w.shape!=(len(source),len(source)) or len(candidate)!=len(source):raise ValueError('Identity graph axes differ')
    total=np.asarray(w.sum(axis=1)).ravel();known=candidate!='unresolved'
    known_w=np.asarray(w@known.astype(float)).ravel();out=[]
    for label in np.unique(candidate[known]):
        members=candidate==label
        same=np.asarray(w@members.astype(float)).ravel()
        ds,counts=np.unique(source[members],return_counts=True)
        eligible_sources=ds[counts>=minimum_source_cells]
        # Only same-type neighbors from qualified sources contribute to effective source count.
        by_source=np.column_stack([np.asarray(w@(members&(source==s)).astype(float)).ravel() for s in eligible_sources]) if len(eligible_sources) else np.zeros((len(source),0))
        denom=by_source.sum(axis=1);effective=np.full(len(source),np.nan)
        valid=(denom>0)&np.isin(source,eligible_sources)
        if len(eligible_sources)>=2:
            probs=by_source[valid]/denom[valid,None];effective[valid]=1/np.sum(probs*probs,axis=1)
        for ds in np.unique(source[members]):
            query=members&(source==ds);connected=query&(total>0);known_connected=query&(known_w>0)
            eq=query&np.isfinite(effective)
            out.append({'source':str(ds),'candidate_type_not_truth':str(label),'n_candidate_cells':int(query.sum()),
                'n_connected':int(connected.sum()),'n_with_known_neighbor':int(known_connected.sum()),
                'same_candidate_fraction_among_known_neighbor_weight':float(np.mean(same[known_connected]/known_w[known_connected])) if known_connected.any() else None,
                'unknown_neighbor_weight_fraction':float(np.mean(1-known_w[connected]/total[connected])) if connected.any() else None,
                'n_qualified_sources_for_candidate':len(eligible_sources),
                'effective_source_status':'estimated_conditional_graph_composition' if eq.any() else 'insufficient_shared_type_coverage_or_neighbors',
                'n_effective_source_evaluable':int(eq.sum()),
                'same_candidate_neighbor_effective_sources':float(effective[eq].mean()) if eq.any() else None,
                'identity_accuracy_claim':False,'standard_iLISI_claim':False})
    return out


def rare_neighbor_profiles(graph, source, values):
    w=checked_graph(graph);source=np.asarray(source);values=np.asarray(values,float)
    if w.shape!=(len(source),len(source)) or len(values)!=len(source):raise ValueError('Rare graph axes differ')
    measured=np.isfinite(values);rare=measured&(values>0)
    denom=np.asarray(w@measured.astype(float)).ravel();num=np.asarray(w@rare.astype(float)).ravel()
    out=[]
    for ds in np.unique(source):
        q=(source==ds)&rare;ok=q&(denom>0);allmeasured=(source==ds)&measured
        out.append({'source':str(ds),'n_measured':int(allmeasured.sum()),'n_RNA_candidates':int(q.sum()),
            'n_candidates_with_measured_neighbors':int(ok.sum()),
            'candidate_prevalence':float(rare[allmeasured].mean()) if allmeasured.any() else None,
            'candidate_neighbor_weight_fraction':float(np.mean(num[ok]/denom[ok])) if ok.any() else None,
            'orthogonal_identity_validation':False})
    return out
