"""Pure RAM diagnostics; additive raw-score covariance is not Spearman/causal.

Reference means are estimated within the same eligible case set, not from an
independent cohort. Zero observed subtype contribution when no captured cells
is an accounting identity; it does not impute unmeasured expression or assert
that the biological subtype is absent. No arrays or gene weights are saved.
"""
import numpy as np
import pandas as pd


def signature_diagnostics(positive,queries,universe,max_rank=1500):
    values=np.asarray(positive,float);assert np.isfinite(values).all() and (values>0).all()
    assert np.all(values==np.floor(values)) and len(values)<=universe
    cutoff=min(max_rank,universe);unique,freq=np.unique(values,return_counts=True)
    greater=np.cumsum(freq[::-1])[::-1]-freq;averages=greater+(freq+1)/2
    zero_rank=len(values)+(universe-len(values)+1)/2;scores=[];diagnostics=[]
    for query in queries:
        q=np.asarray(query,float);assert np.isfinite(q).all() and (q>=0).all() and np.all(q==np.floor(q))
        rank=np.full(len(q),float(zero_rank));positive_mask=q>0
        index=np.searchsorted(unique,q[positive_mask]);assert np.all(index<len(unique)) and np.array_equal(unique[index],q[positive_mask])
        rank[positive_mask]=averages[index];deficit=cutoff-np.minimum(rank,cutoff)
        mass=float(deficit.sum());squares=float(deficit@deficit);m=len(q);denom=m*cutoff-m*(m+1)/2;assert denom>0
        scores.append(mass/denom)
        diagnostics.append([int(positive_mask.sum()),int((deficit>0).sum()),mass*mass/squares if squares else 0.,float(deficit.max(initial=0)/mass) if mass else 0.])
    return np.array(scores),np.array(diagnostics)


def additive_decomposition(left_counts,right_counts,left_sums,right_sums,design,left_labels,right_labels):
    lc=np.asarray(left_counts);rc=np.asarray(right_counts);ls=np.asarray(left_sums,float);rs=np.asarray(right_sums,float);z=np.asarray(design,float)
    for counts,sums in [(lc,ls),(rc,rs)]:
        assert counts.shape==sums.shape and counts.ndim==2 and np.issubdtype(counts.dtype,np.integer) and (counts>=0).all()
        assert np.isfinite(sums).all() and (sums>=-1e-12).all() and (sums<=counts+1e-12).all()
        assert np.all(sums[counts==0]==0) and np.all(counts.sum(axis=1)>=20)
    assert len(lc)==len(rc)==len(z)>=20 and len(left_labels)==lc.shape[1] and len(right_labels)==rc.shape[1]
    assert np.linalg.matrix_rank(z)==z.shape[1] and np.linalg.cond(z)<=1e8 and len(z)-z.shape[1]>=10
    q=np.linalg.qr(z,mode='reduced')[0];maximum_error=0.
    def residual(array):
        nonlocal maximum_error
        r=array-q@(q.T@array);other=array-z@np.linalg.lstsq(z,array,rcond=None)[0]
        np.testing.assert_allclose(r,other,atol=1e-12,rtol=1e-10)
        maximum_error=max(maximum_error,float(np.max(np.abs(r-other),initial=0)))
        return r-r.mean(axis=0)
    left_terms=ls/lc.sum(axis=1)[:,None];right_terms=rs/rc.sum(axis=1)[:,None]
    left=left_terms.sum(axis=1);right=right_terms.sum(axis=1)
    lr=residual(left);rr=residual(right);scale=float(np.sqrt(np.mean(lr**2)*np.mean(rr**2)))
    if scale<=1e-20:return dict(status='RAW_SCORE_RESIDUAL_NOT_ESTIMABLE',reason='NO_RESIDUAL_RAW_SCORE_VARIATION'),[]
    covariance=float(np.mean(lr*rr));raw_rho=covariance/scale;references=[]
    for scheme in ['pooled_cell_weighted_reference','equal_observed_case_reference']:
        pieces=[]
        for counts,sums,total in [(lc,ls,left),(rc,rs,right)]:
            if scheme=='pooled_cell_weighted_reference':
                denom=counts.sum(axis=0);mu=np.divide(sums.sum(axis=0),denom,out=np.zeros(counts.shape[1]),where=denom>0)
            else:
                means=np.divide(sums,counts,out=np.zeros_like(sums),where=counts>0);observed=(counts>0).sum(axis=0)
                mu=np.divide(means.sum(axis=0),observed,out=np.zeros(counts.shape[1]),where=observed>0)
            weights=counts/counts.sum(axis=1)[:,None]
            between=weights@mu;within=total-between
            r=residual(np.column_stack([between,within]));np.testing.assert_allclose(r.sum(axis=1),residual(total),atol=1e-12,rtol=1e-10)
            pieces.append(r)
        terms=[]
        for i,lname in enumerate(['composition_reference','weighted_within_fine_deviation']):
            for j,rname in enumerate(['composition_reference','weighted_within_fine_deviation']):
                value=float(np.mean(pieces[0][:,i]*pieces[1][:,j]));terms.append(dict(left_component=lname,right_component=rname,
                    residual_covariance=value,normalized_covariance_contribution=value/scale))
        error=abs(sum(t['residual_covariance'] for t in terms)-covariance);assert error<=1e-12
        references.append(dict(reference=scheme,terms=terms,additive_covariance_absolute_error=error,
            not_causal_mediation_or_percent_variance=True,same_case_reference_not_independent=True))
    lres=residual(left_terms);rres=residual(right_terms);fine=[]
    for i,lname in enumerate(left_labels):
        for j,rname in enumerate(right_labels):
            value=float(np.mean(lres[:,i]*rres[:,j]));fine.append(dict(left_author_fine=lname,right_author_fine=rname,
                residual_covariance=value,normalized_covariance_contribution=value/scale,
                cases_left_fine_atleast20=int((lc[:,i]>=20).sum()),cases_right_fine_atleast20=int((rc[:,j]>=20).sum()),
                cases_both_fine_atleast20=int(((lc[:,i]>=20)&(rc[:,j]>=20)).sum()),
                observed_zero_term_not_true_subtype_absence=True))
    assert abs(sum(t['residual_covariance'] for t in fine)-covariance)<=1e-12
    np.testing.assert_allclose(sum(t['normalized_covariance_contribution'] for t in fine),raw_rho,atol=1e-10)
    return dict(status='DESCRIPTIVE_RAW_SCORE_COVARIANCE_DECOMPOSED',n_cases=len(lc),
        residual_raw_score_Pearson=raw_rho,residual_raw_score_covariance=covariance,
        reference_decompositions=references,independent_LS_maximum_absolute_error=maximum_error,
        fine_pair_covariances_add_to_total=True,not_a_decomposition_of_primary_Spearman=True,
        no_new_P_CI_or_clinical_test=True),fine


def support_profiles(tab,frame,nodes,diagnostics,inputs):
    assert diagnostics.shape==(len(frame),2,4) and np.isfinite(diagnostics).all()
    f=frame.copy();registry=tab.set_index('library_id');f['case']=f.library_id.map(registry.patient_registry)
    f['parent']=f.label.map(inputs.parent);d=diagnostics[:,1]
    for j,name in enumerate(['detected','active','effective','maximum_fraction']):f[name]=d[:,j]
    f.loc[f.active.eq(0),'maximum_fraction']=np.nan
    for k in [1,2,5,10]:f['fraction_active_atleast'+str(k)]=(f.active>=k).astype(float)
    metrics=['detected','active','effective','maximum_fraction']+['fraction_active_atleast'+str(k) for k in [1,2,5,10]]
    output=[]
    for context in inputs.CONTEXTS:
        selected=f if context==inputs.CONTEXTS[0] else f.loc[f.node_id.isin(nodes)]
        for phase in ['pre','post']:
            current=selected.loc[selected.timepoint.eq(phase)]
            for kind,labels in [('parent',inputs.PARENTS),('fine',inputs.REGULATORY+inputs.CONVENTIONAL+inputs.CD8)]:
                for label in labels:
                    group=current.loc[current['parent' if kind=='parent' else 'label'].eq(label)]
                    counts=group.groupby('case').size();qualified=counts[counts>=20].index
                    means=group.loc[group.case.isin(qualified)].groupby('case')[metrics].mean()
                    summaries={}
                    if len(qualified)>=5:
                        for name in metrics:
                            v=means[name].dropna().to_numpy()
                            summaries[name]=dict(defined_case_count=len(v),case_mean_min_q25_median_q75_max=np.quantile(v,[0,.25,.5,.75,1]).tolist() if len(v)>=5 else None)
                    output.append(dict(context=context,phase=phase,reference_kind=kind,reference_label=label,
                        cases_atleast20_observed_cells=len(qualified),source_cell_observations=len(group),metrics=summaries,
                        minimum5_cases_for_distribution=True,maximum_fraction_conditioned_on_positive_score_cells=True,
                        gene_thresholds_descriptive_not_new_state_labels=True))
    assert len(output)==64
    return output
