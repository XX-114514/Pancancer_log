"""Audit source calculations without changing statistics, draws or outputs."""
import numpy as np


def attach(core):
    counts=dict(covariance_calls=0,family_positions=0,bootstrap_means=0,permutation_source_means=0)
    original_covariance=core.edge_covariance
    def checked_covariance(left,right,a,b):
        result=original_covariance(left,right,a,b)
        x,y=left[a],right[b]
        independently_centered=np.einsum('ij,ik->jk',x-x.mean(0),y-y.mean(0))/len(a)
        np.testing.assert_allclose(result,independently_centered,rtol=1e-9,atol=1e-10)
        counts['covariance_calls']+=1
        return result
    core.edge_covariance=checked_covariance
    original_finish=core.finish_family
    def quantile(values,p):
        ordered=sorted(values);index=(len(ordered)-1)*p;low=int(index);fraction=index-low
        return ordered[low] if low==len(ordered)-1 else ordered[low]*(1-fraction)+ordered[low+1]*fraction
    def checked_finish(observed,expected,draws,roi_counts,family,seed):
        result=original_finish(observed,expected,draws,roi_counts,family,seed)
        probabilities=[]
        for index,row in enumerate(result):
            selected=[k for k in range(len(roi_counts)) if roi_counts[k,index]>0]
            n=len(selected);assert n==row['source_case_units']
            if n<core.SOURCE_MIN:
                assert row['effect'] is None and row['P'] is None and row['BH'] is None
                probabilities.append(1.);counts['family_positions']+=1;continue
            values=[float(observed[k,index]-expected[k,index])/int(roi_counts[k,index]) for k in selected]
            point=sum(values)/n
            np.testing.assert_allclose(row['effect'],point,rtol=1e-9,atol=1e-10)
            null=[sum(float(draws[k,j,index]-expected[k,index])/int(roi_counts[k,index]) for k in selected)/n for j in range(draws.shape[1])]
            p=(1+sum(abs(value)>=abs(point)-1e-12 for value in null))/(len(null)+1)
            assert abs(p-row['P'])<1e-12
            counts['permutation_source_means']+=len(null)
            probabilities.append(p)
            if n>=10:
                rng=np.random.default_rng(seed+index)
                indices=rng.integers(0,n,size=(core.BOOTSTRAPS,n))
                boot=[sum(values[int(k)] for k in sample)/n for sample in indices]
                np.testing.assert_allclose(row['case_bootstrap_95_CI'],[quantile(boot,.025),quantile(boot,.975)],rtol=1e-9,atol=1e-10)
                loo=[sum(values[j] for j in range(n) if j!=k)/(n-1) for k in range(n)]
                np.testing.assert_allclose(row['leave_one_case_effect_range'],[min(loo),max(loo)],rtol=1e-9,atol=1e-10)
                counts['bootstrap_means']+=len(boot)
            else:assert row['case_bootstrap_95_CI'] is None
            assert row['positive_case_count']==sum(value>0 for value in values)
            assert row['negative_case_count']==sum(value<0 for value in values)
            counts['family_positions']+=1
        order=sorted(range(len(probabilities)),key=lambda k:probabilities[k]);q=[1.]*len(order);bound=1.
        for position in range(len(order)-1,-1,-1):
            k=order[position];bound=min(bound,probabilities[k]*len(order)/(position+1));q[k]=bound
        for index,row in enumerate(result):
            if row['P'] is not None:assert abs(row['BH']-q[index])<1e-12
        return result
    core.finish_family=checked_finish
    return counts
