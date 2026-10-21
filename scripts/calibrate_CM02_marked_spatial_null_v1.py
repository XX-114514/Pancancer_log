"""Synthetic conditional-null calibration, no human input or selected effects.

Two noise mechanisms, fixed before observing source effects. This assesses the
residual-vector exchangeability approximation, separately from numerical tests.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,time
import numpy as np
import pandas as pd
from CM02_marked_spatial_core_v1 import (rankz,residualize,graph_edges,weighted_corr,
    permutation_groups,permutations,permuted_corr)


def wilson(k,n):
    z=1.959963984540054;p=k/n;d=1+z*z/n
    c=(p+z*z/(2*n))/d;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [float(c-h),float(c+h)]


def main(output):
    assert not output.exists();rng=np.random.default_rng(20261023)
    started=datetime.now(timezone.utc).isoformat();results=[]
    for mechanism in ['independent_homoscedastic','independent_heteroscedastic_spatial_variance']:
        pvalues=[]
        for replicate in range(200):
            observed=[];null=[]
            for case in range(30):
                while True:
                    xy=[rng.uniform(0,100,size=(30,2)) for _ in range(2)]
                    graph,reason=graph_edges(xy[0],xy[1],25)
                    if graph is not None:break
                marks=[];index=[]
                for points in xy:
                    bg=np.column_stack([rankz(points[:,0]),rankz(points[:,1]),
                        rankz(points[:,0]+points[:,1]+rng.normal(0,20,30)),rankz(rng.normal(size=30)),rankz(rng.normal(size=30))])
                    scale=np.ones(30) if mechanism=='independent_homoscedastic' else .25+3*(points[:,0]/100)**2
                    score=bg@np.array([5.,-3.,2.,4.,1.])+rng.normal(size=30)*scale
                    resid,status=residualize(score,bg);assert status=='ESTIMABLE';marks.append(resid[:,None])
                    frame=pd.DataFrame(points,columns=['Location_Center_X','Location_Center_Y'])
                    groups,n=permutation_groups(frame,(np.array([0,0]),np.array([100,100])),2)
                    assert n>=10;index.append(permutations(groups,30,199,rng))
                observed.append(float(weighted_corr(*marks,graph)[0,0]))
                null.append(permuted_corr(*marks,graph,*index)[:,0,0])
            actual=np.mean(observed);draws=np.mean(null,axis=0)
            p=min(1.,2*min((1+np.count_nonzero(draws>=actual))/200,(1+np.count_nonzero(draws<=actual))/200))
            pvalues.append(float(p))
            if (replicate+1)%25==0:
                output.with_name(output.stem+'_progress.json').write_text(json.dumps(dict(status='RUNNING_SYNTHETIC_NULL_CALIBRATION',
                    mechanism=mechanism,replicates_completed=replicate+1,observation_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
        rejected=sum(p<=.05 for p in pvalues);ci=wilson(rejected,200)
        results.append(dict(mechanism=mechanism,replicates=200,cases_per_replicate=30,cells_per_parent_per_ROI=30,
            permutations_per_replicate=199,alpha=.05,rejected=rejected,empirical_typeI=rejected/200,Wilson95=ci,
            coarse_calibration_flag='UPPER95_AT_MOST_0.10' if ci[1]<=.10 else 'APPROXIMATION_REQUIRES_CAUTION',
            inflation_lower95_exceeds_nominal=ci[0]>.05))
    output.write_text(json.dumps(dict(status='COMPLETE_SYNTHETIC_NULL_CALIBRATION',started_UTC=started,
        completed_UTC=datetime.now(timezone.utc).isoformat(),seed=20261023,results=results,
        no_human_or_source_effect_input=True,synthetic_calibration_not_proof_of_actual_residual_exchangeability=True),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);main(p.parse_args().output)
