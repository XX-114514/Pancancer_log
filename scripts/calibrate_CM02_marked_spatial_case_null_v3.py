"""Larger pre-source calibration, retaining every initial synthetic trial.

Fixed 10 batches x 200 experiments per mechanism, no stopping when a flag
passes. Both prior 200-trial counts are independently replayed in batch one.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,time
import numpy as np
import pandas as pd
from CM02_marked_spatial_core_v1 import (rankz,residualize,graph_edges,weighted_corr,
    permutation_groups,permutations,permuted_corr)


def interval(k,n,z=2.5758293035489004):
    p=k/n;d=1+z*z/n;c=(p+z*z/(2*n))/d;h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [float(c-h),float(c+h)]


def experiment(rng,hetero):
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
            scale=.25+3*(points[:,0]/100)**2 if hetero else np.ones(30)
            score=bg@np.array([5.,-3.,2.,4.,1.])+rng.normal(size=30)*scale
            resid,status=residualize(score,bg);assert status=='ESTIMABLE';marks.append(resid[:,None])
            frame=pd.DataFrame(points,columns=['Location_Center_X','Location_Center_Y'])
            groups,n=permutation_groups(frame,(np.array([0,0]),np.array([100,100])),2)
            assert n>=10;index.append(permutations(groups,30,199,rng))
        observed.append(float(weighted_corr(*marks,graph)[0,0]));null.append(permuted_corr(*marks,graph,*index)[:,0,0])
    return np.asarray(observed),np.asarray(null)


def main(output,prior):
    assert not output.exists();legacy=json.loads(prior.read_text());rng=np.random.default_rng(20261023)
    started=datetime.now(timezone.utc).isoformat();monotonic=time.monotonic();flips=[]
    for i in range(30):
        seed=int.from_bytes(hashlib.sha256(('20261022:'+f'synthetic_{i}').encode()).digest()[:8],'little')
        flips.append(np.random.default_rng(seed).choice([-1.,1.],2000))
    flips=np.column_stack(flips)
    mechanisms=['independent_homoscedastic','independent_heteroscedastic_spatial_variance']
    cells={k:[] for k in mechanisms};cases={k:[] for k in mechanisms};moments={k:[] for k in mechanisms};batches=[]
    for batch in range(10):
        for mechanism in mechanisms:
            before_cell=len(cells[mechanism]);before_case=len(cases[mechanism])
            for replicate in range(200):
                observed,null=experiment(rng,mechanism!=mechanisms[0]);actual=observed.mean();draws=null.mean(axis=0)
                cells[mechanism].append(float(min(1.,2*min((1+np.count_nonzero(draws>=actual))/200,(1+np.count_nonzero(draws<=actual))/200))))
                delta=observed-null.mean(axis=1);moments[mechanism].append(float(delta.mean()))
                cases[mechanism].append(float((1+np.count_nonzero(np.abs(flips@delta/30)>=abs(delta.mean())-1e-12))/2001))
            cell_k=sum(p<=.05 for p in cells[mechanism][before_cell:]);case_k=sum(p<=.05 for p in cases[mechanism][before_case:])
            if batch==0:
                old=next(x for x in legacy['results'] if x['mechanism']==mechanism)
                assert cell_k==old['rejected'] and case_k==old['patient_signflip_rejected']
            batches.append(dict(batch=batch+1,mechanism=mechanism,trials=200,cell_rejected=cell_k,case_rejected=case_k))
            output.with_name(output.stem+'_progress.json').write_text(json.dumps(dict(status='RUNNING_EXTENDED_SYNTHETIC_CALIBRATION',
                completed_batches=len(batches),total_batches=20,observation_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
            assert time.monotonic()-monotonic<15*60
    results=[]
    for mechanism in mechanisms:
        c=sum(p<=.05 for p in cells[mechanism]);k=sum(p<=.05 for p in cases[mechanism]);ci=interval(k,2000)
        results.append(dict(mechanism=mechanism,replicates=2000,cell_residual_rejected=c,cell_residual_typeI=c/2000,
            patient_signflip_rejected=k,patient_signflip_typeI=k/2000,patient_signflip_Wilson99=ci,
            patient_signflip_calibration_flag='QUALIFIED_EXTENDED_SYNTHETIC_SCOPE' if ci[0]<=.05 and ci[1]<=.075 else 'REQUIRES_CAUTION',
            mean_case_excess_over_simulations=float(np.mean(moments[mechanism])),
            alpha_cdf={str(a):float(np.mean(np.asarray(cases[mechanism])<=a)) for a in [.01,.05,.10]}))
    output.write_text(json.dumps(dict(status='COMPLETE_EXTENDED_SYNTHETIC_CASE_NULL_CALIBRATION',started_UTC=started,
        completed_UTC=datetime.now(timezone.utc).isoformat(),seed=20261023,results=results,batches=batches,
        all_initial_200_trials_retained_and_replayed=True,no_human_input_or_real_effect_selection=True,
        synthetic_calibration_not_actual_case_independence_or_symmetry_proof=True),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--prior',type=Path,required=True);a=p.parse_args();main(a.output,a.prior)
