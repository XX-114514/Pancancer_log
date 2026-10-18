"""Independent exact-rank/sign and complete random-stream generated-data checks."""
import itertools
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon,binom
from GSE176021_paired_tissue_math_v1 import wilcoxon_signed,rho,permutations,bootstrap_indices,bootstrap_effects,seed,adjust
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def corr(x,y):
    if min(np.ptp(x),np.ptp(y))==0:return np.nan
    a=pd.Series(x).rank().to_numpy();b=pd.Series(y).rank().to_numpy();a-=a.mean();b-=b.mean();return float(np.sum(a*b)/np.sqrt(np.sum(a*a)*np.sum(b*b)))


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);out=root/'GSE176021_paired_tissue_math_check_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);maxerr=0.;streams=0;bootstrap_count=0;exact_count=0
    vectors=[np.arange(1,11.),-np.arange(1,11.),np.array([0.,0.,1.,1.,-1.,2.,-2.,3.,-3.,4.]),np.zeros(10),np.array([-4.,-3.,-2.,-1.,1.,2.,3.,4.,5.,6.])]
    for d in vectors:
        z=wilcoxon_signed(d);nonzero=d[d!=0];r=pd.Series(np.abs(nonzero)).rank().to_numpy();observed=float(np.sum(np.where(nonzero>0,r,-r)));vals=np.array([np.sum(np.array(sign)*r) for sign in itertools.product([-1,1],repeat=len(r))]);assert z['P']==float((np.abs(vals)>=abs(observed)-1e-12).mean()) and z['statistic']==observed;assert np.array_equal(np.sort(vals),np.sort(z['statistics']));exact_count+=len(vals)
        count=(nonzero>0).sum();signp=min(1.,2*binom.cdf(min(count,len(r)-count),len(r),.5)) if len(r) else 1.;assert abs(signp-z['sign_P'])<1e-14
        if len(r)==10 and len(set(abs(nonzero)))==10:assert abs(wilcoxon(d,method='exact').pvalue-z['P'])<1e-14
    rng=np.random.default_rng(2026101906)
    for n in [10,11,15]:
        x=rng.normal(size=n);e=rng.normal(size=n);ix=bootstrap_indices(n);reference=np.random.default_rng(seed('all_primary_paired_cases')).integers(0,n,size=(2000,n));assert np.array_equal(ix,reference)
        for sign in [-1,1]:
            y=sign*x+.3*e;label=str(n)+'|'+str(sign);v=permutations(x,y,label);a=pd.Series(x).rank().to_numpy();b=pd.Series(y).rank().to_numpy();a-=a.mean();b-=b.mean();g=np.random.default_rng(seed(label+'|cochange'));vv=[]
            for start in range(0,19999,1000):
                index=np.argsort(g.random((min(1000,19999-start),n)),axis=1);vv.extend((np.sum(b[index]*a,axis=1)/np.sqrt(np.sum(a*a)*np.sum(b*b))).tolist())
            maxerr=max(maxerr,float(np.max(np.abs(v-np.array(vv)))));assert maxerr<1e-12;streams+=len(v)
            for cochange in [False,True]:
                got=bootstrap_effects(x,y,ix,cochange);want=np.array([corr(x[j],y[j]) if cochange else np.quantile(y[j],.5) for j in reference]);assert np.allclose(got,want,rtol=0,atol=1e-12,equal_nan=True);bootstrap_count+=len(got)
    assert np.isnan(rho(np.ones(10),np.arange(10)))
    from statsmodels.stats.multitest import multipletests
    for by in [False,True]:
        p=np.array([.001,.3,1.,.04,.9,.1]);assert np.allclose(adjust(p,by),multipletests(p,method='fdr_by' if by else 'fdr_bh')[1])
    receipt={'status':'ACCEPTED_EXACT_SIGNED_RANK_SIGN_TEST_TIES_ZEROS_AND_COMPLETE_GENERATED_RANDOM_STREAMS','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'math_sha256':sha(Path(__file__).with_name('GSE176021_paired_tissue_math_v1.py')),'all_exact_sign_statistics':exact_count,'all_permutation_statistics':streams,'all_bootstrap_effects':bootstrap_count,'maximum_random_statistic_error':maxerr,'new_data_or_effects_read':False}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
