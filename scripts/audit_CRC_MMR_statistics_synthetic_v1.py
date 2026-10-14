"""Independent synthetic checks; no patient or learned signature inputs."""
import numpy as np
import statsmodels.api as sm
from scipy.stats import mannwhitneyu

from test_CRC_DIALOGUE_MMR_cohort_v1 import hc3, rank_effect


def main():
    rng=np.random.default_rng(517)
    group=np.array([False]*15+[True]*15)
    y=np.r_[np.arange(15),np.arange(15)+20].astype(float)
    for value in [y,-y,np.r_[np.ones(15),np.ones(15)*2]]:
        u=mannwhitneyu(value[group],value[~group]).statistic
        assert abs(rank_effect(value,group)-(2*u/225-1))<1e-12
    assert rank_effect(y,group)==1 and rank_effect(-y,group)==-1
    cov=rng.normal(size=(30,4))
    X=np.column_stack((np.ones(30),group.astype(float),cov))
    outcome=.7*X[:,1]+cov[:,0]+rng.normal(size=30)
    actual,status=hc3(outcome,X)
    reference=sm.OLS(outcome,X).fit(cov_type='HC3',use_t=True)
    assert status=='estimated'
    assert abs(actual['beta']-reference.params[1])<1e-10
    assert abs(actual['HC3_p']-reference.pvalues[1])<1e-10
    duplicate=X.copy();duplicate[:,3]=duplicate[:,2]
    assert hc3(outcome,duplicate)[1]=='unidentifiable_design_no_covariate_dropped'
    assert hc3(outcome[:12],X[:12])[1]=='insufficient_patients_or_residual_df'
    print('PASS_SYNTHETIC_POSITIVE_NEGATIVE_TIES_HC3_COLLINEAR_AND_COVERAGE_GATES')


if __name__=='__main__':
    main()
