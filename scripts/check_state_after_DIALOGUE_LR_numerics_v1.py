import importlib.util,numpy as np
from statsmodels.stats.multitest import multipletests
p='/tmp/test_state_after_DIALOGUE_LR_response_v1.py';spec=importlib.util.spec_from_file_location('fixed_lr',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);rng=np.random.default_rng(61142);errors=[]
for n in [25,40,60]:
 c=np.column_stack([np.ones(n),rng.normal(size=(n,8))]);x=rng.normal(size=(n,3));a=m.scorr(x,x[:,0],c)[0];b=m.scorr(x,-x[:,0],c)[0];assert abs(a[0]-1)<1e-12 and abs(b[0]+1)<1e-12;assert np.isnan(m.scorr(x,np.ones(n),c)[0]).all()
 y=rng.normal(size=n);v,ex,ey=m.scorr(x,y,c);rx=m.rankdata(x,axis=0);ry=m.rankdata(y)
 for j in range(3):
  ax=rx[:,j]-c@np.linalg.lstsq(c,rx[:,j],rcond=None)[0];ay=ry-c@np.linalg.lstsq(c,ry,rcond=None)[0];ref=ax@ay/(np.linalg.norm(ax)*np.linalg.norm(ay));errors.append(abs(ref-v[j]));assert abs(ref-v[j])<1e-12
 assert np.allclose(m.residual(rx-ex+ex,c),ex,atol=1e-12)
pv=np.array([.1,.001,.03,np.nan,.99,.015]);finite=np.isfinite(pv);assert np.allclose(m.bh(pv)[finite],multipletests(pv[finite],method='fdr_bh')[1]);assert np.isnan(m.bh(pv)[~finite]).all();print('PASS_POSITIVE_NEGATIVE_CONSTANT_QR_OLS_IDENTITY_PERMUTATION_MISSING_BH',max(errors))
