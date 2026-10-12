"""Known-coefficient and rank-failure checks; never biological validation."""
from pathlib import Path
import importlib.util,json
import numpy as np
p=Path(__file__).with_name('fullgene_receiver_target_association_v1.py');s=importlib.util.spec_from_file_location('target_math',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
rng=np.random.default_rng(172620);a=rng.normal(size=(80,8));plan={'resampling':{'repeats':2000,'seed':20261012,'minimum_valid_bootstrap':1800}}
for sign in [1,-1]:
 y=sign*2*a[:,0]+.4*a[:,1];r,X,Y=m.evaluate(a,y,'known_'+str(sign),plan);expected=sign*2*a[:,0].std(ddof=1)/y.std(ddof=1)
 assert np.isclose(r['sender_beta'],expected,atol=1e-12);assert r['residual_permutation_p']<.01;assert r['bootstrap_valid']>=1800
 print(json.dumps({'fixture':'exact_standardized_beta','sign':sign,'beta':r['sender_beta'],'expected':expected,'status':'PASS'}))
bad=a.copy();bad[:,2]=a[:,0];r=m.evaluate(bad,a[:,0],'collinear',plan);assert r['status']=='REQUIRED_DESIGN_NOT_FULL_RANK';print('PASS: collinear required control rejected, no coefficient removal')
