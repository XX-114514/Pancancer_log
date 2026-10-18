"""Targeted small-case checks of reused HC3 and new rank permutations."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata,mannwhitneyu
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_math_v1 import design,independent_design,estimate,evaluate,wild_values,independent_wild,bootstrap_indices,bootstrap,seed,rank_effect,rank_permutations


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2);out=root/'immune_RNA_MCP_new_source_math_check_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    rng=np.random.default_rng(2026101904);maxerr=0.;wild_n=0;boot_n=0
    for slot,(kind,n,sign) in enumerate([('clinical',35,1),('clinical',35,-1),('coordination',40,1),('coordination',40,-1)],1):
        base=rng.normal(size=(n,8));C=np.c_[base,base[:,0],np.ones(n)];names=np.array(['q'+str(i) for i in range(C.shape[1])]);x=np.r_[np.zeros(20),np.ones(15)] if kind=='clinical' else np.repeat([0,1,2,3],[2,17,17,4]);y=sign*4*x+base@rng.normal(size=8)+rng.normal(size=n)
        m,status=design(x,y,C,names,kind);a,st=independent_design(x,y,C,names,kind);assert status==st=='ESTIMABLE'
        g,effect=estimate(m);alt=evaluate(a);assert sign*effect>0;assert np.allclose([effect,g['t']],[alt['effect'],alt['t']],atol=1e-10,rtol=0)
        w=wild_values(g,19999,seed(slot,'math_wild'));_,aw=independent_wild(a,19999,seed(slot,'math_wild'));assert np.allclose(w,aw,rtol=0,atol=1e-10);maxerr=max(maxerr,float(np.max(np.abs(w-aw))));wild_n+=len(w)
        strata=x.astype(int) if kind=='clinical' else np.zeros(n,int);ix=bootstrap_indices(strata,'math'+str(slot));b,reasons=bootstrap(x,y,C,names,kind,ix,m['nuisance_rank']);ab,ar=bootstrap(x,y,C,names,kind,ix,a['nuisance_rank'],independent=True);assert np.array_equal(reasons,ar) and np.allclose(b,ab,rtol=0,atol=1e-9,equal_nan=True);assert np.isfinite(b).sum()>=1900;maxerr=max(maxerr,float(np.nanmax(np.abs(b-ab))));boot_n+=len(b)
        np.savez_compressed(out/('private_synthetic_%d_full_stream.npz'%slot),wild=w,bootstrap=b,bootstrap_reasons=reasons,indices=ix)
        assert design(x[:19],y[:19],C[:19],names,kind)[0] is None
    # Independently compare rank effects and every permutation using pairwise
    # U comparisons (binary) and pandas ranks/einsum (ordinary Spearman).
    rank_n=0
    for binary in [False,True]:
        x=np.array([0,0,0,0,1,1,1,1]) if binary else np.array([0,1,1,2,3,4,4,5]);y=np.array([3,1,1,4,6,2,9,5.],float);sid='rank_binary' if binary else 'rank_ordinary';vals=rank_permutations(x,y,sid,binary);rg=np.random.default_rng(seed(sid,'rank_permutation'));other=[]
        yr=pd.Series(y).rank().to_numpy();xr=pd.Series(x).rank().to_numpy();xr-=xr.mean();centered=yr-yr.mean()
        for start in range(0,19999,1000):
            idx=np.argsort(rg.random((min(1000,19999-start),len(x))),axis=1);v=yr[idx]
            if binary:
                first=v[:,x==1];second=v[:,x==0];u=(first[:,:,None]>second[:,None,:]).sum(axis=(1,2))+.5*(first[:,:,None]==second[:,None,:]).sum(axis=(1,2));other.extend((2*u/(first.shape[1]*second.shape[1])-1).tolist())
            else:other.extend((np.einsum('ij,j->i',centered[idx],xr)/np.sqrt(np.sum(xr**2)*np.sum(centered**2))).tolist())
        assert np.allclose(vals,other,rtol=0,atol=1e-12);rank_n+=len(vals)
        if binary:
            expected=2*mannwhitneyu(y[x==1],y[x==0]).statistic/16-1;assert abs(expected-rank_effect(x,y,True))<1e-14
        else:assert abs(pd.Series(x).corr(pd.Series(y),method='spearman')-rank_effect(x,y))<1e-14
        np.savez_compressed(out/('private_'+sid+'_all_permutations.npz'),statistics=vals)
    x=rng.normal(size=35);C=np.c_[x,rng.normal(size=(35,10))];names=np.array(['c'+str(i) for i in range(C.shape[1])]);assert design(x,rng.normal(size=35),C,names,'clinical')[0] is None
    assert design(rng.normal(size=35),C[:,0],C,names,'clinical')[0] is None
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2
    receipt=dict(status='ACCEPTED_TARGETED_SMALL_CASE_FULL_HC3_BOOTSTRAP_AND_WILCOXON_SPEARMAN_PERMUTATION_MATH',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,files=files,wild_statistics_replayed=wild_n,bootstrap_effects_replayed=boot_n,rank_permutations_replayed=rank_n,maximum_error=maxerr,synthetic_only=True)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'),maximum_error=maxerr)))


if __name__=='__main__':main()
