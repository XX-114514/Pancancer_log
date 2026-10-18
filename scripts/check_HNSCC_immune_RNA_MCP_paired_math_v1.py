"""Complete synthetic positive/negative/gate and random-stream acceptance."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.stats import binom, binomtest

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from HNSCC_immune_RNA_MCP_paired_math_v1 import design, effect, estimate, seed, wild_values, patient_takes, bootstrap, median_shift
from HNSCC_immune_RNA_MCP_paired_audit_math_v1 import reference_model, reference_estimate, reference_wild, reference_bootstrap


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    out=root/'HNSCC_immune_RNA_MCP_paired_math_gate_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_complete_synthetic_math.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,scope='Pure26-subject positive/negative/collinear/constant/nonidentifiable math; all19999 wild/2000 bootstrap per four estimable scenarios. No real paired scores loaded.'),indent=2)+'\n')
    rng=np.random.default_rng(2026101903);names=['treatment_DplusT','HPV_positive']+['q'+str(i) for i in range(8)]
    C=np.column_stack([np.arange(26)%2,(np.arange(26)%5==0).astype(int),rng.normal(size=(26,8))]);x=rng.normal(size=26);rows=[];maxwild=0.;maxboot=0.
    for slot,(direction,duplicate) in enumerate([(1,False),(-1,False),(1,True),(-1,True)]):
        controls=C.copy()
        if duplicate:controls[:,-1]=controls[:,-2]
        y=direction*.85*x+.3*controls[:,2]+rng.normal(size=26)*.3
        m,status=design(x,y,controls,names);ref,rs=reference_model(x,y,controls,names);assert status==rs=='ESTIMABLE'
        g=estimate(m);b,t,r=reference_estimate(ref);assert abs(g['beta']-b)<1e-9 and abs(g['t']-t)<1e-8 and abs(effect(m)-r)<1e-10
        wild=wild_values(g,19999,seed(slot,'wild'));rw=reference_wild(ref['y'],ref['D'],19999,seed(slot,'wild'));err=float(np.max(np.abs(wild-rw)));assert err<1e-8;maxwild=max(maxwild,err)
        indices=patient_takes(26,slot);bv,bs=bootstrap(x,y,controls,names,indices,m['nuisance_rank']);rv,rbs=reference_bootstrap(x,y,controls,names,indices,ref['nuisance_rank']);assert np.array_equal(bs,rbs) and np.allclose(bv,rv,atol=1e-9,rtol=0,equal_nan=True)
        be=float(np.nanmax(np.abs(bv-rv)));maxboot=max(maxboot,be)
        np.savez_compressed(out/('private_synthetic_slot%02d_streams.npz'%slot),x=x,y=y,C=controls,names=names,wild=wild,reference_wild=rw,bootstrap_indices=indices,bootstrap=bv,reference_bootstrap=rv,statuses=bs)
        rows.append(dict(slot=slot,status=status,direction=direction,duplicate_nuisance=duplicate,maximum_wild_error=err,maximum_bootstrap_error=be))
    for label,xx,yy,cc in [('lowN',x[:19],x[:19],C[:19]),('constant',x,np.ones(26),C),('target_in_controls',C[:,2],x,C)]:
        m,s=design(xx,yy,cc,names);r,rs=reference_model(xx,yy,cc,names);assert m is None and r is None and s==rs;rows.append(dict(scenario=label,status=s))
    for d in [np.arange(26)-5.,-(np.arange(26)-5.),np.r_[np.zeros(6),np.arange(1,21)],np.r_[-np.ones(13),np.ones(13)]]:
        s,status=median_shift(d);assert status=='ESTIMABLE';used=d[np.abs(d)>1e-10];p=binomtest(int((used>0).sum()),len(used),.5).pvalue;assert abs(s['p']-p)<1e-14
        N=len(d);k=max(k for k in range(1,N//2+1) if 2*binom.cdf(k-1,N,.5)<=.05);ordered=np.sort(d);assert s['CI_low']==ordered[k-1] and s['CI_high']==ordered[N-k]
    assert median_shift(np.zeros(26))[1]=='FEWER_THAN10_NONZERO_PAIRED_DELTAS'
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    t=dict(status='ACCEPTED_COMPLETE_PAIRED_SYNTHETIC_MATH_SIGN_INTERVALS_AND_ALL_RANDOM_VALUES',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,scenarios=rows,all_wild_values_replayed=79996,all_bootstrap_values_replayed=8000,maximum_wild_error=maxwild,maximum_bootstrap_error=maxboot,no_real_paired_scores_or_effects=True)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
