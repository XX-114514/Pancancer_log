"""Independent numerical replay of all estimable paired tests and BH."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from audit_HNSCC_fixed_MCP_clinical_v1 import ranks


def corr(x,y):
    x=ranks(x);y=ranks(y);x=x-x.mean();y=y-y.mean()
    d=np.linalg.norm(x)*np.linalg.norm(y)
    return float(x@y/d) if d>1e-12 else np.nan


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs=root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1'
    out=root/'HNSCC_fixed_MCP_paired_results_20261014_v1'
    r=pd.read_csv(out/'all12_private_paired_statistics.tsv',sep='\t')
    assert len(r)==12 and not r.duplicated(['program','question']).any()
    scores=pd.read_csv(inputs/'private_patient_library_frozen_scores.tsv',sep='\t')
    pre=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='CRC_training_z')
    post=scores[scores.timepoint.eq('post')].pivot(index='patient_id',columns='endpoint',values='CRC_training_z')
    ix=pre.index.intersection(post.index).sort_values();pre=pre.loc[ix];post=post.loc[ix]
    questions=dict(paired_shift_T=((post.Tnode9-pre.Tnode9).to_numpy(),None),
      paired_shift_Mac=((post.Macnode18-pre.Macnode18).to_numpy(),None),
      paired_shift_equal_mean=((post.mean(axis=1)-pre.mean(axis=1)).to_numpy(),None),
      pre_cross_type=(pre.Tnode9.to_numpy(),pre.Macnode18.to_numpy()),
      post_cross_type=(post.Tnode9.to_numpy(),post.Macnode18.to_numpy()),
      change_cross_type=((post.Tnode9-pre.Tnode9).to_numpy(),(post.Macnode18-pre.Macnode18).to_numpy()))
    for row in r.itertuples():
        if row.program=='MCP2':assert pd.isna(row.p) and row.status=='unsupported_original_program_slot';continue
        x,y=questions[row.question];n=len(x);assert n==row.n_paired_patients
        effect=float(np.sum(x)/n) if y is None else corr(x,y)
        assert np.isclose(effect,row.effect,atol=1e-12)
        rng=np.random.default_rng(20261014+int(hashlib.sha256(row.question.encode()).hexdigest()[:8],16))
        extreme=0;boots=[]
        for _ in range(20000):
            trial=float(x@rng.choice([-1,1],n)/n) if y is None else corr(x,y[rng.permutation(n)])
            extreme+=abs(trial)>=abs(effect)-1e-12
        assert np.isclose((extreme+1)/20001,row.p,atol=1e-14,rtol=0)
        for _ in range(2000):
            bx=rng.integers(0,n,n);b=float(np.sum(x[bx])/n) if y is None else corr(x[bx],y[bx])
            if np.isfinite(b):boots.append(b)
        assert len(boots)==row.valid_bootstrap
        assert np.allclose(np.quantile(boots,[.025,.975]),[row.lower,row.upper],rtol=1e-9,atol=1e-9)
    p=r.p.fillna(1).to_numpy();expected=[min(1,min(value*12/(p<=value).sum() for value in p if value>=own)) for own in p]
    valid=r.p.notna().to_numpy();assert np.allclose(np.asarray(expected)[valid],r.BH_twelve_slot_family.to_numpy()[valid])
    report=dict(status='PASS_ALL6_PAIRED_EFFECTS_PERMUTATIONS_BOOTSTRAPS_AND12_SLOT_BH',
      finished_UTC=datetime.now(timezone.utc).isoformat(),result_sha256=hashlib.sha256((out/'all12_private_paired_statistics.tsv').read_bytes()).hexdigest())
    (out/'private_independent_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));print(r.to_json(orient='records'))


if __name__=='__main__':main()
