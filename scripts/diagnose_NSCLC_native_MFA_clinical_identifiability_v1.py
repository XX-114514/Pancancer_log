"""Diagnose nuisance-span aliasing without any factor-outcome association."""
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.linalg import null_space

from NSCLC_native_MFA_interpretation_inputs_v1 import load, quota_gate, sha, VIEWS
from preflight_NSCLC_native_MFA_clinical_coverage_v1 import COMMON


def main():
    data=load();root=data['root'];pre=root/'NSCLC_native_MFA_clinical_coverage_20261016_v1'
    terminal=json.loads((pre/'private_terminal.json').read_text())
    path=pre/'private_all208_clinical_eligibility.tsv'
    assert sha(path)==terminal['files'][path.name]['sha256']
    out=root/'NSCLC_native_MFA_clinical_identifiability_20261016_v1';assert not out.exists();quota_gate(root);out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        source_clinical_terminal_sha256=sha(pre/'private_terminal.json'),
        question='Locate all fixed-control aliases and test whether the response contrast itself lies outside the nuisance span. No factor score enters this diagnostic; no clinical association coefficient or P.',
        scope='All4 views and all3 predeclared partitions, all original clinical and own-view quality columns retained.',
        decomposition='SVD relative tolerance1e-12; report rank gaps, null relations, response residual norm and leverage. A nuisance alias is not automatically non-identification of the clinical contrast.',
        followup='Any later equivalent nuisance-span parametrization requires an explicit pre-effect contract and independent input/calibration checks. No dropping columns or retrospective modification of the prior full-rank gates.')
    (out/'plan_before_nuisance_identifiability_diagnostic.json').write_text(json.dumps(plan,indent=2)+'\n')
    c=pd.read_csv(path,sep='\t').set_index('sampleID').loc[data['cases'].sampleID]
    results=[]
    for part in ['discovery','internal_holdout','all_frozen_four_view_cases']:
        original=np.ones(len(c),bool) if part=='all_frozen_four_view_cases' else c.partition.eq(part).to_numpy()
        take=original&c.clinical_eligible.to_numpy(bool);x=c.loc[take,'response'].to_numpy(float);n=len(x)
        for vi,view in enumerate(VIEWS):
            cov=np.column_stack([c.loc[take,COMMON].to_numpy(float),data['prepared'][vi]['quality'][take]])
            names=['intercept','response']+COMMON+data['prepared'][vi]['covariate_names'].astype(str).tolist()
            sd=cov.std(axis=0);cov=(cov-cov.mean(axis=0))/np.where(sd>1e-12,sd,1)
            nuisance=np.column_stack([np.ones(n),cov]);u,s,_=np.linalg.svd(nuisance,full_matrices=False)
            rank=int((s>s[0]*1e-12).sum());basis=u[:,:rank];res=x-basis@(basis.T@x)
            full=np.column_stack([np.ones(n),x,cov]);null=null_space(full,rcond=1e-12)
            relationships=[]
            for j in range(null.shape[1]):
                v=null[:,j];v*=1 if v[np.argmax(abs(v))]>=0 else -1
                relationships.append({name:float(value) for name,value in zip(names,v) if abs(value)>.01})
            norm=float(np.linalg.norm(res));identified=norm>1e-8*np.sqrt(n)
            q=np.column_stack([basis,res/norm]) if identified else basis
            results.append(dict(partition=part,view=view,cases=n,original_nuisance_columns=nuisance.shape[1],nuisance_rank=rank,
                retained_relative_singular_min=float(s[rank-1]/s[0]),discarded_relative_singular_max=float(s[rank]/s[0]) if rank<len(s) else None,
                full_design_null_relations=relationships,response_residual_norm=norm,response_contrast_identified=identified,
                equivalent_span_design_rank=q.shape[1],residual_df=n-q.shape[1],maximum_leverage=float(np.max(np.sum(q*q,axis=1))),
                factor_score_or_clinical_effect_computed=False))
    (out/'private_all12_nuisance_span_diagnostics.json').write_text(json.dumps(results,indent=2)+'\n')
    report=dict(status='COMPLETE_ALL12_NUISANCE_SPAN_AND_RESPONSE_IDENTIFIABILITY_DIAGNOSTICS_NO_FACTOR_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        any_original_control_deleted=False,clinical_tests_executed=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
