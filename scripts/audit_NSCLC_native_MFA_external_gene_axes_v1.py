"""Independent augmented-QR projection and complete physical mapping replay."""
import json
from collections import Counter
from datetime import datetime,timezone

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from NSCLC_native_MFA_clinical_inputs_v1 import matched_quota
from external_functional_coverage_inputs_v1 import read_inputs


def main():
    data=load();root=data['root'];src=root/'NSCLC_native_MFA_external_gene_qualification_20261016_v1'
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    proof=src/'private_independent_gene_qualification_audit.json';assert not proof.exists();matched_quota(root)
    sources=read_inputs(root)['source_axes'];table=pd.read_csv(src/'private_all2440_native_factor_source_gene_coverage.tsv',sep='\t')
    assert len(table)==2440 and not table.duplicated(['arm','dataset','view','factor']).any()
    checked=0;max_error=0.
    for vi,view in enumerate(VIEWS):
        w=data['model']['weights_'+str(vi)];tau=data['model']['precision_'+str(vi)];sw=np.sqrt(tau)
        augmented=np.vstack([sw[:,None]*w,np.eye(10)]);q,r=np.linalg.qr(augmented,mode='reduced')
        coefficient=solve_triangular(r,q[:len(w)].T*sw)
        genes=data['genes'][data['prepared'][vi]['selected_native_columns']]
        for (arm,dataset),axis in sorted(sources.items()):
            chosen=[next(iter(axis['lookup'][g])) if g in axis['lookup'] and len(axis['lookup'][g])==1 else None for g in genes]
            count=Counter(v for v in chosen if v is not None)
            valid=np.array([v is not None and count[v]==1 for v in chosen])
            for k in range(10):
                factor='NSCLC_MFA_20261016_F%02d'%(k+1)
                row=table[(table.arm==arm)&(table.dataset==dataset)&(table['view']==view)&(table.factor==factor)]
                assert len(row)==1;row=row.iloc[0];v=coefficient[k];total=np.sum(abs(v));pos=v>0;neg=v<0
                a=np.sum(abs(v[valid]))/total if total>1e-12 else np.nan
                pm=np.sum(v[pos]);nm=np.sum(-v[neg]);pc=np.sum(v[valid&pos])/pm if pm>1e-12 else np.nan;nc=np.sum(-v[valid&neg])/nm if nm>1e-12 else np.nan
                observed=np.array([row.absolute_coefficient_coverage,row.positive_mass_coverage,row.negative_mass_coverage]);expected=np.array([a,pc,nc])
                assert np.allclose(observed,expected,atol=1e-9,rtol=1e-9,equal_nan=True)
                finite=np.isfinite(observed)&np.isfinite(expected)
                if finite.any():max_error=max(max_error,float(np.max(abs(observed[finite]-expected[finite]))))
                passed=bool(valid.sum()>=500 and a>=.95 and (pm<=1e-12 or pc>=.9) and (nm<=1e-12 or nc>=.9))
                assert int(row.uniquely_mapped_features)==valid.sum() and bool(row.physical_coverage_pass)==passed
                assert row.source_role==('DISCOVERY_NOT_EXTERNAL' if dataset=='GSE243013' else 'OVERLAP_AND_PATIENT_INDEPENDENCE_UNREVIEWED')
                checked+=1
    report=dict(status='PASS_ALL2440_SOURCE_FACTOR_VIEW_GENE_COVERAGE_INDEPENDENT_AUGMENTED_QR_MAPPING_REPLAY',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),all_records_checked=checked,maximum_coefficient_mass_coverage_error=max_error,
        expression_or_external_clinical_effects_computed=False,compatible_lineages_or_patients_accepted=False)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
