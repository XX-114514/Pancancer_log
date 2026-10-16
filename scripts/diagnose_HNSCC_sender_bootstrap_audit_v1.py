"""Identify numerical identifiability disagreement; no new hypotheses or P."""
import hashlib,json,os
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from threadpoolctl import threadpool_limits
from test_HNSCC_frozen_sender_association_v1 import fitted
from audit_HNSCC_frozen_sender_association_v1 import fit


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=root/'HNSCC_sender_bootstrap_audit_diagnosis_20261016_v1';assert not out.exists();out.mkdir()
    result=root/'HNSCC_frozen_sender_association_results_20261016_v1'
    progress=result/'private_independent_audit_progress.json'
    checked=json.loads(progress.read_text())['checked_primary_models'] if progress.exists() else 0
    stats=pd.read_csv(result/'all400_private_primary_association_statistics.tsv',sep='\t');row=stats[stats.p.notna()].iloc[checked]
    prep=root/'HNSCC_frozen_sender_association_preflight_20261016_v1'
    config=json.loads((prep/'contract_before_any_association_effects.json').read_text());frames=pd.read_csv(prep/'private_all_planned_patient_designs.tsv',sep='\t')
    sub=frames[frames.ligand.eq(row.ligand)&frames.sender_group.eq(row.sender_group)&frames.included_in_full_design].sort_values('patient_id')
    fields=['exposure','response']+config['categorical_nuisance']+config['numeric_nuisance'];data=sub[fields].to_numpy(float)
    strata=data[:,2]+2*data[:,3];groups=[np.flatnonzero(strata==v) for v in sorted(set(strata))]
    rng=np.random.default_rng(20261016+int(hashlib.sha256((row.ligand+'|'+row.sender_group).encode()).hexdigest()[:8],16))
    for _ in range(80):
        for group in groups:
            for j in range(250):rng.permutation(group)
    counts=Counter();diagnostics=[]
    with threadpool_limits(limits=1):
        for iteration in range(2000):
            ix=np.concatenate([rng.choice(group,len(group),replace=True) for group in groups]);a=data[ix]
            try:fitted(a);p='valid'
            except (ValueError,np.linalg.LinAlgError) as e:p=str(e)
            try:independent=fit(a);q='valid'
            except (ValueError,np.linalg.LinAlgError) as e:q=str(e)
            counts[(p,q)]+=1
            if (p=='valid') != (q=='valid'):
                if q=='valid':
                    D=independent[2];x=D[:,-1];N=D[:,:-1]
                    residual=x-N@np.linalg.lstsq(N,x,rcond=None)[0]
                    diagnostics.append(dict(iteration=iteration,producer=p,auditor=q,
                        exposure_residual_squared_norm=float(residual@residual),condition=float(np.linalg.cond(D))))
    table=pd.DataFrame([dict(producer=p,auditor=q,draws=n) for (p,q),n in counts.items()])
    table.to_csv(out/'private_bootstrap_acceptance_crosscheck.tsv',sep='\t',index=False)
    pd.DataFrame(diagnostics).to_csv(out/'private_disagreement_numerical_diagnostics.tsv',sep='\t',index=False)
    report=dict(status='BOOTSTRAP_AUDIT_ACCEPTANCE_DISAGREEMENT_LOCALIZED_NO_NEW_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        previously_checked_primary_models=checked,ligand=row.ligand,sender_group=row.sender_group,
        primary_reported_valid=int(row.valid_bootstrap),acceptance_crosscheck=table.to_dict('records'),
        disagreements=len(diagnostics),maximum_disagreement_residual_norm=max((r['exposure_residual_squared_norm'] for r in diagnostics),default=None),
        primary_statistics_changed=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
