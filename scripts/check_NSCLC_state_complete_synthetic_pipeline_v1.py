"""Full synthetic runtime exercise; same wild/bootstrap/R paths, no human data.

Outputs must be labeled synthetic and are never scientific source evidence.
This checks the governance path and all fitted branches before RAM handoff.
"""
from datetime import datetime,timezone
from pathlib import Path
import argparse
import copy
import hashlib
import itertools
import json
import os
import numpy as np
from test_NSCLC_observed_T_state_statistics_v1 import fixture
import NSCLC_observed_T_state_statistics_v1 as statistics


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    task=args.output;task.mkdir(exist_ok=False);(task/'native_measurement_v1').mkdir()
    rng=np.random.default_rng(42097);original=fixture();base={r['case_index']:r for r in original['rows'] if r['QC']==statistics.MAIN_QC and r['parent']==statistics.MAIN_PARENT}
    for i,r in base.items():
        response=r['response'];latent=.7*response+rng.normal(0,.4)
        R=50+i%20;E=100+i%30;O=120+i%37;T=R+E+O
        rf=np.clip(.3+.25*latent+rng.normal(0,.07,size=3),.03,.9)
        ef=np.clip(np.array([.3+.25*latent,.3+.3*latent,.7-.25*latent])+rng.normal(0,.07,size=3),.03,.95)
        r['metrics']=[T,R,E,0,O]+np.r_[np.round(rf*R),np.round(ef*E)].astype(int).tolist()+[0,0]
        r['clinical_eligible']=i<80;r['cancer_type']='LUAD' if i<40 else 'LUSC';r['clinical'][3]=float(r['cancer_type']=='LUSC')
        r['parent_counts']=dict(R=R,E=E,O=O)
        for name,n in [('R',R),('E',E),('O',O)]:r['quality'][name]=[rng.uniform(6,8),rng.uniform(5,6),float(np.log1p(n)),n/T]
        r['other_T_states']=np.round(np.clip(rng.normal(.3,.1,size=3),.01,.9)*O).astype(int).tolist()
    payload=dict(rows=[])
    for qc in ['native_QC500_250_MT20','all_retained_native']:
        for parent in ['dual_regulatory_evidence','costate_regulatory_sensitivity']:
            for i in sorted(base):
                row=copy.deepcopy(base[i]);row.update(QC=qc,parent=parent);payload['rows'].append(row)
    code=Path(__file__).parent
    names=['NSCLC_observed_T_state_statistics_v1.py','GSE288199_state_pair_math_v1.py','GSE288199_state_pair_R_audit_v1.py','observed_T_marker_evidence_v1.py']
    # Synthetic Cartesian exercise, not the privately selected scientific pairs.
    pairs=[list(v)+[None] for v in itertools.product(statistics.STATES[:3],statistics.STATES[3:])]
    plan=dict(synthetic_only=True,created_UTC=datetime.now(timezone.utc).isoformat(),
        Rscript_path=os.environ['PANCANCER_R_BLOCK_REFERENCE'],fixed_pairs=pairs,
        code_SHA256={n:hashlib.sha256((code/n).read_bytes()).hexdigest() for n in names})
    (task/'private_frozen_NSCLC_state_statistics_v1.json').write_text(json.dumps(plan,indent=2)+'\n')
    statistics.run(payload,task,plan)
    out=task/'state_statistics_v1';read=lambda n:json.loads((out/n).read_text())
    terminal=read('private_terminal.json');assert terminal['status'].startswith('COMPLETE_')
    assert len(read('private_all6_clinical_state_associations.json'))==6
    assert len(read('private_pair_state_coordination.json'))==9 and len(read('private_point_sensitivities.json'))==129
    audit=read('private_all_real_design_R_random_stream_audit.json');assert len(audit)==33
    assert sum(r['full_wild_statistics_checked']==19999 for r in audit)==15
    assert all(r['full_bootstrap_statistics_checked']==2000 for r in audit)
    result=dict(status='PASS_FULL_SYNTHETIC_ONLY_PIPELINE_ALL33_FITTED_BRANCHES',
        completed_UTC=datetime.now(timezone.utc).isoformat(),actual_human_data_used=False,
        entire_primary_wild_and_all_bootstrap_R_replayed=True,working_null_gate_status=terminal['status'],
        no_limited_resampling_shortcut=True)
    (task/'synthetic_complete_runtime_acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
