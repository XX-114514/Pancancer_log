"""All pre-effect eligible pair cases, all disjoint functions; no main P read."""
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_external_programme_coordination_math_v1 import prepare
from test_NSCLC_external_programme_coordination_v1 import read_inputs, QCOLS


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 32*1024**2)
    out = root/'NSCLC_external_disjoint_function_coupling_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), budget=budget,
                scope='All pre-effect eligible record slots/all50 functions; identical paired-library patient averages, active cases and full own/partner quality/category span.',
                effects='Own and cross-end programme/function partial ranks and function/function partial ranks; ULM, target mean and nontarget mean retained. No P/FDR or primary-family changes.',
                support='Pooled expressed-target counts/target means reported; low support cannot substantiate activation. Missing readouts remain unavailable, no case rescue.',
                limitations='Gene disjointness is not independent patients or modality; shared RNA denominators/sample biology can couple readouts. Pair-conditioned same-view functions differ with partner QC/cases; all rows retained.',
                main_coordination_results_or_clinical_values_loaded=False)
    (out/'plan_before_all_functional_coupling_descriptions.json').write_text(json.dumps(plan, indent=2)+'\n')
    functional = root/'NSCLC_external_disjoint_functional_RNA_20261016_v1'
    proof = functional/'private_independent_disjoint_functional_RNA_audit.json'
    h = root/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_external_disjoint_functional_RNA_v1'
    receipt = json.loads((h/'launch_receipt.json').read_text())
    deadline = time.monotonic()+3600
    while not proof.exists():
        p = Path('/proc')/str(receipt['pid'])/'stat'
        assert p.exists(), 'Bound functional audit missing without proof'
        f = p.read_text().rsplit(')', 1)[1].split()
        assert f[0] != 'Z' and f[19] == str(receipt['proc_start_ticks'])
        assert time.monotonic()<deadline, 'Observation expired; do not infer termination'
        time.sleep(3)
    audit = json.loads(proof.read_text())
    assert audit['status'].startswith('PASS_ALL61_UNIVERSES3050_GATES6760_')
    assert audit['producer_terminal_sha256'] == sha(functional/'private_terminal.json')
    source, _, _, slots, cases = read_inputs(root)
    src = root/'NSCLC_external_programme_coordination_20261016_v1'
    scope = pd.read_csv(src/'private_all15600_primary_test_eligibility.tsv', sep='\t')
    eligible = scope[scope.test_index.ge(0)]
    with np.load(functional/'private_all6760_by50_disjoint_functional_RNA.npz', allow_pickle=False) as z:
        names = z['functions'].astype(str).tolist()
        values = {k:z[k] for k in ['ULM_score','target_mean','nontarget_mean','expressed_target_genes']}
    members = pd.read_csv(source/'private_same_library_pair_membership.tsv.gz', sep='\t')
    rid = members.row_id.to_numpy(int)
    count = np.bincount(rid, minlength=len(cases))
    aggregate = {k:np.full((len(cases),2,50),np.nan) for k in values}
    for side, col in enumerate(['unit_A','unit_B']):
        units = members[col].to_numpy(int)
        for key in values:
            for j in range(50):
                x = values[key][units,j]
                number = np.bincount(rid,weights=np.isfinite(x).astype(int),minlength=len(cases))
                total = np.bincount(rid,weights=np.nan_to_num(x),minlength=len(cases))
                good = (count>0)&(number==count)
                aggregate[key][good,side,j] = total[good]/count[good]
    index = slots.set_index('slot')
    rows = []
    for item in eligible.itertuples():
        meta = index.loc[item.slot]
        g = cases.iloc[int(meta.row_start):int(meta.row_stop)]
        if item.family=='reviewed_patients':g=g[g.reviewed_patient]
        p = prepare(g.raw_A.to_numpy(),g.raw_B.to_numpy(),g[QCOLS].to_numpy(),g.block_key.to_numpy(str))
        g = g.iloc[np.flatnonzero(p['keep'])]
        assert len(g)==item.active_units
        ids=g.row_id.to_numpy(int);q=p['Q']
        fields=[dict(record=item.record,family=item.family,slot=item.slot,factor=meta.factor,
                     view_A=meta.view_A,view_B=meta.view_B,dataset=meta.dataset,cancer_label=meta.cancer_label,
                     material_class=meta.material_class,phase_key=meta.phase_key,reference=name,
                     active_units=len(g),reviewed_units=int(g.reviewed_patient.sum()),quality_rank=p['design_rank']) for name in names]
        for key in ['ULM_score','target_mean','nontarget_mean']:
            residual=[];norms=[]
            for side in range(2):
                x=aggregate[key][ids,side];good=np.isfinite(x).all(axis=0)
                r=np.full_like(x,np.nan);r[:,good]=rankdata(x[:,good],axis=0,method='average')
                r[:,good]-=q@(q.T@r[:,good]);n=np.linalg.norm(r,axis=0);n[n<1e-10]=np.nan
                residual.append(r);norms.append(n)
            for j in range(50):
                combinations=[('programmeA_functionA',p['u'],residual[0][:,j],p['norm_u']*norms[0][j]),
                              ('programmeB_functionB',p['v'],residual[1][:,j],p['norm_v']*norms[1][j]),
                              ('programmeA_functionB',p['u'],residual[1][:,j],p['norm_u']*norms[1][j]),
                              ('programmeB_functionA',p['v'],residual[0][:,j],p['norm_v']*norms[0][j]),
                              ('functionA_functionB',residual[0][:,j],residual[1][:,j],norms[0][j]*norms[1][j])]
                for tag,a,b,den in combinations:
                    fields[j][key+'_'+tag+'_partial_rho']=float(np.clip(a@b/den,-1,1)) if np.isfinite(den) else np.nan
                fields[j][key+'_same_cases_both_ends_available']=bool(np.isfinite(norms[0][j]) and np.isfinite(norms[1][j]))
        for j in range(50):
            for side,label in enumerate(['A','B']):
                positive=aggregate['expressed_target_genes'][ids,side,j];target=aggregate['target_mean'][ids,side,j]
                good=np.isfinite(positive).all() and np.isfinite(target).all()
                fields[j]['median_expressed_targets_'+label]=float(np.median(positive)) if good else np.nan
                fields[j]['median_target_mean_'+label]=float(np.median(target)) if good else np.nan
                fields[j]['low_foreground_support_'+label]=bool(good and (np.median(positive)<10 or np.median(target)<=0))
            rows.append(fields[j])
        if item.test_index%10==0:
            (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),completed_record_index=int(item.test_index),total_records=len(eligible)))+'\n')
    table=pd.DataFrame(rows)
    assert len(table)==len(eligible)*50 and not table.duplicated(['record','reference']).any()
    table.to_csv(out/'private_all_fixed_disjoint_function_couplings.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
    terminal=dict(status='COMPLETE_ALL_ELIGIBLE_RECORDS50_DISJOINT_FUNCTION_COUPLING_DESCRIPTIONS_AUDIT_PENDING_NO_P',
                  finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),functional_audit_sha256=sha(proof),
                  scope_sha256=sha(src/'private_all15600_primary_test_eligibility.tsv'),
                  files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
                  main_coordination_results_or_clinical_values_loaded=False,new_P=False)
    assert sum(v['bytes'] for v in terminal['files'].values())<budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
