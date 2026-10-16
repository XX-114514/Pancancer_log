"""Diagnose first unaccepted bootstrap disagreement; no result selection."""
import json
import os
from datetime import datetime,timezone
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

import NSCLC_external_programme_coordination_math_v1 as production
from audit_NSCLC_external_programme_coordination_v1 import calculate
from test_NSCLC_external_programme_coordination_v1 import read_inputs,QCOLS
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def exact_rank(matrix):
    a=[[Fraction(int(v)) for v in row] for row in matrix]
    rank=0
    for col in range(matrix.shape[1]):
        pivot=next((j for j in range(rank,len(a)) if a[j][col]),None)
        if pivot is None:continue
        a[rank],a[pivot]=a[pivot],a[rank]
        div=a[rank][col];a[rank]=[v/div for v in a[rank]]
        for j in range(rank+1,len(a)):
            factor=a[j][col]
            if factor:a[j]=[u-factor*v for u,v in zip(a[j],a[rank])]
        rank+=1
        if rank==len(a):break
    return rank


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);resource(root,2*1024**2)
    h=root/'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_external_programme_coordination_v1'
    receipt=json.loads((h/'launch_receipt.json').read_text())
    assert not (Path('/proc')/str(receipt['pid'])/'stat').exists(),'Audit still live, do not diagnose as stopped'
    assert 'line 144' in (h/'private_job_stdout_stderr.log').read_text()
    last=json.loads((h/'private_audit_checkpoint.json').read_text())['completed_tests']
    out=root/'NSCLC_external_bootstrap_numeric_diagnosis_20261016_v2';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_numeric_diagnosis.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),first_unaccepted_after_checkpoint=last,scope='First pending audit records only, inspect numerical mismatch and exact rational span. No main P/source-effect export or clinical selection.'))+'\n')
    _,_,_,slots,cases=read_inputs(root);index=slots.set_index('slot')
    src=root/'NSCLC_external_programme_coordination_20261016_v1'
    scope=pd.read_csv(src/'private_all15600_primary_test_eligibility.tsv',sep='\t');scope=scope[scope.test_index.ge(0)].set_index('test_index')
    saved=np.load(src/'private_all_bootstrap_and_LOO.npz',allow_pickle=False)['bootstrap']
    diagnosis=[]
    for i in range(last,min(last+10,len(scope))):
        item=scope.loc[i];meta=index.loc[item.slot];g=cases.iloc[int(meta.row_start):int(meta.row_stop)]
        if item.family=='reviewed_patients':g=g[g.reviewed_patient]
        p=production.prepare(g.raw_A.to_numpy(),g.raw_B.to_numpy(),g[QCOLS].to_numpy(),g.block_key.to_numpy(str))
        rng=np.random.default_rng(production.seed(item.record,'bootstrap'))
        groups=[np.flatnonzero(p['blocks']==s) for s in sorted(set(p['blocks']))]
        replay=np.full(2000,np.nan);independent=replay.copy();examples=[]
        for b in range(2000):
            take=np.concatenate([ids[rng.integers(len(ids),size=len(ids))] for ids in groups])
            a,bb,q,blocks=p['a'][take],p['b'][take],p['quality'][take],p['blocks'][take]
            try:
                pp=production.prepare(a,bb,q,blocks);replay[b]=production.effect(pp)
            except ValueError:pp=None
            try:
                ii=calculate(a,bb,q,blocks);independent[b]=ii['rho']
            except ValueError:ii=None
            agree=np.isclose(replay[b],independent[b],equal_nan=True,atol=1e-9,rtol=1e-9)
            if not agree and len(examples)<3:
                qr=rankdata(q,axis=0)*2;ar=rankdata(a)*2;br=rankdata(bb)*2
                z=np.column_stack([np.ones(len(a)),*[(blocks==s).astype(int) for s in sorted(set(blocks))],qr]).astype(int)
                er=exact_rank(z);era=exact_rank(np.column_stack([z,ar.astype(int)]));erb=exact_rank(np.column_stack([z,br.astype(int)]))
                info=dict(draw=b,draw_unique_input_units=int(np.unique(take).size),exact_design_rank=er,
                          exact_rank_with_A=era,exact_rank_with_B=erb,
                          production_valid=pp is not None,independent_valid=ii is not None,
                          saved_value=float(saved[i,b]) if np.isfinite(saved[i,b]) else None,
                          production_replay=float(replay[b]) if np.isfinite(replay[b]) else None,
                          independent_replay=float(independent[b]) if np.isfinite(independent[b]) else None)
                if pp is not None:
                    singular=np.linalg.svd(pp['design'],compute_uv=False);positive=singular[:pp['design_rank']]
                    info.update(production_rank=pp['design_rank'],production_norm_A=pp['norm_u'],production_norm_B=pp['norm_v'],retained_condition_number=float(positive[0]/positive[-1]))
                if ii is not None:info.update(independent_rank=int(ii['rank']),independent_residual_norms=np.linalg.norm(ii['residual'],axis=0).tolist())
                examples.append(info)
        assert np.allclose(replay,saved[i],equal_nan=True,atol=1e-12,rtol=1e-12),'Production no longer replays original saved stream'
        valid=np.isfinite(replay)&np.isfinite(independent)
        record=dict(test_index=i,production_valid_draws=int(np.isfinite(replay).sum()),independent_valid_draws=int(np.isfinite(independent).sum()),
                    missing_mask_disagreements=int((np.isfinite(replay)!=np.isfinite(independent)).sum()),
                    maximum_common_finite_difference=float(np.max(abs(replay[valid]-independent[valid]))) if valid.any() else None,
                    examples=examples)
        diagnosis.append(record)
        if examples:break
    proof=out/'private_first_bootstrap_disagreement_diagnosis.json'
    report=dict(status='COMPLETE_ORIGINAL_BOOTSTRAP_REPLAY_AND_EXACT_SPAN_DIAGNOSIS_NO_PRODUCTION_CHANGE',
                finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),diagnosis=diagnosis,
                primary_scope_or_thresholds_changed=False,source_identity_patient_rows_or_gene_lists_exported=False)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    public=dict(status=report['status'],finished_UTC=report['finished_UTC'],diagnosis_sha256=sha(proof),
                records_checked=len(diagnosis),disagreeing_draw_examples=sum(len(r['examples']) for r in diagnosis),
                exact_span_A_constant_examples=sum(e['exact_rank_with_A']==e['exact_design_rank'] for r in diagnosis for e in r['examples']),
                exact_span_B_constant_examples=sum(e['exact_rank_with_B']==e['exact_design_rank'] for r in diagnosis for e in r['examples']))
    print(json.dumps(public))


if __name__=='__main__':main()
