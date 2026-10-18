"""Independently accept all RNA adjustments and source consistency rows.

Scientific tables and original forest plots stay private. No row export.
"""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_NSCLC_external_programme_coordination_v2 import read_inputs,QCOLS
from audit_NSCLC_external_programme_coordination_v3 import calculate


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);resource(root,1048576)
    src=root/'NSCLC_external_programme_coordination_review_20261016_v3'
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    primary=root/'NSCLC_external_programme_coordination_20261016_v2'
    assert sha(primary/'private_independent_coordination_audit.json')==terminal['coordination_audit_sha256']
    source,_,_,slots,cases=read_inputs(root)
    members=pd.read_csv(source/'private_same_library_pair_membership.tsv.gz',sep='\t')
    diagnostic=root/'NSCLC_external_public_RNA_diagnostics_20261016_v1'
    assert sha(diagnostic/'private_independent_public_RNA_diagnostic_audit.json')==terminal['diagnostic_audit_sha256']
    values=np.load(diagnostic/'private_all6760_public_RNA_diagnostics.npz',allow_pickle=False)['values']
    aggregation=csr_matrix((np.ones(len(members)),(members.row_id,np.arange(len(members)))),shape=(len(cases),len(members)))
    counts=np.asarray(aggregation.sum(axis=1)).ravel();panels=np.full((len(cases),4),np.nan)
    for side,offset in [('A',0),('B',2)]:
        for col in range(2):
            x=values[members['unit_'+side].to_numpy(int),col]
            finite=np.asarray(aggregation@np.isfinite(x).astype(int)).ravel()
            summed=np.asarray(aggregation@np.nan_to_num(x)).ravel()
            good=(counts>0)&(finite==counts);panels[good,offset+col]=summed[good]/counts[good]
    original=pd.read_csv(primary/'private_all_estimated_programme_coordination.tsv',sep='\t')
    full=pd.read_csv(src/'private_all_primary_and_sensitivity_review.tsv',sep='\t')
    pd.testing.assert_frame_equal(original,full[original.columns],check_dtype=False)
    lookup=slots.set_index('slot');maximum=0.
    for row in full.itertuples():
        meta=lookup.loc[row.slot];g=cases.iloc[int(meta.row_start):int(meta.row_stop)]
        if row.family=='reviewed_patients':g=g[g.reviewed_patient]
        p=calculate(g.raw_A.to_numpy(),g.raw_B.to_numpy(),g[QCOLS].to_numpy(),g.block_key.to_numpy(str))
        frame=g.iloc[np.flatnonzero(p['keep'])];extra=panels[frame.row_id.to_numpy(int)]
        assert len(frame)==row.active_units
        for name,columns in [('RNA_sex',[0,2]),('IEG_stress',[1,3]),('both_panels',[0,1,2,3])]:
            q=np.column_stack([p['quality'],extra[:,columns]])
            try:rho=calculate(p['a'],p['b'],q,p['blocks'])['rho'];status='ESTIMABLE'
            except ValueError as error:rho=np.nan;status=str(error)
            observed=getattr(row,name+'_rho')
            assert np.isclose(rho,observed,atol=1e-9,rtol=1e-9,equal_nan=True)
            assert getattr(row,name+'_status')==status
            assert getattr(row,name+'_complete_case_units')==np.isfinite(extra[:,columns]).all(axis=1).sum()
            attenuated=bool(np.isfinite(rho) and (abs(rho)<.5*abs(row.partial_rank_rho) or rho*row.partial_rank_rho<0))
            assert getattr(row,name+'_attenuation_or_reversal')==attenuated
            if np.isfinite(rho):maximum=max(maximum,abs(rho-observed))
    keys=['family','factor','view_A','view_B','material_group','dataset']
    expected=full.sort_values(['phase_priority','active_units','reviewed_units','context_id','arm'],ascending=[True,False,False,True,True]).drop_duplicates(keys).copy()
    observed=pd.read_csv(src/'private_effect_blind_one_context_per_source.tsv',sep='\t')
    # Producer's extra phase-key tie-breaker has no effect in this frozen input.
    assert set(expected.record)==set(observed.record)
    pd.testing.assert_frame_equal(expected.set_index('record').sort_index(),observed.set_index('record').sort_index()[expected.columns.drop('record')],check_dtype=False)
    expected['support']=expected.BH_q.lt(.05)&((expected.CI_low.gt(0)&expected.partial_rank_rho.gt(0))|(expected.CI_high.lt(0)&expected.partial_rank_rho.lt(0)))
    assert np.array_equal(expected.set_index('record').support.sort_index(),observed.set_index('record').BH_CI_support.sort_index())
    summary=pd.read_csv(src/'private_all_programme_pair_material_source_consistency.tsv',sep='\t').set_index(keys[:-1])
    checked=0
    for levels,g in expected.groupby(keys[:-1],sort=True):
        actual=summary.loc[levels];pos=g[g.support&g.partial_rank_rho.gt(0)];neg=g[g.support&g.partial_rank_rho.lt(0)]
        if len(pos) and len(neg):grade='DIRECTION_HETEROGENEOUS'
        elif max(len(pos),len(neg))>=2:grade='CROSS_SOURCE_CONSISTENCY_CANDIDATE' if levels[0]=='reviewed_patients' else 'LIBRARY_LEVEL_CROSS_SOURCE_LEAD_PATIENT_REVIEW_REQUIRED'
        elif max(len(pos),len(neg))==1:grade='SOURCE_LIMITED_CANDIDATE' if levels[0]=='reviewed_patients' else 'LIBRARY_LEVEL_SOURCE_LEAD_PATIENT_REVIEW_REQUIRED'
        else:grade='INSUFFICIENT_STATISTICAL_SUPPORT'
        assert actual.grade==grade and actual.unique_sources==g.dataset.nunique()==len(g)
        assert actual.BH_CI_supported_positive_sources==len(pos) and actual.BH_CI_supported_negative_sources==len(neg)
        assert actual.positive_source_effects==g.partial_rank_rho.gt(0).sum() and actual.negative_source_effects==g.partial_rank_rho.lt(0).sum()
        assert actual.recorded_cancer_label_strings==g.cancer_label.nunique()
        weights=np.maximum(g.active_units.to_numpy()-g.quality_rank.to_numpy()-3,1)
        effects=np.arctanh(np.clip(g.partial_rank_rho.to_numpy(),-1+1e-12,1-1e-12))
        center=float(np.tanh(sum(weights*effects)/sum(weights)))
        leaves=[float(np.tanh((sum(weights*effects)-weights[j]*effects[j])/(sum(weights)-weights[j]))) for j in range(len(g))] if len(g)>1 else []
        want=[center,g.partial_rank_rho.min(),g.partial_rank_rho.max(),min(leaves) if leaves else np.nan,max(leaves) if leaves else np.nan]
        have=[actual.descriptive_weighted_rho,actual.minimum_source_rho,actual.maximum_source_rho,actual.leave_one_source_min,actual.leave_one_source_max]
        assert np.allclose(want,have,atol=1e-12,rtol=1e-12,equal_nan=True)
        checked+=1
    assert checked==len(summary) and len(terminal['figures'])==10
    report=dict(status='PASS_ALL630_RNA_BACKGROUND_ADJUSTMENTS_EFFECT_BLIND_SOURCE_SELECTION_AND198_CONSISTENCY_SUMMARIES',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),review_terminal_sha256=sha(src/'private_terminal.json'),primary_result_audit_sha256=terminal['coordination_audit_sha256'],RNA_adjustment_records_checked=len(full)*3,consistency_rows_checked=checked,maximum_OLS_difference=maximum,extra_phase_key_sort_has_no_effect_on_frozen_input=True,visual_review_still_required=True,clinical_or_direct_communication_validation=False)
    proof=src/'private_independent_context_review_audit.json';assert not proof.exists();proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
