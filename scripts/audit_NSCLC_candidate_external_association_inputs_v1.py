"""Complete new paired-unit reconstruction and frozen source-model eligibility.

Reuses the independently accepted rank kernel for its design/zero guards.
No external ligand-response effects, P, resampling or clinical outcomes.
"""
import json
import os
import time
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha,VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_candidate_external_associations_v1 import OUT,SCOPE_SHA,QCOLS
from prepare_NSCLC_external_programme_pair_tests_v1 import phase_key,packed,CONTEXT,BLOCK
from NSCLC_external_lineage_inputs_v2 import load as source_metadata
from preflight_NSCLC_continuous_prior_native_association_v1 import design_gate
import NSCLC_external_programme_coordination_math_v2 as kernel


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);resource(root,1024**2)
    launch=json.loads((root/'HNSCC_receiver_followup_job_control_20261016_v2/prepare_NSCLC_candidate_external_associations_v1/launch_receipt.json').read_text())
    assert launch['source_sha256']=='b7944d3f6eb03057977f5debc289028c5e19309733e7a3423f1b9972d891a66c'
    out=root/OUT;deadline=time.monotonic()+3600
    while not (out/'private_terminal.json').exists():
        p=Path('/proc')/str(launch['pid'])/'stat';assert p.exists(),'Bound input producer absent without terminal'
        f=p.read_text().rsplit(')',1)[1].split();assert f[0]!='Z' and f[19]==str(launch['proc_start_ticks'])
        assert time.monotonic()<deadline,'Observation expired; no restart inferred';time.sleep(3)
    terminal=json.loads((out/'private_terminal.json').read_text());assert terminal['source_sha256']==launch['source_sha256']
    for name,item in terminal['files'].items():assert sha(out/name)==item['sha256']
    scope_path=root/'NSCLC_candidate_external_association_scope_20261018_v1.json';assert sha(scope_path)==SCOPE_SHA
    scope=json.loads(scope_path.read_text());candidates=pd.DataFrame(scope['candidates'])
    e=pd.read_csv(out/'private_all_planned_model_eligibility.tsv.gz',sep='\t',keep_default_na=False)
    cases=pd.read_csv(out/'private_all_case_measurement_units.tsv.gz',sep='\t',keep_default_na=False)
    library=pd.read_csv(out/'private_all_library_measurement_eligibility.tsv.gz',sep='\t',keep_default_na=False)
    excluded=pd.read_csv(out/'private_excluded_mixed_case_blocks.tsv.gz',sep='\t',keep_default_na=False)
    assert len(e)==scope['planned_positions'] and e.model.is_unique
    assert cases.case_row.tolist()==list(range(len(cases))) and not cases.duplicated(['record','context_id','resampling_block']).any()
    assert len(library)==len(candidates)*1690 and not library.duplicated(['record','library_index']).any()
    raw=root/'NSCLC_candidate_external_exact_RNA_20261018_v1'
    assert sha(raw/'private_independent_exact_panel_audit.json')==scope['raw_RNA_audit_sha256']
    with np.load(raw/'private_all_frozen_pool_exact_panel_measurements.npz',allow_pickle=False) as f:z={k:f[k] for k in f.files}
    masked=root/'NSCLC_candidate_external_original_masked_scores_20261018_v1';assert sha(masked/'private_terminal.json')==scope['masked_score_terminal_sha256']
    with np.load(masked/'private_all_candidate_library_masked_scores.npz',allow_pickle=False) as f:response={k:f[k] for k in f.files}
    assay=root/'NSCLC_external_fixed_gene_programme_assay_20261016_v1'
    units=pd.read_csv(assay/'private_all6760_programme_units_metadata_quality_TMM.tsv.gz',sep='\t',keep_default_na=False)
    with np.load(assay/'private_all6760_fixed_programme_indices.npz',allow_pickle=False) as f:quality=f['quality']
    lib=units.iloc[::4].copy();lib.index=np.arange(1690)
    metadata=source_metadata();assert metadata['lib'].library_id.tolist()==lib.library_id.tolist()
    denominator=metadata['lib'].annotated_scope_cells.to_numpy(float)
    fraction=pd.to_numeric(units.diagnostic_fraction_within_parent,errors='coerce').to_numpy(float)
    lib['phase_key']=[phase_key(r) for r in lib.itertuples()]
    lib['block_key']=[packed([getattr(r,c) for c in BLOCK]) for r in lib.itertuples()]
    import hashlib
    lib['context_id']=['PC_'+hashlib.sha256(packed([str(getattr(r,c)) for c in CONTEXT]).encode()).hexdigest()[:16] for r in lib.itertuples()]
    curated=json.loads((root/'NSCLC_native_interface_mechanistic_scope_20261018_v1.json').read_text())
    for field in ['ligand','response','receiver_own_ligand','fraction_within_parent_sender','fraction_within_parent_receiver','fraction_whole_sender','fraction_whole_receiver','sender_cells_sum','receiver_cells_sum','sender_ligand_positive']+QCOLS:
        cases[field]=pd.to_numeric(cases[field],errors='raise')
    max_error=0.;checked_cases=set();expected_mixed=set();counter=0;checked_models=set()
    for ri,item in enumerate(candidates.itertuples()):
        sv,rv=VIEWS.index(item.sender_view),VIEWS.index(item.receiver_view);gi=list(z['genes']).index(item.ligand)
        recs=['_'.join(v) for v in curated['ligand_specs'][item.ligand]['complexes']];rcols=[list(z['complexes']).index(v) for v in recs]
        for ctx in scope['contexts']:
            members=lib[lib.context_id.eq(ctx['context_id'])];ids=members.index.to_numpy(int);su=ids*4+sv;ru=ids*4+rv
            flags=dict(sender_RNA_measurement_available=units.iloc[su].RNA_measurement_available.to_numpy(bool),receiver_RNA_measurement_available=units.iloc[ru].RNA_measurement_available.to_numpy(bool),
                sender_cells_ge20=z['n_cells'][su]>=20,receiver_cells_ge20=z['n_cells'][ru]>=20,literal_ligand_measured=z['measured'][su,gi],
                receiver_masked_score_available=response['measurement_available'][ri,ids],receiver_whole_source_gate=response['original_full_source_gene_gate'][ri,ids],
                receiver_disjoint_source_gate=response['original_disjoint_source_gate'][ri,ids],any_complete_canonical_receptor_measured=z['complex_measured'][ru][:,rcols].any(axis=1),
                both_full_QC_finite=np.isfinite(quality[su]).all(axis=1)&np.isfinite(quality[ru]).all(axis=1),both_within_parent_fractions_finite=np.isfinite(fraction[su])&np.isfinite(fraction[ru]))
            available=np.logical_and.reduce(list(flags.values()))
            observed=library[library.record.eq(item.record)&library.context_id.eq(ctx['context_id'])].set_index('library_index').loc[ids]
            assert np.array_equal(observed.measurement_and_receiver_gates_pass,available) and not observed.sender_whole_programme_gate_required.any()
            for key,value in flags.items():assert np.array_equal(observed[key],value)
            for pos in range(len(ids)):assert observed.failed_measurement_flags.iloc[pos]=='|'.join(k for k,v in flags.items() if not v[pos])
            actual=cases[cases.record.eq(item.record)&cases.context_id.eq(ctx['context_id'])];expected_ids=[]
            for patient,g in members.iloc[np.flatnonzero(available)].groupby('resampling_block',sort=True):
                if g.block_key.nunique()!=1:
                    expected_mixed.add((item.record,ctx['context_id'],str(patient)));continue
                li=g.index.to_numpy(int);a,b=li*4+sv,li*4+rv
                row=actual[actual.resampling_block.eq(str(patient))];assert len(row)==1;row=row.iloc[0]
                expected_ids.append(int(row.case_row));checked_cases.add(int(row.case_row))
                assert json.loads(row.library_indices_json)==li.tolist() and row.paired_libraries==len(li)
                assert row.block_key==g.block_key.iloc[0] and bool(row.reviewed_patient)==bool(g.reviewed_patient.iloc[0])
                fields=dict(ligand=np.mean(z['mean_logCP10K'][a,gi]),response=item.direction*np.mean(response['response'][ri,li]),
                    receiver_own_ligand=np.mean(z['mean_logCP10K'][b,gi]),fraction_within_parent_sender=np.mean(fraction[a]),fraction_within_parent_receiver=np.mean(fraction[b]),
                    fraction_whole_sender=np.mean(z['n_cells'][a]/denominator[li]),fraction_whole_receiver=np.mean(z['n_cells'][b]/denominator[li]),
                    sender_cells_sum=np.sum(z['n_cells'][a]),receiver_cells_sum=np.sum(z['n_cells'][b]),sender_ligand_positive=np.sum(z['positive_counts'][a,gi]),
                    **{f'Q_sender_{j}':np.mean(quality[a,j]) for j in range(5)},**{f'Q_receiver_{j}':np.mean(quality[b,j]) for j in range(5)})
                for key,value in fields.items():
                    error=abs(float(row[key])-float(value));assert error<1e-9;max_error=max(max_error,error)
                for j,col in enumerate(rcols):
                    assert str(row[f'receptor_{j}_measured'])==str(bool(z['complex_measured'][b,col].all()))
                    for label,src in [('positive','receptor_positive'),('positive2','receptor_positive2')]:
                        value=z[src][b,col].sum();actual_value=pd.to_numeric(pd.Series([row[f'receptor_{j}_{label}']]),errors='coerce').iloc[0]
                        assert (pd.isna(actual_value) and np.isnan(value)) or actual_value==value
            assert set(expected_ids)==set(actual.case_row)
            for family in scope['families']:
                g=actual if family=='all_units' else actual[actual.reviewed_patient]
                candidate_count=len(g)
                if len(g):g=g[g.block_key.map(g.block_key.value_counts()).ge(2)]
                model=item.record+'|'+ctx['context_id']+'|'+family
                row=e[e.model.eq(model)];assert len(row)==1;row=row.iloc[0];checked_models.add(model)
                assert row.candidate_units==candidate_count and row.active_units==len(g)
                assert row.active_reviewed_units==int(g.reviewed_patient.sum())
                assert json.loads(row.active_case_rows_json)==g.case_row.tolist() and json.loads(row.canonical_receptors_json)==recs
                reasons=[];qualified=[];rank=df=0
                if ctx['in_model_source']:reasons.append('IN_MODEL_SOURCE_NOT_EXTERNAL')
                if len(g)<20:reasons.append('FEWER_THAN20_INFORMATIVE_UNITS')
                if len(g):
                    if g.sender_ligand_positive.sum()*100<g.sender_cells_sum.sum() or g.sender_ligand_positive.gt(0).sum()<10:reasons.append('SENDER_POOLED1PCT_OR10_UNIT_RNA_GATE_FAILED')
                    for j,rec in enumerate(recs):
                        measured=g[f'receptor_{j}_measured'].astype(str).eq('True').all()
                        pos=pd.to_numeric(g[f'receptor_{j}_positive'],errors='coerce')
                        if measured and pos.sum()*100>=g.receiver_cells_sum.sum() and pos.gt(0).sum()>=10:qualified.append(rec)
                    if not qualified:reasons.append('NO_COMPLETE_RECEPTOR_POOLED1PCT_AND10_UNIT_SUPPORT')
                    try:
                        p=kernel.prepare(g.ligand.to_numpy(float),g.response.to_numpy(float),g[QCOLS].to_numpy(float),g.block_key.to_numpy(str),min_n=20,min_df=8)
                        design_gate(p);assert p['keep'].all();rank=p['design_rank'];df=p['residual_df']
                    except ValueError as err:reasons.append('DESIGN:'+str(err))
                assert json.loads(row.RNA_qualified_receptors_json)==qualified and row.nuisance_rank==rank and row.residual_df==df
                assert row.status==(';'.join(reasons) if reasons else 'ELIGIBLE_FIXED_EXTERNAL_PRIMARY')
                assert row.test_index==(-1 if reasons else counter)
                if not reasons:counter+=1
    assert checked_cases==set(cases.case_row) and checked_models==set(e.model)
    assert expected_mixed==set(zip(excluded.record,excluded.context_id,excluded.resampling_block))
    frozen=json.loads((out/'private_frozen_families_before_effects.json').read_text())
    assert frozen['eligible_family_sizes']=={f:int((e.family.eq(f)&e.test_index.ge(0)).sum()) for f in scope['families']}
    report=dict(status='PASS_COMPLETE_NEW_LIBRARY_PAIRS_CASE_AGGREGATES_AND_ALL_FROZEN_MODEL_ELIGIBILITY',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),scope_sha256=SCOPE_SHA,
        producer_terminal_sha256=sha(out/'private_terminal.json'),all_planned_models_checked=len(e),all_case_rows_checked=len(cases),
        all_library_measurement_flags_checked=len(library),maximum_case_aggregate_error=max_error,
        limits='All pair/case numeric aggregation and eligibility replayed; same previously accepted rank kernel reused for design/zero guards. No inference/random streams or biological validation.',new_effects_P_or_clinical_tests=False)
    p=out/'private_independent_external_association_input_audit.json';assert not p.exists();p.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(status=report['status'],completed_UTC=report['completed_UTC'],proof_sha256=sha(p))))


if __name__=='__main__':main()
