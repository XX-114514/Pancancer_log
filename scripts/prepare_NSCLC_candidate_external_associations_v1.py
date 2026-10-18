"""Freeze complete source-model eligibility before external ligand effects/P.

Rebuild paired libraries; do not inherit an unused sender-programme gene gate.
"""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha,VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_external_programme_pair_tests_v1 import phase_key,packed,CONTEXT,BLOCK
from NSCLC_external_lineage_inputs_v2 import load as source_metadata
from preflight_NSCLC_continuous_prior_native_association_v1 import design_gate
import NSCLC_external_programme_coordination_math_v2 as kernel

SCOPE_SHA='d4141c94fe11d2fbccd764fbc2bff503d8a068e7e32132cb6f5c790bb158d752'
OUT='NSCLC_candidate_external_association_inputs_20261018_v1'
QCOLS=[f'Q_{side}_{j}' for side in ['sender','receiver'] for j in range(5)]+['fraction_within_parent_sender','fraction_within_parent_receiver']


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,32*1024**2)
    path=root/'NSCLC_candidate_external_association_scope_20261018_v1.json';assert sha(path)==SCOPE_SHA
    scope=json.loads(path.read_text());candidates=pd.DataFrame(scope['candidates'])
    raw=root/'NSCLC_candidate_external_exact_RNA_20261018_v1'
    proof=raw/'private_independent_exact_panel_audit.json';assert sha(proof)==scope['raw_RNA_audit_sha256']
    assert json.loads(proof.read_text())['producer_terminal_sha256']==sha(raw/'private_terminal.json')
    manifest=json.loads((raw/'private_terminal.json').read_text())
    for name in ['private_all_frozen_pool_exact_panel_measurements.npz','private_all6760_exact_panel_unit_axis.tsv.gz']:
        assert sha(raw/name)==manifest['files'][name]['sha256']
    with np.load(raw/'private_all_frozen_pool_exact_panel_measurements.npz',allow_pickle=False) as f:z={k:f[k] for k in f.files}
    masked=root/'NSCLC_candidate_external_original_masked_scores_20261018_v1'
    assert sha(masked/'private_terminal.json')==scope['masked_score_terminal_sha256']
    manifest=json.loads((masked/'private_terminal.json').read_text())
    for name in ['private_all_candidate_library_masked_scores.npz','private_fixed_candidate_axis.tsv']:
        assert sha(masked/name)==manifest['files'][name]['sha256']
    with np.load(masked/'private_all_candidate_library_masked_scores.npz',allow_pickle=False) as f:response={k:f[k] for k in f.files}
    assert response['records'].tolist()==candidates.record.tolist()
    assay=root/'NSCLC_external_fixed_gene_programme_assay_20261016_v1'
    proof=assay/'private_independent_fixed_programme_assay_audit.json'
    assert sha(proof)=='d745f978950b0882f5ff30f3c975cd6c1e626296b2a55e4ab4464b081d44d936'
    assert json.loads(proof.read_text())['producer_terminal_sha256']==sha(assay/'private_terminal.json')
    manifest=json.loads((assay/'private_terminal.json').read_text())
    for name in ['private_all6760_programme_units_metadata_quality_TMM.tsv.gz','private_all6760_fixed_programme_indices.npz']:
        assert sha(assay/name)==manifest['files'][name]['sha256']
    units=pd.read_csv(assay/'private_all6760_programme_units_metadata_quality_TMM.tsv.gz',sep='\t',keep_default_na=False)
    with np.load(assay/'private_all6760_fixed_programme_indices.npz',allow_pickle=False) as f:quality=f['quality']
    assert np.array_equal(units.n_cells,z['n_cells']) and units.unit_index.tolist()==list(range(6760))
    lib=units.iloc[::4].copy();lib.index=np.arange(1690)
    metadata=source_metadata()
    assert metadata['hashes']==json.loads((raw/'private_terminal.json').read_text())['frozen_inputs']['metadata']
    assert metadata['lib'].library_id.tolist()==lib.library_id.tolist()
    assert metadata['lib'].dataset.tolist()==lib.dataset.tolist() and metadata['lib'].arm.tolist()==lib.arm.tolist()
    lib['phase_key']=[phase_key(x) for x in lib.itertuples()]
    lib['block_key']=[packed([getattr(x,c) for c in BLOCK]) for x in lib.itertuples()]
    lib['context_id']=['PC_'+hashlib.sha256(packed([str(getattr(x,c)) for c in CONTEXT]).encode()).hexdigest()[:16] for x in lib.itertuples()]
    assert set(lib.context_id)=={c['context_id'] for c in scope['contexts']}
    assert pd.api.types.is_bool_dtype(lib.reviewed_patient)
    fraction=pd.to_numeric(units.diagnostic_fraction_within_parent,errors='coerce').to_numpy(float)
    denominator=pd.to_numeric(metadata['lib'].annotated_scope_cells,errors='raise').to_numpy(float)
    assert (denominator>0).all()
    out=root/OUT;assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),scope_sha256=SCOPE_SHA,
        budget=budget,kernel_sha256=sha(kernel.__file__),QCOLS=QCOLS,
        frozen_source_scope=scope,old_sender_programme_gene_gate_inherited=False,
        external_ligand_programme_effects_or_P_computed=False)
    (out/'plan_before_external_association_qualification.json').write_text(json.dumps(plan,indent=2)+'\n')
    curated=json.loads((root/'NSCLC_native_interface_mechanistic_scope_20261018_v1.json').read_text())
    assert sha(root/'NSCLC_native_interface_mechanistic_scope_20261018_v1.json')=='b87fa57fb8246056d36c619c0edaa86dcd0b0d3d830d670524c1d8c6d3144550'
    genes=list(z['genes']);complexes=list(z['complexes']);records=[];cases=[];libraries=[];mixed=[];counter=0
    for ri,item in enumerate(candidates.itertuples()):
        sv,rv=VIEWS.index(item.sender_view),VIEWS.index(item.receiver_view);gi=genes.index(item.ligand)
        recs=['_'.join(x) for x in curated['ligand_specs'][item.ligand]['complexes']];rcols=[complexes.index(x) for x in recs]
        for ctx in scope['contexts']:
            members=lib[lib.context_id.eq(ctx['context_id'])];ids=members.index.to_numpy(int);su=ids*4+sv;ru=ids*4+rv
            flags=dict(sender_RNA_measurement_available=units.iloc[su].RNA_measurement_available.to_numpy(bool),
                receiver_RNA_measurement_available=units.iloc[ru].RNA_measurement_available.to_numpy(bool),
                sender_cells_ge20=z['n_cells'][su]>=20,receiver_cells_ge20=z['n_cells'][ru]>=20,
                literal_ligand_measured=z['measured'][su,gi],receiver_masked_score_available=response['measurement_available'][ri,ids],
                receiver_whole_source_gate=response['original_full_source_gene_gate'][ri,ids],
                receiver_disjoint_source_gate=response['original_disjoint_source_gate'][ri,ids],
                any_complete_canonical_receptor_measured=z['complex_measured'][ru][:,rcols].any(axis=1),
                both_full_QC_finite=np.isfinite(quality[su]).all(axis=1)&np.isfinite(quality[ru]).all(axis=1),
                both_within_parent_fractions_finite=np.isfinite(fraction[su])&np.isfinite(fraction[ru]))
            available=np.logical_and.reduce(list(flags.values()))
            for pos,(l,ok) in enumerate(zip(ids,available)):
                libraries.append(dict(record=item.record,context_id=ctx['context_id'],library_index=int(l),
                    measurement_and_receiver_gates_pass=bool(ok),sender_whole_programme_gate_required=False,
                    **{key:bool(value[pos]) for key,value in flags.items()},
                    failed_measurement_flags='|'.join(key for key,value in flags.items() if not value[pos])))
            base=[]
            for patient,g in members.iloc[np.flatnonzero(available)].groupby('resampling_block',sort=True):
                li=g.index.to_numpy(int);a,b=li*4+sv,li*4+rv
                assert g.reviewed_patient.nunique()==1
                if g.block_key.nunique()!=1:
                    mixed.append(dict(record=item.record,context_id=ctx['context_id'],resampling_block=str(patient),reason='MIXED_TISSUE_TREATMENT_CAPTURE_EXVIVO_BLOCKS'));continue
                row=dict(case_row=len(cases)+len(base),record=item.record,context_id=ctx['context_id'],resampling_block=str(patient),
                    reviewed_patient=bool(g.reviewed_patient.iloc[0]),block_key=g.block_key.iloc[0],paired_libraries=len(li),library_indices_json=json.dumps(li.tolist()),
                    ligand=float(z['mean_logCP10K'][a,gi].mean()),response=float(item.direction*response['response'][ri,li].mean()),receiver_own_ligand=float(z['mean_logCP10K'][b,gi].mean()),
                    fraction_within_parent_sender=float(fraction[a].mean()),fraction_within_parent_receiver=float(fraction[b].mean()),
                    fraction_whole_sender=float(np.mean(z['n_cells'][a]/denominator[li])),fraction_whole_receiver=float(np.mean(z['n_cells'][b]/denominator[li])),
                    sender_cells_sum=int(z['n_cells'][a].sum()),receiver_cells_sum=int(z['n_cells'][b].sum()),sender_ligand_positive=int(z['positive_counts'][a,gi].sum()),
                    **{f'Q_sender_{j}':float(quality[a,j].mean()) for j in range(5)},**{f'Q_receiver_{j}':float(quality[b,j].mean()) for j in range(5)})
                for j,col in enumerate(rcols):
                    row[f'receptor_{j}_measured']=bool(z['complex_measured'][b,col].all())
                    row[f'receptor_{j}_positive']=float(z['receptor_positive'][b,col].sum())
                    row[f'receptor_{j}_positive2']=float(z['receptor_positive2'][b,col].sum())
                base.append(row)
            cases.extend(base);frame=pd.DataFrame(base)
            for family in scope['families']:
                g=frame if family=='all_units' or frame.empty else frame[frame.reviewed_patient]
                if not g.empty:
                    counts=g.block_key.value_counts();g=g[g.block_key.map(counts).ge(2)]
                model=item.record+'|'+ctx['context_id']+'|'+family
                row=dict(model=model,record=item.record,context_id=ctx['context_id'],family=family,**{k:ctx[k] for k in ['arm','dataset','cancer_label','sample_class','material_class','phase_key','in_model_source']},
                    candidate_units=len(frame) if family=='all_units' else sum(x['reviewed_patient'] for x in base),active_units=len(g),active_reviewed_units=int(g.reviewed_patient.sum()) if len(g) else 0,
                    active_case_rows_json=json.dumps(g.case_row.tolist()) if len(g) else '[]',canonical_receptors_json=json.dumps(recs),RNA_qualified_receptors_json='[]',
                    pooled_sender_positive_fraction=np.nan,sender_positive_units=0,nuisance_rank=0,residual_df=0,test_index=-1,status='')
                reasons=[]
                if ctx['in_model_source']:reasons.append('IN_MODEL_SOURCE_NOT_EXTERNAL')
                if len(g)<scope['min_informative_units']:reasons.append('FEWER_THAN20_INFORMATIVE_UNITS')
                if len(g):
                    denominator_sender=int(g.sender_cells_sum.sum());fraction_sender=float(g.sender_ligand_positive.sum()/denominator_sender)
                    row.update(pooled_sender_positive_fraction=fraction_sender,sender_positive_units=int(g.sender_ligand_positive.gt(0).sum()))
                    if fraction_sender<.01 or row['sender_positive_units']<10:reasons.append('SENDER_POOLED1PCT_OR10_UNIT_RNA_GATE_FAILED')
                    qualified=[]
                    for j,rec in enumerate(recs):
                        if g[f'receptor_{j}_measured'].all() and g[f'receptor_{j}_positive'].sum()/g.receiver_cells_sum.sum()>=.01 and g[f'receptor_{j}_positive'].gt(0).sum()>=10:qualified.append(rec)
                    row['RNA_qualified_receptors_json']=json.dumps(qualified)
                    if not qualified:reasons.append('NO_COMPLETE_RECEPTOR_POOLED1PCT_AND10_UNIT_SUPPORT')
                    try:
                        p=kernel.prepare(g.ligand.to_numpy(float),g.response.to_numpy(float),g[QCOLS].to_numpy(float),g.block_key.to_numpy(str),min_n=20,min_df=8)
                        design_gate(p);assert p['keep'].all()
                        row.update(nuisance_rank=p['design_rank'],residual_df=p['residual_df'])
                    except ValueError as e:reasons.append('DESIGN:'+str(e))
                row['status']=';'.join(reasons) if reasons else 'ELIGIBLE_FIXED_EXTERNAL_PRIMARY'
                if not reasons:row['test_index']=counter;counter+=1
                records.append(row)
    eligibility=pd.DataFrame(records);case_frame=pd.DataFrame(cases) if cases else pd.DataFrame(columns=['case_row','record','context_id'])
    assert len(eligibility)==scope['planned_positions'] and eligibility.model.is_unique
    assert case_frame.case_row.tolist()==list(range(len(case_frame)))
    artifacts={'private_all_planned_model_eligibility.tsv.gz':eligibility,'private_all_case_measurement_units.tsv.gz':case_frame,
        'private_all_library_measurement_eligibility.tsv.gz':pd.DataFrame(libraries),
        'private_excluded_mixed_case_blocks.tsv.gz':pd.DataFrame(mixed,columns=['record','context_id','resampling_block','reason'])}
    for name,table in artifacts.items():table.to_csv(out/name,sep='\t',index=False,compression='gzip')
    sizes={f:int(((eligibility.family==f)&eligibility.test_index.ge(0)).sum()) for f in scope['families']}
    (out/'private_frozen_families_before_effects.json').write_text(json.dumps(dict(scope_sha256=SCOPE_SHA,all_planned_positions=len(eligibility),eligible_family_sizes=sizes,no_effects_or_P=True),indent=2)+'\n')
    files=[p for p in out.iterdir() if p.is_file()];assert sum(p.stat().st_size for p in files)<32*1024**2
    terminal=dict(status='COMPLETE_ALL_EXTERNAL_MODELS_AND_NEW_PAIR_UNITS_FROZEN_BEFORE_EFFECTS_AUDIT_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),scope_sha256=SCOPE_SHA,
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},new_effects_P_or_clinical_outcomes=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],completed_UTC=terminal['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
