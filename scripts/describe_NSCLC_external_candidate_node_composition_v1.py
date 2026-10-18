"""Post-result same-case full frozen-node composition diagnostic for ALL models.

No new P/CI, source selection, node removal or primary result change.
"""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_external_lineage_inputs_v2 import load as metadata_inputs
from NSCLC_native_MFA_interpretation_inputs_v1 import sha,VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_NSCLC_candidate_external_associations_v1 import inputs
from prepare_NSCLC_candidate_external_associations_v1 import QCOLS
import NSCLC_candidate_external_association_math_v1 as math


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    src=root/'NSCLC_candidate_external_association_results_20261018_v1'
    proof=src/'private_independent_external_candidate_result_audit.json'
    assert sha(proof)=='dc9e287dea77c2044f03f4b7087c3cc1462436f2d913088315e33736dd0a818b'
    terminal=json.loads((src/'private_terminal.json').read_text());name='private_all_estimable_external_effects_and_sensitivities.tsv.gz'
    assert sha(src/name)==terminal['files'][name]['sha256']
    result=pd.read_csv(src/name,sep='\t');_,_,contract,_,cases,_=inputs(root)
    candidates=pd.DataFrame(contract['candidates']).set_index('record')
    metadata=metadata_inputs();lib=metadata['lib'];counts=metadata['counts'];nu=len(lib)*4
    expected_root=root/'NSCLC_candidate_external_exact_RNA_20261018_v1'
    with np.load(expected_root/'private_all_frozen_pool_exact_panel_measurements.npz',allow_pickle=False) as z:expected=z['n_cells']
    raw_terminal=json.loads((expected_root/'private_terminal.json').read_text())
    assert metadata['hashes']==raw_terminal['frozen_inputs']['metadata']
    out=root/'NSCLC_external_candidate_node_composition_diagnostic_20261018_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        accepted_result_audit_sha256=sha(proof),metadata_sha256=metadata['hashes'],
        scope='Every estimable original model in both families, exact original case/library membership and full primary controls. Post-result explanatory diagnostic, not prospective confirmation.',
        fractions='All original frozen annotation nodes in both selected pools, including constant/zero columns. Exact library pool cell denominator, then equal same-case library average.',
        full_span='Retain all columns plus original twelve continuous covariates and full blocks. Same rank/machine/stability/df gates; no variable removal when unidentifiable.',
        interpretation='Nodes are frozen subtype/marker partitions, not formally assigned cell states. Attenuation can reflect real ecology and is not automatically artifact. No new P/CI or primary modification.')
    (out/'plan_before_all_model_composition_diagnostic.json').write_text(json.dumps(plan,indent=2)+'\n')
    lookup={(r.arm,r.dataset,r.library_id):i for i,r in enumerate(lib.itertuples())}
    matrices={};reconstructed=np.zeros(nu,np.int64)
    for arm,cat in metadata['catalogs'].items():
        cat=cat.reset_index(drop=True);node_to_col={x:i for i,x in enumerate(cat.node_id)}
        values=np.zeros((1690,len(cat)),np.int64)
        for row in counts[counts.arm.eq(arm)].itertuples():
            values[lookup[(arm,row.dataset,row.library_id)],node_to_col[row.node_id]]+=int(row.n_cells)
        for vi in range(4):
            take=cat.frozen_role.eq(vi).to_numpy();ids=lib.index[lib.arm.eq(arm)].to_numpy(int)
            reconstructed[ids*4+vi]=values[ids][:,take].sum(axis=1)
        matrices[arm]=(cat,values)
    assert np.array_equal(reconstructed,expected)
    rows=[]
    for item in result.itertuples():
        candidate=candidates.loc[item.record];sv,rv=VIEWS.index(candidate.sender_view),VIEWS.index(candidate.receiver_view)
        g=cases.iloc[json.loads(item.active_case_rows_json)];cat,values=matrices[item.arm]
        columns=[cat.frozen_role.eq(v).to_numpy() for v in [sv,rv]]
        extra=[]
        for case in g.itertuples():
            li=np.array(json.loads(case.library_indices_json),int)
            assert lib.loc[li,'arm'].eq(item.arm).all() and lib.loc[li,'dataset'].eq(item.dataset).all()
            parts=[]
            for vi,take in zip([sv,rv],columns):
                denominator=expected[li*4+vi]
                assert (denominator>=20).all() and np.array_equal(values[li][:,take].sum(axis=1),denominator)
                fractions=values[li][:,take]/denominator[:,None]
                assert np.allclose(fractions.sum(axis=1),1,atol=1e-12)
                parts.append(fractions.mean(axis=0))
            extra.append(np.concatenate(parts))
        extra=np.stack(extra);q=g[QCOLS].to_numpy(float);x,y,blocks=g.ligand.to_numpy(float),g.response.to_numpy(float),g.block_key.to_numpy(str)
        base=math.prepare(x,y,q,blocks);base_effect=math.base.effect(base)
        assert abs(base_effect-item.partial_rank_rho)<1e-10
        adjusted_q=np.column_stack([q,extra]);effect=np.nan;rank=df=0
        try:
            p=math.prepare(x,y,adjusted_q,blocks);effect=math.base.effect(p);rank=p['design_rank'];df=p['residual_df'];status='ESTIMABLE'
        except ValueError as e:status=str(e)
        rows.append(dict(model=item.model,record=item.record,family=item.family,arm=item.arm,dataset=item.dataset,
            cancer_label=item.cancer_label,material_class=item.material_class,phase_key=item.phase_key,
            exact_same_case_units=len(g),extra_frozen_node_columns=extra.shape[1],
            constant_node_columns=int(np.ptp(extra,axis=0).eq(0).sum()) if isinstance(extra,pd.DataFrame) else int((np.ptp(extra,axis=0)==0).sum()),
            original_partial_rho=base_effect,full_node_composition_adjusted_rho=effect,status=status,nuisance_rank=rank,residual_df=df,
            signed_change=effect-base_effect,absolute_magnitude_change=abs(effect)-abs(base_effect),
            same_sign_as_original=bool(np.isfinite(effect) and np.sign(effect)==np.sign(base_effect))))
    table=pd.DataFrame(rows);assert table.model.tolist()==result.model.tolist()
    table.to_csv(out/'private_all_original_models_same_case_node_composition.tsv.gz',sep='\t',index=False,compression='gzip')
    files=[p for p in out.iterdir() if p.is_file()]
    terminal=dict(status='COMPLETE_ALL_ORIGINAL_MODELS_POST_RESULT_SAME_CASE_FULL_NODE_COMPOSITION_NO_NEW_P',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),result_audit_sha256=sha(proof),
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},
        independently_raw_replayed=False,new_P_CI_or_primary_changes=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],completed_UTC=terminal['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
