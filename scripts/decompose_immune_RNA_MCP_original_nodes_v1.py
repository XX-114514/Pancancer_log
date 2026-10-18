"""Frozen discovery-node composition baseline; complete six-slot follow-up."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_immune_RNA_MCP_source_coordination_v1 import fit_rank_partial, resample, adjust_family


def node_baseline(frame):
    """Equal-case WLS within transformation; independently solve full FE design."""
    train = frame[~frame.dataset.eq('GSE305078')].copy()
    assert train['sample'].nunique() == 38
    nodes = sorted(train.node_id.unique())
    reference = train.groupby('node_id').selected_cells.sum().sort_values(
        ascending=False, kind='stable').index[0]
    other = [n for n in nodes if n != reference]
    cases = sorted(train['sample'].unique())
    case_index = pd.Index(cases).get_indexer(train['sample'])
    node_index = pd.Index(nodes).get_indexer(train.node_id)
    one_node = np.eye(len(nodes))[node_index]
    one_case = np.eye(len(cases))[case_index]
    size = train.groupby('sample').selected_cells.transform('sum').to_numpy(float)
    w = train.selected_cells.to_numpy(float) / size
    fractions = one_case.T @ (w[:, None] * one_node)
    assert np.allclose(fractions.sum(axis=1), 1.)
    y = train.mean_score.to_numpy(float)
    ymean = one_case.T @ (w * y)
    cols = pd.Index(nodes).get_indexer(other)
    X = (one_node - fractions[case_index])[:, cols]
    A = X * np.sqrt(w[:, None]); b = (y-ymean[case_index])*np.sqrt(w)
    rank = np.linalg.matrix_rank(A); condition = float(np.linalg.cond(A))
    if rank != len(other) or condition >= 1e8:
        return None, dict(status='UNIDENTIFIABLE_DISCOVERY_NODE_BASELINE', rank=int(rank),
            columns=len(other), condition=condition)
    beta = np.zeros(len(nodes)); beta[cols] = np.linalg.lstsq(A, b, rcond=None)[0]
    full = np.column_stack([one_case, one_node[:, cols]]) * np.sqrt(w[:, None])
    assert np.linalg.matrix_rank(full) == full.shape[1]
    independent = np.linalg.lstsq(full, y*np.sqrt(w), rcond=None)[0][len(cases):]
    error = float(np.max(np.abs(beta[cols]-independent)))
    assert error < 1e-8
    # Fix only the arbitrary additive origin; all source ranks are invariant to it.
    beta -= np.mean(fractions @ beta)
    result = pd.Series(beta, index=nodes)
    counts = train.groupby('node_id').agg(training_cells=('selected_cells','sum'),
        training_cases=('sample','nunique')).reindex(nodes)
    counts['coefficient'] = result
    counts.index.name = 'node_id'
    return counts.reset_index(), dict(status='ESTIMABLE', training_cases=38,
        training_nodes=len(nodes), columns=len(other), rank=int(rank), condition=condition,
        independent_full_design_max_error=error,
        rare_nodes_below20_cells_or3_cases=int(((counts.training_cells<20)|(counts.training_cases<3)).sum()))


def point(x, y, C):
    fitted, reason = fit_rank_partial(x, y, C)
    return (fitted['effect'] if fitted is not None else np.nan), reason


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 64*1024**2)
    output = root/'immune_RNA_MCP_original_node_decomposition_20261018_v1'
    assert not output.exists(); output.mkdir(mode=0o700)
    loc = root/'immune_RNA_MCP_original_node_localization_20261018_v1'
    assert sha(loc/'private_terminal.json') == 'fc0675de85f514ca94b73d063e0e6d1ebbb9f27269061edf9a3b7b439e3e77a3'
    for name,digest in json.loads((loc/'private_terminal.json').read_text())['files'].items():
        assert sha(loc/name) == digest
    source = root/'immune_RNA_MCP_source_coordination_20261018_v1'
    assert sha(source/'private_terminal.json') == '72fdf3928ac53318e2c5b4db162d4e9c3c432b76f5105191e463a02ee73d7b23'
    for name,digest in json.loads((source/'private_terminal.json').read_text())['files'].items():
        assert sha(source/name) == digest
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        budget=budget, source_acceptance_sha256='27efdf31dfee1a3937de6790a03947ff13ab7af594f8a74ddd3e5ffeeaf11a6a',
        scope='All two fixed anchor programmes, three capture views and six endpoint pairs; no model refit or identity promotion.',
        baseline='Discovery38 only: equal-case weighted case fixed effects plus all observed original node effects. Independent full FE solution. Unseen heldout node contribution is missing, never zero.',
        primary='TNBC source follow-up: original score rank association adjusting original eight controls plus both frozen discovery-node composition predictions. >=20 complete units and >=10 full residual df, no covariate deletion.',
        inference=dict(seed=2026101809, residual_permutations=19999, block_bootstrap=2000,
            family='six additional follow-up tests, separate BH6/BY6; same heldout source, not independent validation'),
        sensitivities='Same-case original eight-control effect; residual-expression effect; complete full-native node-mixture effect. Point-only, no new P/CI/BH.',
        limits='Discovery-node coefficient uncertainty not included in conditional intervals. Sparse training-node coefficients and unknown heldout nodes explicit. Old labels may contain unresolved states or finer identities; residual is not proof of a formal state. TNBC units are samples pending reviewed patient identity.')
    (output/'plan_before_any_node_baseline_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    m = pd.read_csv(loc/'private_selected_case_node_score_means.tsv',sep='\t',dtype={'sample':str,'node_id':str})
    full = pd.read_csv(loc/'private_all2970_full_view_node_compositions.tsv',sep='\t',dtype={'sample':str,'node_id':str})
    metadata = pd.read_csv(source/'private_all_case_view_controls.tsv',sep='\t',dtype={'sample':str})
    assert not metadata.duplicated(['sample','view']).any()
    contract = json.loads((root/'immune_RNA_MCP_postfit_contract_20261018_v1/plan_before_learned_outputs_or_projection.json').read_text())
    scores=[]; parameters=[]; designs=[]
    for (view,program),g in m.groupby(['view','program'],sort=True):
        fitted, proof = node_baseline(g)
        designs.append(dict(view=view,program=program,**proof))
        if fitted is not None:
            parameters.append(fitted.assign(view=view,program=program))
            effects = fitted.set_index('node_id').coefficient
        for (sample,dataset),one in g.groupby(['sample','dataset'],sort=True):
            n = int(one.selected_cells.sum())
            total = float(np.average(one.mean_score,weights=one.selected_cells))
            part = one.node_id.map(effects) if fitted is not None else pd.Series(np.nan,index=one.index)
            unseen = int(one.loc[part.isna(),'selected_cells'].sum())
            comp = float(np.average(part,weights=one.selected_cells)) if part.notna().all() else np.nan
            f = full[full['sample'].eq(sample)&full.view.eq(view)&full.full_node_cells.gt(0)]
            assert len(f)>0 and np.isclose(f.fraction.sum(),1.)
            fp = f.node_id.map(effects) if fitted is not None else pd.Series(np.nan,index=f.index)
            full_unseen = int(f.loc[fp.isna(),'full_node_cells'].sum())
            full_comp = float(np.sum(fp*f.fraction)) if fp.notna().all() else np.nan
            scores.append(dict(sample=sample,dataset=dataset,view=view,program=program,
                selected_cells=n,score=total,composition_prediction=comp,
                residual_expression=total-comp,full_composition_prediction=full_comp,
                unseen_selected_node_cells=unseen,unseen_full_node_cells=full_unseen))
    score = pd.DataFrame(scores);assert len(score)==396
    pd.concat(parameters).to_csv(output/'private_frozen_discovery_node_coefficients.tsv',sep='\t',index=False)
    pd.DataFrame(designs).to_csv(output/'private_all6_discovery_node_designs.tsv',sep='\t',index=False)
    score.to_csv(output/'private_all396_case_view_decompositions.tsv',sep='\t',index=False)
    rows=[]; sensitivities=[]; descriptions=[]
    for slot,item in enumerate(contract['primary_source_association']['slots']):
        a,b,program=item['view_A'],item['view_B'],item['program']
        sides=[]
        for view in [a,b]:
            s=score[score.view.eq(view)&score.program.eq(program)]
            t=metadata[metadata.view.eq(view)].merge(s,on=['sample','dataset','view'],validate='one_to_one')
            sides.append(t)
        joined=sides[0].merge(sides[1],on=['sample','dataset','partition','reviewed_patient'],suffixes=('_A','_B'),validate='one_to_one')
        oldcontrols=[c+'_'+side for side in ['A','B'] for c in ['fraction','logN','cellQ','geneQ']]
        controls=oldcontrols+['composition_prediction_A','composition_prediction_B']
        held=joined[joined.dataset.eq('GSE305078')]
        complete=np.isfinite(held[['score_A','score_B']+controls].to_numpy(float)).all(axis=1)
        d=held[complete];x=d.score_A.to_numpy();y=d.score_B.to_numpy();C=d[controls].to_numpy(float)
        fit,reason=fit_rank_partial(x,y,C)
        old,oldreason=point(x,y,d[oldcontrols].to_numpy(float))
        row=dict(slot=slot,program=program,view_A=a,view_B=b,dataset='GSE305078',
            eligible_units=len(held),complete_units=len(d),unavailable_baseline_units=int((~complete).sum()),
            reviewed_patient_units=int(d.reviewed_patient.sum()),status=reason,
            unadjusted_rho=float(spearmanr(x,y).statistic) if len(d)>=3 else np.nan,
            original_adjusted_same_cases=old,original_adjusted_status=oldreason,
            adjusted_for_frozen_node_composition=np.nan,lower=np.nan,upper=np.nan,
            permutation_p=np.nan,valid_bootstrap=0,residual_df=np.nan)
        if fit is not None:
            row.update(adjusted_for_frozen_node_composition=fit['effect'],residual_df=fit['residual_df'])
            result,draws=resample(x,y,C,fit,slot,2026101809);row.update(result)
            np.savez_compressed(output/('private_slot_%02d_resampling.npz'%slot),
                samples=d['sample'].to_numpy(str),x=x,y=y,controls=C,control_names=np.asarray(controls),**draws)
        rows.append(row)
        for mode in ['subtracted_node_expression','full_native_node_composition']:
            if mode=='subtracted_node_expression':
                fields=['residual_expression_A','residual_expression_B']+oldcontrols
            else:
                fields=['score_A','score_B']+oldcontrols+['full_composition_prediction_A','full_composition_prediction_B']
            valid=np.isfinite(held[fields].to_numpy(float)).all(axis=1);z=held[valid]
            effect,status=point(z[fields[0]].to_numpy(),z[fields[1]].to_numpy(),z[fields[2:]].to_numpy(float))
            sensitivities.append(dict(slot=slot,program=program,view_A=a,view_B=b,dataset='GSE305078',mode=mode,
                complete_units=len(z),status=status,point_effect=effect,no_P_or_CI=True))
        for dataset,z in joined.groupby('dataset',sort=True):
            valid=np.isfinite(z[['score_A','score_B']+controls].to_numpy(float)).all(axis=1);z=z[valid]
            effect,status=point(z.score_A.to_numpy(),z.score_B.to_numpy(),z[controls].to_numpy(float))
            descriptions.append(dict(slot=slot,program=program,view_A=a,view_B=b,dataset=dataset,
                complete_units=len(z),status=status,point_effect=effect,fitted_source=dataset!='GSE305078',no_P=True))
    result=pd.DataFrame(rows);assert len(result)==6 and len(sensitivities)==12 and len(descriptions)==18
    p=result.permutation_p.fillna(1).to_numpy(float)
    result['BH6']=adjust_family(p);result['BY6']=adjust_family(p,by=True)
    result.loc[result.permutation_p.isna(),['BH6','BY6']]=np.nan
    result.to_csv(output/'private_all6_composition_followup_tests.tsv',sep='\t',index=False)
    pd.DataFrame(sensitivities).to_csv(output/'private_all12_point_sensitivities.tsv',sep='\t',index=False)
    pd.DataFrame(descriptions).to_csv(output/'private_all18_source_descriptions.tsv',sep='\t',index=False)
    files={p.name:sha(p) for p in sorted(output.iterdir()) if p.is_file()}
    (output/'private_terminal.json').write_text(json.dumps(dict(
        status='COMPLETE_ALL6_NODE_COMPOSITION_FOLLOWUP_TESTS_NUMERIC_REVIEW_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files,
        case_view_decompositions=396,primary_followup_slots=6,point_sensitivities=12,
        source_descriptions=18,estimable_slots=int(result.status.eq('ESTIMABLE').sum()),
        identity_promoted=False,directional_or_clinical_validation=False),indent=2)+'\n')
    print('ALL_NODE_COMPOSITION_FOLLOWUP_SLOTS_COMPLETE_REVIEW_PENDING')


if __name__=='__main__':main()
