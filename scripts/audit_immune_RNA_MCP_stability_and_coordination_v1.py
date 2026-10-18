"""Independent all-slot matching, source designs and complete random-stream replay."""
import itertools
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def equal(a, b, tolerance=1e-9):
    assert (pd.isna(a) and pd.isna(b)) or (np.isfinite(a) and np.isfinite(b) and abs(a-b) <= tolerance), 'NUMERIC_MISMATCH'


def rho(x, y):
    if len(x) < 3 or min(np.ptp(x), np.ptp(y)) == 0:
        return np.nan
    a = pd.Series(x).rank().to_numpy(); b = pd.Series(y).rank().to_numpy()
    a -= a.mean(); b -= b.mean()
    return float(a @ b / np.sqrt((a @ a)*(b @ b)))


def independent_fit(x, y, C):
    n = len(x)
    if n < 20: return None, 'FEWER_THAN20_COMPLETE_UNITS'
    if not np.isfinite(np.column_stack([x, y, C])).all(): return None, 'NONFINITE_COMPLETE_DESIGN'
    if min(np.ptp(x), np.ptp(y)) == 0: return None, 'CONSTANT_SCORE'
    rank_c = pd.DataFrame(C).rank().to_numpy()
    sd = rank_c.std(axis=0, ddof=1)
    if np.any(sd == 0): return None, 'CONSTANT_REQUIRED_COVARIATE_NO_AUTOMATIC_DELETION'
    D = np.column_stack([np.ones(n), (rank_c-rank_c.mean(axis=0))/sd])
    singular = np.linalg.svd(D, compute_uv=False)
    rank = np.count_nonzero(singular > singular[0]*max(D.shape)*np.finfo(float).eps)
    if rank != D.shape[1]: return None, 'COLLINEAR_REQUIRED_DESIGN_NO_AUTOMATIC_DELETION'
    condition = singular[0]/singular[-1]
    if not condition < 1e8: return None, 'ILL_CONDITIONED_REQUIRED_DESIGN'
    df = n-D.shape[1]-1
    if df < 10: return None, 'FEWER_THAN10_FULL_MODEL_RESIDUAL_DF'
    rx0 = pd.Series(x).rank().to_numpy(); ry0 = pd.Series(y).rank().to_numpy()
    fit_y = D @ np.linalg.lstsq(D, ry0, rcond=None)[0]
    rx = rx0-D @ np.linalg.lstsq(D, rx0, rcond=None)[0]; ry = ry0-fit_y
    norm = np.sqrt((rx @ rx)*(ry @ ry))
    if norm < 1e-12: return None, 'SCORE_EXPLAINED_BY_REQUIRED_COVARIATES'
    return dict(effect=float(rx @ ry/norm), D=D, rx=rx, ry=ry, fit_y=fit_y,
        residual_df=df, condition=float(condition)), 'ESTIMABLE'


def family_q(values, by=False):
    p = np.asarray(values, float); m = len(p)
    factor = sum(1/i for i in range(1, m+1)) if by else 1
    return np.asarray([min(1., min(m*factor*b/np.count_nonzero(p<=b) for b in p[p>=a])) for a in p])


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    cp=root/'immune_RNA_MCP_postfit_contract_20261018_v1/plan_before_learned_outputs_or_projection.json'
    assert sha(cp)=='1172d887cffb7c3b0c28a21583154fbe3718a1f13d79edb85500343ead5e2b09'
    contract=json.loads(cp.read_text())
    execution=root/'immune_RNA_MCP_accepted_postfit_execution_20261018_v1/private_terminal.json'
    execution_proof=json.loads(execution.read_text())
    assert execution_proof['status']=='FIXED_STABILITY_AND_SOURCE_INFERENCE_DONE_REAL_NUMERIC_ACCEPTANCE_PENDING'
    for dirname,expected in [
            ('immune_RNA_MCP_discovery_scope_20261018_v1','0f4833a07bd5cdf5e8f6806d06805777dddf78adad60a6151eca991dfdda6196'),
            ('immune_RNA_MCP_native_expression_inputs_20261018_v1','2474355a003d96232b752ad2c77ba92a506f4c68829c7e0adf434045cf03d656')]:
        parent=root/dirname
        assert sha(parent/'private_terminal.json')==expected
        for name,item in json.loads((parent/'private_terminal.json').read_text())['files'].items():
            assert sha(parent/name)==item['sha256']
    stability=root/'immune_RNA_MCP_seed_stability_20261018_v1'
    source=root/'immune_RNA_MCP_source_coordination_20261018_v1'
    for parent in [stability,source]:
        terminal=json.loads((parent/'private_terminal.json').read_text())
        for name,digest in terminal['files'].items(): assert sha(parent/name)==digest
    assert [s['terminal_sha256'] for s in execution_proof['stages']]==[
        sha(stability/'private_terminal.json'),sha(source/'private_terminal.json')]
    projection=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection'
    acceptance=root/'immune_RNA_MCP_frozen_projection_audit_20261018_v1/private_terminal.json'
    proof=json.loads(acceptance.read_text())
    assert proof['status']=='ALL_PLANNED_SLOTS_AND_ALL_SUPPORTED_PROJECTION_VALUES_ACCEPTED'
    assert sha(projection.parent/'private_terminal.json')==proof['producer_terminal_sha256']
    for name,digest in json.loads((projection.parent/'private_terminal.json').read_text())['files'].items():
        assert sha(projection.parent/name)==digest
    output=root/'immune_RNA_MCP_stability_coordination_audit_20261018_v1'
    assert not output.exists();output.mkdir(mode=0o700)
    (output/'plan_before_independent_all_slot_replay.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        stability_terminal_sha256=sha(stability/'private_terminal.json'),source_terminal_sha256=sha(source/'private_terminal.json'),
        scope='All54 comparisons/24 matches/2 engineering labels by source-weighted ranks, sparse overlap cosine and exhaustive assignment. All6 heldout tests,6 reviewed sensitivities,12 discovery descriptions. Original input metadata, all19999 residual permutations and all2000 bootstrap values per available slot; independent OLS, pandas ranks/SVD and direct family-q formula.'),indent=2)+'\n')
    score_path=projection/'private_all_case_mean_scores.tsv'
    scores=pd.read_csv(score_path,sep='\t',dtype={'sample':str}) if score_path.exists() else pd.DataFrame(
        columns=['fit_id','view','program','partition','sample','dataset','score'])
    scores['score']=scores.score.astype(float)
    wp=projection/'private_frozen_gene_weights.tsv'
    weights=pd.read_csv(wp,sep='\t') if wp.exists() else pd.DataFrame(columns=['fit_id','view','program','gene','weight'])
    full=pd.read_csv(stability/'private_all54_seed_rank_source_pair_comparisons.tsv',sep='\t')
    matches=pd.read_csv(stability/'private_all24_joint_matches_by_view.tsv',sep='\t')
    classes=pd.read_csv(stability/'private_all2_anchor_engineering_classification.tsv',sep='\t')
    assert len(full)==54 and len(matches)==24 and len(classes)==2
    views=list(contract['view_to_official_name']);anchor=contract['anchor'];pair_keys=[]
    for item in contract['planned_fits']:
        fit_id='K%d_seed%d_%s'%(item['k'],item['seed'],item['source_mode'])
        if fit_id==anchor: continue
        pair_keys.extend((fit_id,a,b,v) for a,b,v in itertools.product(['MCP1','MCP2'],
            ['MCP%d'%i for i in range(1,item['k']+1)],views))
    assert set(map(tuple,full[['comparator','anchor_program','comparator_program','view']].to_numpy()))==set(pair_keys)
    for r in full.itertuples(index=False):
        a=scores[scores.fit_id.eq(anchor)&scores.view.eq(r.view)&scores.program.eq(r.anchor_program)&scores.partition.eq('discovery')]
        b=scores[scores.fit_id.eq(r.comparator)&scores.view.eq(r.view)&scores.program.eq(r.comparator_program)&scores.partition.eq('discovery')]
        paired=a.merge(b,on=['sample','dataset'],suffixes=('_A','_B'),validate='one_to_one')
        expected=np.nan;by_source={}
        if len(paired)==38:
            for s,g in paired.groupby('dataset'):
                by_source[s]=(len(g)-1,rho(g.score_A.to_numpy(),g.score_B.to_numpy()))
            if len(by_source)==2 and all(np.isfinite(v[1]) for v in by_source.values()):
                expected=sum(n*v for n,v in by_source.values())/sum(n for n,v in by_source.values())
            else:
                by_source={}
        equal(expected,r.within_source_centred_rank_correlation)
        if by_source:
            equal(by_source['GSE286827'][1],r.HNSCC_rank_correlation)
            equal(by_source['GSE294109'][1],r.ccRCC_rank_correlation)
        wa=weights[weights.fit_id.eq(anchor)&weights.view.eq(r.view)&weights.program.eq(r.anchor_program)].set_index('gene').weight
        wb=weights[weights.fit_id.eq(r.comparator)&weights.view.eq(r.view)&weights.program.eq(r.comparator_program)].set_index('gene').weight
        common=wa.index.intersection(wb.index)
        cosine=float(wa.loc[common]@wb.loc[common]/np.sqrt((wa@wa)*(wb@wb))) if len(wa) and len(wb) else np.nan
        equal(cosine,r.signed_weight_cosine)
    for fit_id,g in full.groupby('comparator'):
        candidates=sorted(g.comparator_program.unique());cost={}
        for a,b in itertools.product(['MCP1','MCP2'],candidates):
            z=g[g.anchor_program.eq(a)&g.comparator_program.eq(b)]
            valid=z.status.eq('AVAILABLE').all()
            mean=z.within_source_centred_rank_correlation.mean() if valid else np.nan
            cost[a,b]=-abs(mean) if np.isfinite(mean) else 1e6
        selected=matches[matches.comparator.eq(fit_id)]
        pairs=selected[['anchor_program','comparator_program']].drop_duplicates()
        assert len(pairs)==2 and pairs.anchor_program.nunique()==2 and pairs.comparator_program.nunique()==2
        actual=sum(cost[a,b] for a,b in pairs.to_numpy())
        minimum=min(cost['MCP1',a]+cost['MCP2',b] for a,b in itertools.permutations(candidates,2))
        equal(actual,minimum)
        for r in selected.itertuples(index=False):
            z=g[g.anchor_program.eq(r.anchor_program)&g.comparator_program.eq(r.comparator_program)]
            valid=z.status.eq('AVAILABLE').all();mean=z.within_source_centred_rank_correlation.mean() if valid else np.nan
            sign=1 if not np.isfinite(mean) or mean>=0 else -1
            assert r.global_orientation==sign and r.matching_available==valid
            original=z[z.view.eq(r.view)].iloc[0]
            equal(r.aligned_correlation,sign*original.within_source_centred_rank_correlation)
            equal(r.aligned_weight_cosine,sign*original.signed_weight_cosine)
    for r in classes.itertuples(index=False):
        g=matches[matches.anchor_program.eq(r.program)&matches.comparator_role.eq('seed_stability')]
        assert len(g)==6
        available=g.matching_available.all() and np.isfinite(g.aligned_correlation).all() and np.isfinite(g.aligned_weight_cosine).all()
        stable=available and g.aligned_correlation.ge(.70).all() and g.aligned_weight_cosine.ge(.50).all()
        assert r.all_seed_views_available==available and r.prespecified_engineering_stability_pass==stable
    mdrows=[]
    case_scope=pd.read_csv(root/'immune_RNA_MCP_discovery_scope_20261018_v1/private_metadata_only_case_source_partition.tsv',sep='\t')
    reviewed=case_scope.set_index('fit_sample').reviewed_patient
    for view,part in itertools.product(views,['discovery','TNBC_holdout']):
        md=pd.read_csv(root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_inputs'/(view+'_'+part+'_metadata.tsv'),sep='\t',dtype={'sample':str})
        for sample,g in md.groupby('sample',sort=True):
            mdrows.append(dict(sample=sample,dataset=g.dataset.iloc[0],partition=part,view=view,
                reviewed_patient=bool(reviewed.loc[sample]),fraction=g.fraction.iloc[0],logN=g.logN.iloc[0],
                cellQ=float(np.average(g.cellQ)),geneQ=float(np.average(g.geneQ)),selected_cells=len(g)))
    metadata=pd.DataFrame(mdrows);assert len(metadata)==198
    saved=pd.read_csv(source/'private_all_case_view_controls.tsv',sep='\t',dtype={'sample':str}).set_index(['sample','view'])
    for r in metadata.itertuples(index=False):
        s=saved.loc[(r.sample,r.view)]
        assert s.dataset==r.dataset and s.partition==r.partition and s.reviewed_patient==r.reviewed_patient
        for c in ['fraction','logN','cellQ','geneQ','selected_cells']:equal(getattr(r,c),s[c])
    main=pd.read_csv(source/'private_all6_primary_source_tests.tsv',sep='\t')
    sensitivity=pd.read_csv(source/'private_all6_reviewed_patient_sensitivity.tsv',sep='\t')
    discovery=pd.read_csv(source/'private_all12_fitted_source_descriptions.tsv',sep='\t')
    assert len(main)==len(sensitivity)==6 and len(discovery)==12
    actual_perms=0;actual_boot=0;maximum_error=0.
    for slot,item in enumerate(contract['primary_source_association']['slots']):
        sides=[]
        for view in [item['view_A'],item['view_B']]:
            t=metadata[metadata.view.eq(view)].copy()
            s=scores[scores.fit_id.eq(anchor)&scores.view.eq(view)&scores.program.eq(item['program'])]
            sides.append(t.merge(s[['sample','score']],on='sample',how='left',validate='one_to_one'))
        pair=sides[0].merge(sides[1],on=['sample','dataset','partition','reviewed_patient'],suffixes=('_A','_B'),validate='one_to_one')
        cols=[c+'_'+side for side in ['A','B'] for c in ['fraction','logN','cellQ','geneQ']]
        held=pair[pair.partition.eq('TNBC_holdout')]
        for label,table in [('main_sample_cases',held),('reviewed_patient_only',held[held.reviewed_patient])]:
            rows=main if label=='main_sample_cases' else sensitivity
            r=rows[rows.slot.eq(slot)].iloc[0]
            valid=np.isfinite(table[['score_A','score_B']+cols].to_numpy(float)).all(axis=1);d=table[valid]
            x=d.score_A.to_numpy();y=d.score_B.to_numpy();C=d[cols].to_numpy(float)
            f,reason=independent_fit(x,y,C)
            assert r.status==reason and r.complete_units==len(d) and r.eligible_units==len(table)
            assert r.reviewed_patient_units==d.reviewed_patient.sum()
            equal(r.unadjusted_rho,rho(x,y))
            equal(r.adjusted_rho,f['effect'] if f is not None else np.nan)
            if f is None:
                assert pd.isna(r.permutation_p) and pd.isna(r.lower) and pd.isna(r.upper)
                if label=='main_sample_cases':
                    assert not (source/('private_slot_%02d_resampling.npz'%slot)).exists()
            if label!='main_sample_cases' or f is None:continue
            equal(r.design_condition,f['condition'],1e-6);assert r.residual_df==f['residual_df']
            z=np.load(source/('private_slot_%02d_resampling.npz'%slot),allow_pickle=False)
            assert z['samples'].tolist()==d['sample'].tolist()
            assert np.allclose(z['x'],x) and np.allclose(z['y'],y) and np.allclose(z['controls'],C)
            n=len(d);seed=contract['primary_source_association']['seed']
            rng=np.random.default_rng(np.random.SeedSequence([seed,slot,0]))
            indices=np.stack([rng.permutation(n) for _ in range(19999)])
            assert np.array_equal(z['permutation_indices'],indices)
            outcomes=f['fit_y'][:,None]+f['ry'][indices].T
            residual=outcomes-f['D']@np.linalg.lstsq(f['D'],outcomes,rcond=None)[0]
            perm=f['rx']@residual/np.sqrt((f['rx']@f['rx'])*np.sum(residual**2,axis=0))
            error=float(np.max(np.abs(perm-z['permutation_effects'])));assert error<1e-9
            maximum_error=max(maximum_error,error);actual_perms+=len(perm)
            p=(1+np.count_nonzero(np.abs(perm)>=abs(f['effect'])-1e-12))/20000
            equal(p,r.permutation_p)
            rng=np.random.default_rng(np.random.SeedSequence([seed,slot,1]))
            bootstrap=rng.integers(0,n,size=(2000,n),dtype=np.int16)
            assert np.array_equal(z['bootstrap_indices'],bootstrap)
            values=[];reasons=[]
            for take in bootstrap:
                b,reason=independent_fit(x[take],y[take],C[take]);reasons.append(reason)
                values.append(b['effect'] if b is not None else np.nan)
            values=np.asarray(values);assert np.array_equal(z['bootstrap_reasons'],reasons)
            assert np.allclose(z['bootstrap_effects'],values,atol=1e-9,rtol=1e-9,equal_nan=True)
            keep=np.isfinite(values);assert r.valid_bootstrap==keep.sum()
            ci=np.quantile(values[keep],[.025,.975]) if keep.sum()>=1000 else [np.nan,np.nan]
            equal(ci[0],r.lower);equal(ci[1],r.upper);actual_boot+=2000
        for s,g in pair[pair.partition.eq('discovery')].groupby('dataset',sort=True):
            r=discovery[discovery.slot.eq(slot)&discovery.dataset.eq(s)].iloc[0]
            valid=np.isfinite(g[['score_A','score_B']+cols].to_numpy(float)).all(axis=1);d=g[valid]
            x=d.score_A.to_numpy();y=d.score_B.to_numpy();C=d[cols].to_numpy(float)
            f,reason=independent_fit(x,y,C)
            assert r.eligible_units==len(g) and r.complete_units==len(d) and r.adjusted_status==reason
            equal(r.unadjusted_rho,rho(x,y));equal(r.adjusted_rho,f['effect'] if f is not None else np.nan)
            leave=[rho(np.delete(x,i),np.delete(y,i)) for i in range(len(d))]
            leave=[v for v in leave if np.isfinite(v)]
            equal(r.leave_one_unadjusted_min,min(leave) if leave else np.nan)
            equal(r.leave_one_unadjusted_max,max(leave) if leave else np.nan)
    p=main.permutation_p.fillna(1).to_numpy()
    assert np.allclose(p,main.family_adjustment_p)
    for name,by in [('BH6',False),('BY6',True)]:
        q=family_q(p,by);q[main.permutation_p.isna()]=np.nan
        assert np.allclose(q,main[name],equal_nan=True)
    record=dict(status='ACCEPTED_ALL_STABILITY_SOURCE_DESIGNS_AND_COMPLETE_RANDOM_STREAMS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        stability_terminal_sha256=sha(stability/'private_terminal.json'),source_terminal_sha256=sha(source/'private_terminal.json'),
        pair_view_comparisons=54,matched_view_rows=24,primary_slots=6,reviewed_sensitivities=6,
        discovery_descriptions=12,all_replayed_permutation_values=actual_perms,all_replayed_bootstrap_values=actual_boot,
        maximum_permutation_difference=maximum_error,not_independent_biological_or_clinical_validation=True)
    (output/'private_terminal.json').write_text(json.dumps(record,indent=2)+'\n')
    print('ALL_STABILITY_SOURCE_AND_RANDOM_STREAM_VALUES_ACCEPTED')


if __name__=='__main__':main()
