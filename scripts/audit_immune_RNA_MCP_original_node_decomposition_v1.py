"""Replay discovery-node baselines and every follow-up random value independently."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from audit_immune_RNA_MCP_stability_and_coordination_v1 import independent_fit, family_q, rho, equal


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,8*1024**2)
    loc=root/'immune_RNA_MCP_original_node_localization_20261018_v1'
    output=root/'immune_RNA_MCP_original_node_decomposition_20261018_v1'
    terminal=json.loads((output/'private_terminal.json').read_text())
    assert terminal['source_sha256']=='aabbbda88d497cf52da175463b2770761f73c9cf733b3ad38828cfdcb343dd6d'
    for name,digest in terminal['files'].items():assert sha(output/name)==digest
    audit=root/'immune_RNA_MCP_original_node_decomposition_audit_20261018_v1'
    assert not audit.exists();audit.mkdir(mode=0o700)
    (audit/'plan_before_full_decomposition_numeric_replay.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,
        producer_terminal_sha256=sha(output/'private_terminal.json'),
        scope='All6 discovery-node models reconstructed from accepted anchor cell scores and original node mapping by full QR fixed effects. All396 decompositions,6 source tests,12 point sensitivities,18 descriptions,119994 permutations and12000 bootstraps; no scientific interpretation before acceptance.'),indent=2)+'\n')
    units=pd.read_csv(root/'existing_immune_fine_native_measurement_scope_20261016_v1/private_all_frozen_T_myeloid_node_units.tsv',sep='\t')
    mapping=np.load(root/'existing_immune_fine_native_measurement_scope_20261016_v1/private_original_cell_to_fine_node_unit.npy',mmap_mode='r')
    selection=pd.read_csv(root/'immune_RNA_MCP_native_expression_inputs_20261018_v1/private_frozen_original_cell_selection.tsv',sep='\t')
    lookup=units.set_index('unit_index').node_id
    selection['node_id']=lookup.loc[mapping[selection.global_row.to_numpy(int)]].to_numpy()
    cell=pd.read_csv(root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_all_cell_score_reconstruction.tsv.gz',sep='\t')
    cell=cell[cell.fit_id.eq('K2_seed20261018_source_adjusted')].merge(selection[['cell','view','node_id']],on=['cell','view'],validate='many_to_one')
    assert len(cell)==38036
    grouped=cell.groupby(['sample','dataset','view','program','node_id']).agg(
        selected_cells=('score','size'),mean_score=('score','mean')).reset_index()
    native=pd.read_csv(loc/'private_all2970_full_view_node_compositions.tsv',sep='\t')
    values=[];parameter_checks=0
    saved_parameters=pd.read_csv(output/'private_frozen_discovery_node_coefficients.tsv',sep='\t')
    for (view,program),g in grouped.groupby(['view','program'],sort=True):
        t=g[~g.dataset.eq('GSE305078')].copy();cases=sorted(t['sample'].unique());nodes=sorted(t.node_id.unique())
        reference=t.groupby('node_id').selected_cells.sum().sort_values(ascending=False,kind='stable').index[0]
        others=[n for n in nodes if n!=reference]
        one=pd.get_dummies(pd.Categorical(t['sample'],categories=cases),dtype=float).to_numpy()
        node=pd.get_dummies(pd.Categorical(t.node_id,categories=others),dtype=float).to_numpy()
        weight=t.selected_cells.to_numpy(float)/t.groupby('sample').selected_cells.transform('sum').to_numpy(float)
        design=np.column_stack([one,node])*np.sqrt(weight[:,None])
        assert np.linalg.matrix_rank(design)==design.shape[1]
        Q,R=np.linalg.qr(design,mode='reduced')
        coefficients=np.linalg.solve(R,Q.T@(t.mean_score.to_numpy(float)*np.sqrt(weight)))
        beta=pd.Series(0.,index=nodes);beta.loc[others]=coefficients[len(cases):]
        origin=np.sum(weight*t.node_id.map(beta))/len(cases);beta-=origin
        s=saved_parameters[saved_parameters.view.eq(view)&saved_parameters.program.eq(program)].set_index('node_id')
        assert set(s.index)==set(beta.index)
        assert np.allclose(beta.loc[s.index],s.coefficient,atol=1e-8,rtol=0)
        parameter_checks+=len(beta)
        for (sample,dataset),h in g.groupby(['sample','dataset'],sort=True):
            b=h.node_id.map(beta);n=h.selected_cells.to_numpy(float);missing=b.isna()
            total=float(h.mean_score.to_numpy()@n/n.sum())
            comp=float(b.to_numpy()@n/n.sum()) if not missing.any() else np.nan
            f=native[native['sample'].eq(sample)&native.view.eq(view)&native.full_node_cells.gt(0)]
            fb=f.node_id.map(beta)
            fc=float(np.sum(fb*f.fraction)) if fb.notna().all() else np.nan
            values.append(dict(sample=sample,dataset=dataset,view=view,program=program,score=total,
                composition_prediction=comp,residual_expression=total-comp,full_composition_prediction=fc,
                unseen_selected_node_cells=int(n[missing].sum()),unseen_full_node_cells=int(f.loc[fb.isna(),'full_node_cells'].sum())))
    recomputed=pd.DataFrame(values);assert len(recomputed)==396
    saved=pd.read_csv(output/'private_all396_case_view_decompositions.tsv',sep='\t').set_index(['sample','view','program'])
    for r in recomputed.itertuples(index=False):
        s=saved.loc[(r.sample,r.view,r.program)];assert s.dataset==r.dataset
        for c in ['score','composition_prediction','residual_expression','full_composition_prediction','unseen_selected_node_cells','unseen_full_node_cells']:
            equal(getattr(r,c),s[c],1e-8)
    metadata=pd.read_csv(root/'immune_RNA_MCP_source_coordination_20261018_v1/private_all_case_view_controls.tsv',sep='\t')
    main=pd.read_csv(output/'private_all6_composition_followup_tests.tsv',sep='\t')
    sensitivities=pd.read_csv(output/'private_all12_point_sensitivities.tsv',sep='\t')
    descriptions=pd.read_csv(output/'private_all18_source_descriptions.tsv',sep='\t')
    assert len(main)==6 and len(sensitivities)==12 and len(descriptions)==18
    permutation_values=0;bootstrap_values=0;maximum_error=0.
    for r in main.itertuples(index=False):
        sides=[]
        for view in [r.view_A,r.view_B]:
            s=recomputed[recomputed.view.eq(view)&recomputed.program.eq(r.program)]
            sides.append(metadata[metadata.view.eq(view)].merge(s,on=['sample','dataset','view'],validate='one_to_one'))
        both=sides[0].merge(sides[1],on=['sample','dataset','partition','reviewed_patient'],suffixes=('_A','_B'),validate='one_to_one')
        base=[c+'_'+side for side in ['A','B'] for c in ['fraction','logN','cellQ','geneQ']]
        controls=base+['composition_prediction_A','composition_prediction_B']
        held=both[both.dataset.eq('GSE305078')]
        d=held[np.isfinite(held[['score_A','score_B']+controls].to_numpy(float)).all(axis=1)]
        x=d.score_A.to_numpy();y=d.score_B.to_numpy();C=d[controls].to_numpy(float)
        fit,status=independent_fit(x,y,C);assert status==r.status
        assert len(d)==r.complete_units and len(held)-len(d)==r.unavailable_baseline_units
        equal(r.unadjusted_rho,rho(x,y))
        original,_=independent_fit(x,y,d[base].to_numpy(float))
        equal(r.original_adjusted_same_cases,original['effect'] if original else np.nan)
        equal(r.adjusted_for_frozen_node_composition,fit['effect'] if fit else np.nan)
        if fit:
            z=np.load(output/('private_slot_%02d_resampling.npz'%r.slot),allow_pickle=False)
            assert z['samples'].tolist()==d['sample'].tolist()
            assert np.allclose(z['x'],x) and np.allclose(z['y'],y) and np.allclose(z['controls'],C)
            rng=np.random.default_rng(np.random.SeedSequence([2026101809,r.slot,0]))
            expected=np.stack([rng.permutation(len(d)) for _ in range(19999)])
            assert np.array_equal(expected,z['permutation_indices'])
            perm=np.empty(19999)
            for start in range(0,19999,1000):
                ys=fit['fit_y'][:,None]+fit['ry'][expected[start:start+1000]].T
                residual=ys-fit['D']@np.linalg.lstsq(fit['D'],ys,rcond=None)[0]
                perm[start:start+residual.shape[1]]=fit['rx']@residual/(np.linalg.norm(fit['rx'])*np.linalg.norm(residual,axis=0))
            error=float(np.max(np.abs(perm-z['permutation_effects'])));assert error<1e-9
            maximum_error=max(maximum_error,error);permutation_values+=len(perm)
            equal(r.permutation_p,(1+np.sum(np.abs(perm)>=abs(fit['effect'])-1e-12))/20000)
            rng=np.random.default_rng(np.random.SeedSequence([2026101809,r.slot,1]))
            takes=rng.integers(0,len(d),size=(2000,len(d)),dtype=np.int16)
            assert np.array_equal(takes,z['bootstrap_indices'])
            boot=[];reasons=[]
            for take in takes:
                f,reason=independent_fit(x[take],y[take],C[take]);boot.append(f['effect'] if f else np.nan);reasons.append(reason)
            boot=np.array(boot);assert np.allclose(boot,z['bootstrap_effects'],atol=1e-9,rtol=0,equal_nan=True)
            assert reasons==z['bootstrap_reasons'].tolist();valid=np.isfinite(boot)
            assert valid.sum()==r.valid_bootstrap;ci=np.quantile(boot[valid],[.025,.975])
            equal(r.lower,ci[0]);equal(r.upper,ci[1]);bootstrap_values+=len(boot)
        for s in sensitivities[sensitivities.slot.eq(r.slot)].itertuples(index=False):
            fields=(['residual_expression_A','residual_expression_B']+base if s.mode=='subtracted_node_expression'
                else ['score_A','score_B']+base+['full_composition_prediction_A','full_composition_prediction_B'])
            z=held[np.isfinite(held[fields].to_numpy(float)).all(axis=1)]
            f,status=independent_fit(z[fields[0]].to_numpy(),z[fields[1]].to_numpy(),z[fields[2:]].to_numpy(float))
            assert status==s.status and len(z)==s.complete_units;equal(s.point_effect,f['effect'] if f else np.nan)
        for s in descriptions[descriptions.slot.eq(r.slot)].itertuples(index=False):
            z=both[both.dataset.eq(s.dataset)];z=z[np.isfinite(z[['score_A','score_B']+controls].to_numpy(float)).all(axis=1)]
            f,status=independent_fit(z.score_A.to_numpy(),z.score_B.to_numpy(),z[controls].to_numpy(float))
            assert status==s.status and len(z)==s.complete_units;equal(s.point_effect,f['effect'] if f else np.nan)
    p=main.permutation_p.fillna(1).to_numpy()
    for col,by in [('BH6',False),('BY6',True)]:
        expected=family_q(p,by=by)
        for i,r in enumerate(main.itertuples(index=False)):equal(getattr(r,col),expected[i] if np.isfinite(r.permutation_p) else np.nan)
    (audit/'private_terminal.json').write_text(json.dumps(dict(
        status='ACCEPTED_ALL_NODE_BASELINES_DECOMPOSITIONS_AND_COMPLETE_FOLLOWUP_RANDOM_STREAMS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        producer_terminal_sha256=sha(output/'private_terminal.json'),node_parameter_checks=parameter_checks,
        case_view_decompositions=396,followup_tests=6,point_sensitivities=12,source_descriptions=18,
        permutation_values_replayed=permutation_values,bootstrap_values_replayed=bootstrap_values,
        max_permutation_error=maximum_error,scientific_or_clinical_validation=False),indent=2)+'\n')
    print('ALL_NODE_COMPOSITION_AND_COMPLETE_RANDOM_STREAMS_INDEPENDENTLY_ACCEPTED')


if __name__=='__main__':main()
