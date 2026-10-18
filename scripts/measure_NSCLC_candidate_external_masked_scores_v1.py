"""Original LR-excluded frozen receiver indices for every fixed external record.

Reuse accepted counts and TMM effective libraries. No new inference or refit.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import load,sha,VIEWS
from describe_NSCLC_native_interface_mechanistic_followup_v1 import read_scope,SCOPE_SHA
from prepare_NSCLC_fixed_programme_continuous_priors_v1 import read_rds,PRIORS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,16*1024**2)
    scope,candidates=read_scope(dict(root=root));data=load()
    assay=root/'NSCLC_external_fixed_gene_programme_assay_20261016_v1'
    proof=assay/'private_independent_fixed_programme_assay_audit.json'
    assert sha(proof)=='d745f978950b0882f5ff30f3c975cd6c1e626296b2a55e4ab4464b081d44d936'
    assert json.loads(proof.read_text())['producer_terminal_sha256']==sha(assay/'private_terminal.json')
    at=json.loads((assay/'private_terminal.json').read_text())
    for name in ['private_all6760_fixed_programme_indices.npz','private_all6760_programme_units_metadata_quality_TMM.tsv.gz']:
        assert sha(assay/name)==at['files'][name]['sha256']
    with np.load(assay/'private_all6760_fixed_programme_indices.npz',allow_pickle=False) as f:
        full=f['raw'];whole_gate=f['source_gene_gate']
    units=pd.read_csv(assay/'private_all6760_programme_units_metadata_quality_TMM.tsv.gz',sep='\t',keep_default_na=False)
    assert len(units)==6760 and units.unit_index.tolist()==list(range(6760))
    assert units['view'].tolist()==list(np.tile(VIEWS,1690))
    cal=root/'NSCLC_disjoint_receiver_universe_calibration_20261018_v2'
    proof2=cal/'private_independent_disjoint_universe_audit.json'
    assert sha(proof2)=='52dfdcf10bf4873a2ca83f9da0b8a044f4b087a7c6b75735e0691c50b8f717b9'
    assert json.loads(proof2.read_text())['producer_terminal_sha256']==sha(cal/'private_terminal.json')
    ct=json.loads((cal/'private_terminal.json').read_text());name='private_all97600_source_disjoint_receiver_calibrations.tsv.gz'
    assert sha(cal/name)==ct['files'][name]['sha256']
    gates=pd.read_csv(cal/name,sep='\t').set_index(['slot_index','arm','dataset'])
    assert gates.index.is_unique
    counts_root=root/'NSCLC_native_MFA_external_frozen_role_count_inputs_20261016_v2'
    count_proof=counts_root/'private_independent_native_count_input_audit.json'
    assert sha(count_proof)=='7bd028cc1c89ac7324ac25c65042b5b184cb377162a38f9544caf7bb24048c04'
    assert json.loads(count_proof.read_text())['producer_terminal_sha256']==sha(counts_root/'private_terminal.json')
    count_terminal=json.loads((counts_root/'private_terminal.json').read_text())
    name='private_all61_count_archive_manifest.tsv'
    assert sha(counts_root/name)==count_terminal['files'][name]['sha256']
    archives=pd.read_csv(counts_root/name,sep='\t')
    prior=root/'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    path=prior/'lr_network_human_21122021.rds';assert sha(path)==PRIORS[path.name]
    flat=read_rds(path)
    path=root/'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv'
    assert sha(path)=='df574533d866438e646423f29bd555cb283173dbceadfe9f8264e1622a86fd06'
    reference=pd.read_csv(path,sep='\t',keep_default_na=False)
    omitted={}
    for ligand in set(candidates.ligand):
        strings={ligand,*flat.loc[flat['from'].eq(ligand),'to'],*reference.loc[reference.ligand.eq(ligand),'receptor']}
        omitted[ligand]=set(strings)|{g for s in strings for g in str(s).split('_')}
    out=root/'NSCLC_candidate_external_original_masked_scores_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),scope_sha256=SCOPE_SHA,
        budget=budget,assay_audit_sha256=sha(proof),disjoint_universe_audit_sha256=sha(proof2),count_audit_sha256=sha(count_proof),
        candidates_sha256=scope['candidate_table_sha256'],model_sha256=sha(data['model_path']),
        definition='Same original literal ligand plus ALL original flat/reference receptor components removed. Do not silently add later curated receptor components or refit/renormalize coefficients.',
        external_assay='Original accepted common-native pseudobulk log2(1+CPM), original accepted TMM effective libraries and frozen A/RMS. Source quality belongs to later full association model. Extension assay, not exact full-author latent identity.',
        scope='Every original candidate across all1690 libraries; measurable receiver pools regardless of gene-universe gate retained with separate gate flags, empty/unavailable pools remain NaN.',
        numeric_checks='All finite receiver full indices replay original accepted assay; direct kept-gene projection equals full index minus excluded contribution AND independently augmented least-squares frozen projection. 1e-9 absolute/relative tolerance.',
        new_normalization_refit_associations_P_or_clinical_outcomes=False)
    (out/'plan_before_any_external_masked_score.json').write_text(json.dumps(plan,indent=2)+'\n')
    response=np.full((len(candidates),1690),np.nan);available=np.zeros_like(response,bool);qualified=np.zeros_like(response,bool)
    source_gate=np.zeros_like(response,bool);disjoint_gate=np.zeros_like(response,bool)
    projection=[];alternate=[]
    for vi in range(4):
        w=data['model']['weights_'+str(vi)];tau=data['model']['precision_'+str(vi)];rms=float(data['prepared'][vi]['training_view_RMS'])
        a=np.linalg.solve(np.eye(10)+w.T@(tau[:,None]*w),w.T*tau)/rms
        aug=np.vstack([np.eye(10),np.sqrt(tau)[:,None]*w]);rhs=np.vstack([np.zeros((10,2000)),np.diag(np.sqrt(tau))])
        other=np.linalg.lstsq(aug,rhs,rcond=None)[0]/rms
        assert np.allclose(a,other,atol=1e-12,rtol=1e-10)
        projection.append(a);alternate.append(other)
    maximum_full_error=maximum_subtraction_error=maximum_augmented_error=0.
    for source in archives.itertuples():
        path=counts_root/source.archive;assert sha(path)==source.sha256
        with np.load(path,allow_pickle=False) as f:counts=f['counts'];columns=f['common_native_columns'];unit=f['global_unit_indices']
        physical=np.full(31831,-1,int);physical[columns]=np.arange(len(columns))
        assert units.iloc[unit].arm.eq(source.arm).all() and units.iloc[unit].dataset.eq(source.dataset).all()
        for ri,item in enumerate(candidates.itertuples()):
            vi=VIEWS.index(item.receiver_view);k=int(item.factor.rsplit('F',1)[1])-1
            g=gates.loc[(item.slot_index,source.arm,source.dataset)]
            assert g['view']==item.receiver_view and g['factor']==item.factor and int(g.direction)==int(item.direction)
            all_ids=unit[unit%4==vi];libraries=all_ids//4
            source_gate[ri,libraries]=whole_gate[all_ids,k]
            disjoint_gate[ri,libraries]=bool(g.source_disjoint_receiver_universe_gate)
            assert np.all(whole_gate[all_ids,k]==bool(g.original_full_source_universe_gate))
            take=(unit%4==vi)&np.isfinite(full[unit,k]);ids=unit[take]
            if not len(ids):continue
            norm=pd.to_numeric(units.iloc[ids].effective_common_library_size).to_numpy(float)
            assert np.isfinite(norm).all() and (norm>0).all()
            selected=data['prepared'][vi]['selected_native_columns'];exists=physical[selected]>=0
            value=np.log2(1+counts[take][:,physical[selected[exists]]].astype(float)*1e6/norm[:,None])
            names=data['genes'][selected[exists]];keep=~np.isin(names,list(omitted[item.ligand]))
            a=projection[vi][k,exists];b=alternate[vi][k,exists]
            reconstructed=value@a;direct=value[:,keep]@a[keep]
            subtracted=full[ids,k]-value[:,~keep]@a[~keep];other=value[:,keep]@b[keep]
            for expected,observed in [(full[ids,k],reconstructed),(direct,subtracted),(direct,other)]:
                assert np.allclose(expected,observed,atol=1e-9,rtol=1e-9)
            maximum_full_error=max(maximum_full_error,float(np.max(abs(full[ids,k]-reconstructed))))
            maximum_subtraction_error=max(maximum_subtraction_error,float(np.max(abs(direct-subtracted))))
            maximum_augmented_error=max(maximum_augmented_error,float(np.max(abs(direct-other))))
            response[ri,ids//4]=direct;available[ri,ids//4]=True
    qualified=available&source_gate&disjoint_gate
    assert np.isfinite(response[available]).all() and np.isnan(response[~available]).all()
    for _,group in candidates.groupby('slot_index'):
        rows=[candidates.index.get_loc(i) for i in group.index]
        for ri in rows[1:]:assert np.allclose(response[ri],response[rows[0]],equal_nan=True)
    np.savez_compressed(out/'private_all_candidate_library_masked_scores.npz',
        records=candidates.record.to_numpy(str),slot_index=candidates.slot_index.to_numpy(int),direction=candidates.direction.to_numpy(int),
        library_index=np.arange(1690),response=response,measurement_available=available,
        original_full_source_gene_gate=source_gate,original_disjoint_source_gate=disjoint_gate,measurement_and_source_gates=qualified)
    candidates[['record','slot_index','sender_view','receiver_view','factor','direction','ligand']].to_csv(out/'private_fixed_candidate_axis.tsv',sep='\t',index=False)
    files=[p for p in out.iterdir() if p.is_file()]
    assert sum(p.stat().st_size for p in files)<16*1024**2
    terminal=dict(status='COMPLETE_ALL_FIXED_CANDIDATE_EXTERNAL_ORIGINAL_MASKED_INDICES_WITH_DUAL_FORMULA_CHECKS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        maximum_full_assay_error=maximum_full_error,maximum_subtraction_error=maximum_subtraction_error,maximum_augmented_error=maximum_augmented_error,
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},
        standalone_independent_external_validation=False,no_new_P_or_score_definition_change=True)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],completed_UTC=terminal['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
