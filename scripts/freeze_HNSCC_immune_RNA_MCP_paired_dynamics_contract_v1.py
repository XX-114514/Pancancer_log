"""Freeze all paired questions using identity/design only; no RNA scores/outcomes."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    qualification=root/'HNSCC_immune_RNA_MCP_paired_materials_20261019_v1'
    assert sha(qualification/'private_independent_metadata_acceptance_20261019_v1.json')=='8dc0f5b74ea8cf41ff10a62eb5c804064a0dc558ec90e787f4e9cbac5b4461c9'
    pairs=pd.read_csv(qualification/'private_all_exact_pre_post_patient_material_pairs.tsv',sep='\t');pairs=pairs[~pairs.post_used_in_original_discovery].sort_values('patient_registry')
    assert len(pairs)==26 and pairs.patient_registry.is_unique
    path=root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    assert sha(path)=='9942706870518c20c87a45fed44caff6bec22a1a512b6f194dc17d406df5cd22'
    a=pd.read_csv(path,sep='\t',keep_default_na=False,usecols=['library_id','patient_registry','timepoint_registry','neo_type','HPV/p16'])
    a['sampleID']='GSE286827::'+a.library_id.astype(str);assert a.sampleID.is_unique
    rows=[]
    for r in pairs.itertuples(index=False):
        g=a[a.patient_registry.eq(r.patient_registry)]
        pre=g[g.sampleID.eq(r.pre_sampleID)];post=g[g.sampleID.eq(r.post_sampleID)]
        assert len(pre)==len(post)==1 and pre.timepoint_registry.eq('pre').all() and post.timepoint_registry.eq('post').all()
        assert g.neo_type.nunique()==1 and g['HPV/p16'].nunique()==1
        treatment=pre.neo_type.iloc[0];hpv=pre['HPV/p16'].iloc[0]
        complete=treatment in ['D','D+T'] and hpv in ['negative','positive']
        rows.append(dict(patient_registry=r.patient_registry,pre_sampleID=r.pre_sampleID,post_sampleID=r.post_sampleID,
            treatment_DplusT=int(treatment=='D+T') if treatment in ['D','D+T'] else None,
            HPV_positive=int(hpv=='positive') if hpv in ['negative','positive'] else None,complete_static_controls=complete))
    controls=pd.DataFrame(rows)
    views=['R_FOXP3_capture','E_effector_capture','Myeloid'];programmes=['MCP1','MCP2']
    shifts=[dict(slot=i,kind='paired_median_shift',program=p,view=v) for i,(p,v) in enumerate((p,v) for p in programmes for v in views)]
    cochange=[dict(slot=i+6,kind='paired_cochange',program=p,view_A=views[j],view_B=views[k]) for i,(p,j,k) in enumerate((p,j,k) for p in programmes for j,k in [(0,1),(0,2),(1,2)])]
    out=root/'HNSCC_immune_RNA_MCP_paired_dynamics_contract_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        paired_metadata_acceptance_sha256=sha(qualification/'private_independent_metadata_acceptance_20261019_v1.json'),author_identity_design_source_sha256=sha(path),
        measurement_producer_source_sha256='21c7010b54590a7b86c3a46a6d55cd2ef7a066b6ace6887012da3c04cca4462d',measurement_auditor_source_sha256='1ebb1b8172ea491ebd15fc5673470508ddf7f076759735d53436c950815e18b2',
        prerequisite='All paired native raw/count/score/domain/official anchor independent acceptance BEFORE computing any difference/model. Wait original same source-bound audit PID, never restart because observation times out.',
        slots=shifts+cochange,patients=26,materials=52,difference='post minus pre, original38 discovery case/view/program SD(ddof0), no new scales or learned features',
        shift=dict(effect='Median paired delta',P='Two-sided exact binomial sign test, numerical zeros abs(standardized_delta)<=1e-10 omitted only from nonzero P denominator',minimum_pairs=20,minimum_nonzero=10,CI='All paired deltas including zeros; largest order k with 2*BinomCDF(k-1,n,.5)<=.05, [delta[k-1],delta[n-k]], at least95% distribution-free median interval',point_only=['mean_delta','paired_rank_biserial']),
        cochange=dict(effect='Partial Spearman correlation of both paired deltas',controls=['treatment_DplusT','HPV_positive']+[f+'_'+s for s in ['A','B'] for f in ['cellQ','geneQ','fraction','logN']],continuous_controls='post minus pre, averaged ranks; categories retain binary coding',
            minimum_pairs=20,minimum_full_residual_df=10,relative_rank_threshold=1e-12,minimum_retained_singular_ratio=1e-8,maximum_condition=1e8,maximum_leverage=.999999,expected_bootstrap_nuisance_rank='fixed original full nuisance rank',residual_norm='greater than1e-8*max(1,norm(rawranktarget))',full_nuisance_span_reconstruction_required=True,no_covariate_deletion=True,
            test='HC3 studentized restricted-null HC3 Rademacher wild19999',bootstrap='Unstratified whole patient, pre/post together;2000 draws, ranks rebuilt, at least1800 valid for95CI; every invalid reason retained',leave_one='All complete patients, full original control set/rank',seed=2026101903),
        multiplicity='All12 BH/BY, unavailable main P=1 internally with display missing; six cochange99% simultaneous ClopperPearson MC upper P plus exact signP then all12 upper BH. No family reduction.',
        finite_null_gate=dict(draws_per_error_design=256,wild_per_draw=999,error_models=['normal','unit_variance_t3','normal_times_exp(.6_standardized_first_continuous_control)'],alpha=.05,maximum_rejection_rate=.10,seed=2026101904,any_failure_stops_real_effects=True,not_universal_typeI_guarantee=True),
        cochange_point_sensitivities=['unadjusted','treatment_HPV_only','both_RNA_quality_changes_only','all_captured_mean','fixed_stress_complement','plus_both_frozen_original_node_composition_changes'],
        shift_point_sensitivities=['all_captured_mean','fixed_stress_complement'],point_sensitivities_new_P=False,
        node_composition='Frozen original discovery-node coefficients, primary sampled-view proportions; unseen node contribution missing, not zero, and sensitivity unavailable. Original coefficients not learned from post.',
        candidate_gate='Cochange main BH12 and upperBH12<.05,95CI excludes0, all LOO estimable/same sign; all sensitivities/domain/identity limits reported. Shift evidence not substitute for two-view cochange.',
        limitations='Source/internal and conditional fixed weights; baseline/patients overlap discovery; time/procedure/therapy confounded, no untreated controls, no independent clinical or direct LR mechanism validation. No new clinical outcome model or MCP2-MCP1 composite.',
        no_RNA_score_difference_or_outcome_fields_read=True)
    (out/'plan_before_any_paired_difference_or_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    controls.to_csv(out/'private_all26_static_identity_design_controls.tsv',sep='\t',index=False)
    summary=dict(dataset='GSE286827',reviewed_paired_patients=26,paired_materials=52,complete_static_controls=int(controls.complete_static_controls.sum()),fixed_primary_questions=12,classification='Internal longitudinal programme study; no effects computed')
    (out/'private_source_design_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    terminal=dict(status='FROZEN_ALL12_PAIRED_DYNAMIC_QUESTIONS_NO_SCORES_OR_EFFECTS_READ',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],terminal_sha256=sha(out/'private_terminal.json'),complete_static_controls=summary['complete_static_controls'])))


if __name__=='__main__':main()
