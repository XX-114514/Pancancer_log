"""All source calibration aggregates and rule replay; no human-array export."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os

import numpy as np
import pandas as pd


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root=Path(os.environ['PANCANCER_RAM_OUTPUT']);out=Path(os.environ['PANCANCER_RAM_AUDIT_OUTPUT']);producer=Path(os.environ['PANCANCER_MIXTURE_PRODUCER'])
    assert out.parent==Path('/dev/shm') and not out.exists()
    t=json.loads((root/'private_terminal.json').read_text());assert sha(producer)==t['source_SHA256']
    for name,v in t['files'].items():assert sha(root/name)==v['SHA256']
    contexts=json.loads((root/'private_all8_mixture_contexts.json').read_text())
    m=pd.DataFrame(json.loads((root/'private_all_source_target_recovery_metrics.json').read_text()))
    a=pd.DataFrame(json.loads((root/'private_all_source_target_applicability.json').read_text()))
    b=pd.DataFrame(json.loads((root/'private_all_source_target_binned_calibration.json').read_text()))
    p=pd.DataFrame(json.loads((root/'private_full_NNLS_numeric_proofs.json').read_text()))
    coverage=json.loads((root/'private_all_source_class_profile_coverage.json').read_text())
    plan=json.loads((root/'private_plan_before_expression_and_mixture_effects.json').read_text())['plan']
    assert len(contexts)==8 and len({(x['source'],x['scaling']) for x in contexts})==8
    assert not m.duplicated(['source','scaling','regime','target']).any()
    assert not a.duplicated(['source','scaling','target']).any() and not b.duplicated(['source','scaling','target','bin']).any()
    targets=plan['reference_classes']+list(plan['marginal_targets'])+['paired_R_only_E_minimum']
    assert len(targets)==17 and len(coverage)==40
    prior_root=Path(os.environ['PANCANCER_PRIOR_REFERENCE_OUTPUT'])
    assert sha(prior_root/'private_terminal.json')==plan['prior_reference_terminal_SHA256']
    prior_contexts=json.loads((prior_root/'private_all8_reference_contexts.json').read_text())
    for c in coverage:
        task='T_capture_state' if c['candidate'] in plan['prior_reference_plan']['T_classes'] else 'coarse_lineage'
        original=next(x for x in prior_contexts if x['task']==task and x['heldout_source']==c['source'])
        support=next(x for x in original['class_support'] if x['candidate']==c['candidate'])
        assert c['selected_cells']==support['heldout_cells'] and c['material_units']==support['heldout_material_units']
        assert c['source_profile_eligible']==(c['selected_cells']>=20 and c['material_units']>=2)
    lookup=m.set_index(['source','scaling','regime','target'])
    checked=0
    for x in contexts:
        if x['status']!='COMPLETE_KNOWN_RNA_MIXTURE_CALIBRATION':continue
        assert x['full_basis_classes']==10 and x['genes']==25188 and x['basis_condition']<=1e4
        proof=p[p.source.eq(x['source'])&p.scaling.eq(x['scaling'])];assert len(proof)==1
        q=proof.iloc[0];assert q.total_fits==30+610+x['library_stress_mixtures']
        assert q.full_objective_KKT_max<1e-7 and q.direct_coefficient_max_error<1e-6 and q.direct_full_feature_checks>=15
        for regime,n in [('dense',200),('sparse',200),('T_dilution',200),('pure_profile_anchors',10),('combined_mixed',600),('library_stress',x['library_stress_mixtures'])]:
            for target in targets:
                row=lookup.loc[(x['source'],x['scaling'],regime,target)];assert row.mixtures==n
                assert 0<=row.MAE<=row.RMSE+1e-12<=1+1e-12 and abs(row.bias)<=row.MAE+1e-12
                assert abs(row.estimate_mean-row.true_mean-row.bias)<1e-12
                assert 0<=row.positive_RNA_ge001_mixtures<=n and 0<=row.designed_zero_mixtures<=n
                if regime=='library_stress':
                    assert row.material_units==x['library_stress_material_units']
                    assert 0<=row.positive_RNA_ge001_material_units<=row.material_units
                if pd.notna(row.Pearson):assert -1-1e-12<=row.Pearson<=1+1e-12
                if pd.notna(row.Spearman):assert -1-1e-12<=row.Spearman<=1+1e-12
                if regime=='pure_profile_anchors':
                    expected=.1 if target in plan['reference_classes'] else len(plan['marginal_targets'][target])/10 if target in plan['marginal_targets'] else 0
                    assert abs(row.true_mean-expected)<1e-12
                    if target=='paired_R_only_E_minimum':assert pd.isna(row.Pearson) and row.designed_zero_mixtures==10
                checked+=1
        for target in targets:
            combined=lookup.loc[(x['source'],x['scaling'],'combined_mixed',target)]
            parts=[lookup.loc[(x['source'],x['scaling'],r,target)] for r in ['dense','sparse','T_dilution']]
            for key in ['true_mean','estimate_mean','MAE','bias']:
                assert abs(float(combined[key])-float(np.mean([row[key] for row in parts])))<1e-12
            assert abs(combined.RMSE**2-np.mean([row.RMSE**2 for row in parts]))<1e-12
            bins=b[b.source.eq(x['source'])&b.scaling.eq(x['scaling'])&b.target.eq(target)]
            assert len(bins)==10 and bins.mixtures.sum()==600
            take=bins.mixtures.gt(0)
            assert abs(np.average(bins.loc[take,'true_mean'],weights=bins.loc[take,'mixtures'])-combined.true_mean)<1e-12
            assert abs(np.average(bins.loc[take,'estimate_mean'],weights=bins.loc[take,'mixtures'])-combined.estimate_mean)<1e-12
            zero=lookup.loc[(x['source'],x['scaling'],'pure_profile_anchors',target)]
            eligible=combined.positive_RNA_ge001_mixtures>=30 and combined.MAE<=.05 and pd.notna(combined.Pearson) and combined.Pearson>=.7 and pd.notna(zero.designed_zero_max_contribution) and zero.designed_zero_max_contribution<=.1
            record=a[a.source.eq(x['source'])&a.scaling.eq(x['scaling'])&a.target.eq(target)];assert len(record)==1
            assert record.iloc[0].status==('EXPLORATORY_MEASUREMENT_CANDIDATE' if eligible else 'RECOVERY_LIMITATION')
    for (source,regime,target),g in m.groupby(['source','regime','target']):
        assert len(g)==2
        for key in ['mixtures','true_mean','positive_RNA_ge001_mixtures','designed_zero_mixtures']:
            assert abs(float(g.iloc[0][key])-float(g.iloc[1][key]))<1e-12
    spec=importlib.util.spec_from_file_location('frozen_mixture_producer',producer);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    truth=np.array([0,0,.2,.5,1.]);estimate=np.array([.1,.05,.3,.4,.9]);s=module.summarize(truth,estimate)
    assert abs(s['MAE']-.09)<1e-12 and abs(s['bias']-.01)<1e-12
    assert s['designed_zero_mixtures']==2 and s['designed_zero_leakage_gt005']==.5 and s['designed_zero_max_contribution']==.1
    rt=pd.Series(truth).rank(method='average').to_numpy();re=pd.Series(estimate).rank(method='average').to_numpy();assert abs(np.corrcoef(rt,re)[0,1]-s['Spearman'])<1e-12
    assert module.summarize(np.zeros(4),np.zeros(4))['Pearson'] is None
    out.mkdir(mode=0o700)
    result=dict(status='ACCEPTED_ALL_SOURCE_TARGET_REGIME_CONSERVATION_RULES_MATCHED_TRUTHS_BINS_AND_NUMERIC_PROOFS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),parent_terminal_SHA256=sha(root/'private_terminal.json'),registered_contexts=8,recovery_metric_positions_checked=checked,original40_class_support_exact=True,all_applicability_and_bin_positions_replayed=True,paired_scaling_truth_summaries_identical=True,synthetic_false_contribution_bias_rank_ties_and_zero_variance_checks=True,full_real_RNA_or_all_random_stream_independent_rerun=False,individual_rows_or_profiles_exported=False)
    (out/'private_source_aggregate_acceptance.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
