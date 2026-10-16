"""Wait for complete audit, then review every factor, pair and source context."""
import itertools
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha, VIEWS
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_NSCLC_external_programme_coordination_v1 import read_inputs, QCOLS
from NSCLC_external_programme_coordination_math_v1 import prepare, effect


def wait_owned(root, task, required):
    deadline = time.monotonic()+10800
    receipt = json.loads((root/'HNSCC_receiver_followup_job_control_20261016_v2'/task/'launch_receipt.json').read_text())
    while not required.exists():
        proc = Path('/proc')/str(receipt['pid'])/'stat'
        assert proc.exists(), 'Bound audit missing without proof; inspect original log'
        fields = proc.read_text().rsplit(')', 1)[1].split()
        assert fields[0] != 'Z' and fields[19] == str(receipt['proc_start_ticks'])
        assert time.monotonic() < deadline, 'Observation expired; do not infer termination'
        time.sleep(3)


def material_group(text):
    if text in {'RECORDED_PRIMARY_TUMOR', 'RECORDED_SOLID_TUMOR', 'RECORDED_SOLID_METASTASIS', 'TUMOR_CAPTURE_SITE_UNSPECIFIED'}:
        return 'SOLID_TUMOR_CONTEXTS'
    return str(text)


def secondary(a, b, quality, blocks):
    try:
        return effect(prepare(a, b, quality, blocks)), 'ESTIMABLE'
    except ValueError as e:
        return np.nan, str(e)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 32*1024**2)
    out = root/'NSCLC_external_programme_coordination_review_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), budget=budget,
                scope='All630 estimated slots and all15600 eligibility positions; no primary endpoint, gene, factor or sign changes. Review only after complete real-result audit.',
                diagnostics='On the exact same active cases, add both sex-linked RNA proxies, both stress readouts, or both panels to the full primary nuisance span; no P or altered primary family. Missing panels mean unavailable, no complete-case rescue.',
                source_selection='For each family/factor/pair/material group, exactly one effect-blind context per GSE: phase priority pre_treatment, unknown_recorded_time, on_treatment, post_treatment, other literal phases; then largest active units, largest reviewed units, lexical context and arm. Actual source-arm/phase rows all retained.',
                consistency='Different phases and arms of one GSE count at most once. All-unit and reviewed families are not independent evidence. Descriptive Fisher-z mean weights=max(n-quality_rank-3,1), source range and leave-one-source effects only; no pooled P/CI or untouched-independent validation claim.',
                grading='At least2 distinct reviewed-family GSE with same-direction BH<.05 and bootstrapCI excluding0 and no opposite supported source: cross-source consistency candidate. One source: source-limited candidate. Opposite supported sources: heterogeneous. All-unit-only counterparts remain library-level leads. No direct LR or clinical validation.',
                attenuation='Per-result descriptive flag if absolute adjusted rho <half primary or sign reverses; not a new statistical significance test or automatic artifact diagnosis.',
                recorded_cancer_labels='Literal labels retained; distinct strings not asserted to be different cancers. Blood, products, normal and other materials stay separate. Solid-primary/metastatic/site-unspecified rows remain explicit in selected context.',
                coordination_results_loaded=False)
    (out/'plan_before_any_external_result_review.json').write_text(json.dumps(plan, indent=2)+'\n')
    src = root/'NSCLC_external_programme_coordination_20261016_v1'
    proof = src/'private_independent_coordination_audit.json'
    wait_owned(root, 'audit_NSCLC_external_programme_coordination_v1', proof)
    audit = json.loads(proof.read_text())
    assert audit['status'].startswith('PASS_ALL15600_ELIGIBILITY_ALL_REAL_')
    assert audit['producer_terminal_sha256'] == sha(src/'private_terminal.json')
    diagnostic = root/'NSCLC_external_public_RNA_diagnostics_20261016_v1'
    diag_proof = diagnostic/'private_independent_public_RNA_diagnostic_audit.json'
    wait_owned(root, 'audit_NSCLC_external_public_diagnostic_panels_v1', diag_proof)
    assert json.loads(diag_proof.read_text())['status'].startswith('PASS_ALL61_PANEL_COVERAGE_')
    resource(root, 32*1024**2)
    source, _, _, slots, cases = read_inputs(root)
    result = pd.read_csv(src/'private_all_estimated_programme_coordination.tsv', sep='\t')
    assert len(result) == audit['estimated_tests_checked']
    members = pd.read_csv(source/'private_same_library_pair_membership.tsv.gz', sep='\t')
    with np.load(diagnostic/'private_all6760_public_RNA_diagnostics.npz', allow_pickle=False) as z:
        values = z['values']
    # Exact matched-library aggregation for every frozen pair analysis row.
    panel = np.full((len(cases), 4), np.nan)
    rid = members.row_id.to_numpy(int)
    for side, offset in [('A', 0), ('B', 2)]:
        ids = members['unit_'+side].to_numpy(int)
        number = np.bincount(rid, minlength=len(cases))
        for k in range(2):
            x = values[ids, k]
            finite = np.bincount(rid, weights=np.isfinite(x).astype(int), minlength=len(cases))
            summed = np.bincount(rid, weights=np.nan_to_num(x), minlength=len(cases))
            good = (number > 0) & (finite == number)
            panel[good, offset+k] = summed[good]/number[good]
    additions = []
    scope = slots.set_index('slot')
    for row in result.itertuples():
        base = scope.loc[row.slot]
        g = cases.iloc[int(base.row_start):int(base.row_stop)]
        if row.family == 'reviewed_patients':
            g = g[g.reviewed_patient]
        p = prepare(g.raw_A.to_numpy(), g.raw_B.to_numpy(), g[QCOLS].to_numpy(), g.block_key.to_numpy(str))
        g = g.iloc[np.flatnonzero(p['keep'])]
        assert len(g) == row.active_units
        extra = panel[g.row_id.to_numpy(int)]
        item = dict(record=row.record, active_units=len(g))
        for name, columns in [('RNA_sex', [0, 2]), ('IEG_stress', [1, 3]), ('both_panels', [0, 1, 2, 3])]:
            q = np.column_stack([p['quality'], extra[:, columns]])
            rho, status = secondary(p['a'], p['b'], q, p['blocks'])
            item[name+'_rho'] = rho
            item[name+'_status'] = status
            item[name+'_complete_case_units'] = int(np.isfinite(extra[:, columns]).all(axis=1).sum())
            item[name+'_attenuation_or_reversal'] = bool(np.isfinite(rho) and (abs(rho)<.5*abs(row.partial_rank_rho) or rho*row.partial_rank_rho<0))
        additions.append(item)
    sensitivity = pd.DataFrame(additions)
    sensitivity.to_csv(out/'private_all_fixed_RNA_background_sensitivities.tsv', sep='\t', index=False, na_rep='NA')
    full = result.merge(sensitivity.drop(columns='active_units'), on='record', validate='one_to_one')
    full['material_group'] = full.material_class.map(material_group)
    phase_order = {'pre_treatment': 0, 'unknown_recorded_time': 1, 'on_treatment': 2, 'post_treatment': 3}
    full['phase_priority'] = full.phase_key.map(phase_order).fillna(4)
    keys = ['family', 'factor', 'view_A', 'view_B', 'material_group', 'dataset']
    selected = full.sort_values(['phase_priority', 'phase_key', 'active_units', 'reviewed_units', 'context_id', 'arm'], ascending=[True, True, False, False, True, True]).drop_duplicates(keys)
    selected['BH_CI_support'] = selected.BH_q.lt(.05) & ((selected.CI_low.gt(0) & selected.partial_rank_rho.gt(0)) | (selected.CI_high.lt(0) & selected.partial_rank_rho.lt(0)))
    selected.to_csv(out/'private_effect_blind_one_context_per_source.tsv', sep='\t', index=False, na_rep='NA')
    consistency = []
    for levels, g in selected.groupby(keys[:-1], sort=True):
        family, factor, va, vb, mat = levels
        positive = g[g.BH_CI_support & g.partial_rank_rho.gt(0)]
        negative = g[g.BH_CI_support & g.partial_rank_rho.lt(0)]
        weights = np.maximum(g.active_units.to_numpy()-g.quality_rank.to_numpy()-3, 1)
        zz = np.arctanh(np.clip(g.partial_rank_rho.to_numpy(), -1+1e-12, 1-1e-12))
        mean = float(np.tanh(np.average(zz, weights=weights)))
        leave = [float(np.tanh(np.average(np.delete(zz, j), weights=np.delete(weights, j)))) for j in range(len(g))] if len(g)>1 else []
        if len(positive) and len(negative):
            grade = 'DIRECTION_HETEROGENEOUS'
        elif max(len(positive), len(negative)) >= 2:
            grade = 'CROSS_SOURCE_CONSISTENCY_CANDIDATE' if family == 'reviewed_patients' else 'LIBRARY_LEVEL_CROSS_SOURCE_LEAD_PATIENT_REVIEW_REQUIRED'
        elif max(len(positive), len(negative)) == 1:
            grade = 'SOURCE_LIMITED_CANDIDATE' if family == 'reviewed_patients' else 'LIBRARY_LEVEL_SOURCE_LEAD_PATIENT_REVIEW_REQUIRED'
        else:
            grade = 'INSUFFICIENT_STATISTICAL_SUPPORT'
        consistency.append(dict(family=family, factor=factor, view_A=va, view_B=vb, material_group=mat,
                                unique_sources=len(g), recorded_cancer_label_strings=g.cancer_label.nunique(),
                                positive_source_effects=int(g.partial_rank_rho.gt(0).sum()), negative_source_effects=int(g.partial_rank_rho.lt(0).sum()),
                                BH_CI_supported_positive_sources=len(positive), BH_CI_supported_negative_sources=len(negative),
                                descriptive_weighted_rho=mean, minimum_source_rho=g.partial_rank_rho.min(), maximum_source_rho=g.partial_rank_rho.max(),
                                leave_one_source_min=min(leave) if leave else np.nan, leave_one_source_max=max(leave) if leave else np.nan,
                                grade=grade, clinical_or_direct_LR_validation=False))
    pd.DataFrame(consistency).to_csv(out/'private_all_programme_pair_material_source_consistency.tsv', sep='\t', index=False, na_rep='NA')
    # Every factor, all six pairs, all eligible contexts, no result-selected plots.
    pairs = list(itertools.combinations(VIEWS, 2))
    short = {VIEWS[0]:'CD8', VIEWS[1]:'Treg', VIEWS[2]:'Mac', VIEWS[3]:'B/Plasma'}
    figures = []
    for factor in sorted(full.factor.unique()):
        fig, axes = plt.subplots(2, 3, figsize=(17, 11), constrained_layout=True)
        for ax, (va, vb) in zip(axes.flat, pairs):
            g = full[(full.factor == factor) & (full.family == 'all_units') & (full.view_A == va) & (full.view_B == vb)].sort_values(['material_group', 'dataset', 'phase_key'])
            for j, row in enumerate(g.itertuples()):
                known = row.active_units == row.reviewed_units
                color = '#2166ac' if known else '#888888'
                if np.isfinite(row.CI_low) and np.isfinite(row.CI_high):
                    ax.plot([row.CI_low, row.CI_high], [j, j], color=color, linewidth=1.3)
                supported = row.BH_q < .05 and (row.CI_low > 0 or row.CI_high < 0)
                ax.scatter([row.partial_rank_rho], [j], color=color, marker='*' if supported else 'o', s=75 if supported else 24, zorder=3)
            labels = [f'{r.dataset} | {r.phase_key} | {r.material_group} | n={r.active_units}' for r in g.itertuples()]
            ax.set_yticks(np.arange(len(g)), labels, fontsize=6)
            ax.invert_yaxis()
            ax.axvline(0, color='black', linewidth=.7)
            ax.set_xlim(-1.05, 1.05)
            ax.set_title(short[va]+' / '+short[vb], fontsize=11)
            ax.set_xlabel('Partial rank rho (own/partner QC + category blocks)', fontsize=8)
        fig.suptitle(factor+' | fixed gene programme co-variation; not direct communication\nBlue: all units reviewed patients; grey: library fallback present; star: all-unit BH<0.05 and CI excludes zero', fontsize=12)
        name='private_'+factor+'_all_source_context_forests.png'
        fig.savefig(out/name, dpi=140)
        plt.close(fig)
        figures.append(name)
    full.to_csv(out/'private_all_primary_and_sensitivity_review.tsv', sep='\t', index=False, na_rep='NA')
    report = dict(status='COMPLETE_ALL_EXTERNAL_PROGRAMME_COORDINATION_CONTEXT_REVIEW_AND_FIXED_BACKGROUND_DESCRIPTIONS_VISUAL_REVIEW_PENDING',
                  finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                  coordination_audit_sha256=sha(proof), diagnostic_audit_sha256=sha(diag_proof), figures=figures,
                  files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
                  direct_communication_or_clinical_validation=False)
    assert sum(v['bytes'] for v in report['files'].values()) < budget['output_cap_bytes']
    (out/'private_terminal.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status=report['status'], finished_UTC=report['finished_UTC'], terminal_sha256=sha(out/'private_terminal.json'))))


if __name__ == '__main__':
    main()
