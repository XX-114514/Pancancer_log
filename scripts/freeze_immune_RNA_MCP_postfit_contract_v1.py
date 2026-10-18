"""Freeze output acceptance, seed matching and holdout tests without learned values."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 16 * 1024**2)
    model_plan = root / 'immune_RNA_MCP_official_models_20261018_v1/plan_before_any_official_MCP_fit.json'
    plan = json.loads(model_plan.read_text())
    assert len(plan['fits']) == 5
    output = root / 'immune_RNA_MCP_postfit_contract_20261018_v1'
    assert not output.exists()
    output.mkdir(mode=0o700)
    views = ['R_FOXP3_capture', 'E_effector_capture', 'Myeloid']
    contract = dict(
        frozen_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), budget=budget,
        official_fit_plan_sha256=sha(model_plan),
        learned_genes_scores_effects_or_clinical_outcomes_read=False,
        anchor='K2_seed20261018_source_adjusted',
        planned_fits=plan['fits'], planned_projection_slots=33,
        view_to_official_name={v: v.replace('_', '') for v in views},
        technical_acceptance=dict(
            complete_five_fit_statuses_required=True,
            no_retraining_or_parameter_repair=True,
            score_parity_absolute_tolerance=1e-6,
            unsupported='Keep missing final objects, no positive NNLS coefficients and any reconstruction failures as explicit unavailable slots; never fill scores or missing genes with zero.',
            cell_sample_gene_keys='Every discovery/heldout expression and score key, source labels and train-only means/SD/covariates must match. Independent implementation reconstructs all supported cell scores and case means.',
            old_slot_summary='The original worker convenience counters use underscored names. Official make.cell.type removes underscores; recompute all33 slots using the explicit map without changing or rerunning official models.'),
        projection=dict(
            gene_scores='Use official positive NNLS coefficients with official up/down signs. Fit gene means/SD and nuisance coefficients exclusively in discovery; reproduce every saved official discovery score.',
            unknown_source='Heldout metadata sourceCC stays NA. The discovery source coefficient is unidentifiable for a new source. Subtract only known nuisance/intercept terms; heldout primary scores retain an unknown constant per view/program. Within-source correlations/ranks are invariant to this offset; absolute cross-source levels are not comparable.',
            all_five_models='Project every supported slot, not only the most favourable seed/source mode.',
            retained_source='Entire GSE305078 source,28 eligible sample cases; patient alias unreviewed, not independent patient validation or untouched independent original-CM validation.'),
        seed_matching=dict(
            full_pair_matrix='For every view and programme pair, Spearman across discovery case means after ranking/centering/scaling within each discovery source. Match all programme pairs jointly across three views by Hungarian assignment maximizing absolute mean correlation. One global sign per matched programme; do not independently flip each view.',
            weights='Also report full signed gene-weight cosine matrix on union of candidate genes, missing coefficient zero within measured opportunity only. All candidate features are already common measured.',
            all_slots_retained=True,
            descriptive_stability_rule='Both other K2 seeds must have aligned within-source-centred case-score correlation>=0.70 and signed weight cosine>=0.50 in every view. These are prespecified engineering diagnostics, not calibrated significance/clinical thresholds. Unmatched or unsupported slots remain unstable/unavailable.',
            sensitivities='K3 and source-permissive models matched/reported descriptively; neither substitutes for the fixed primary anchor. No seed/K chosen by holdout P.'),
        primary_source_association=dict(
            slots=[dict(program=p, view_A=views[i], view_B=views[j])
                   for p in ['MCP1', 'MCP2'] for i, j in [(0, 1), (0, 2), (1, 2)]],
            family='All6 anchor K2 heldout view-pair slots, two-sided BH across6; unavailable P=1 only for family adjustment and displayed effect/P=NA. BY additionally descriptive.',
            primary='Partial correlation of ranked per-sample mean frozen scores, controlling both endpoints original full-material view fractions and log cell counts plus mean selected-cell full-RNA/detected-gene quality variables. Rank all continuous variables within heldout source; no automatic covariate deletion.',
            unadjusted='Spearman of same sample scores with same complete-case denominator.',
            prerequisites='At least20 complete sample units, both scores variable, full-rank design, condition<1e8 and residual degrees of freedom>=10. No imputation of missing measurement or unavailable NNLS score.',
            permutation='19999 Freedman-Lane reduced-model residual permutations of ranked outcome at fixed complete-case sample units; re-residualize permuted outcomes on the full nuisance design, compare absolute residual correlations; plus-one correction.',
            interval='2000 sample-block bootstrap draws; rerank and refit same complete design each draw, preserve failed/degenerate draws and require>=1000 valid draws for percentile95% CI.',
            seed=2026101803,
            reviewed_patient_sensitivity='Separate reviewed-patient-only subset, coverage gate unchanged. If none reviewed, explicitly unavailable. Main sample-case inference assumes distinct units; unverified recurrence remains a limitation.',
            discovery='Both discovery-source within-source effects and leave-one-unit ranges are descriptive fitted-data consistency checks, not heldout validation. No discovery P selects candidate slots.',
            clinical='No clinical endpoint tests are in this contract. Exact patient/outcome mapping and a separate effect-blind clinical family are required first.'),
        boundary='MCP coordination is not a direct ligand-receptor interaction, secretion, binding, immune suppression or treatment mechanism. Source/cancer-specific programmes can retain ecological value with their scope stated. All scientific scores/weights/effects remain private.')
    path = output / 'plan_before_learned_outputs_or_projection.json'
    path.write_text(json.dumps(contract, indent=2) + '\n')
    (output / 'private_terminal.json').write_text(json.dumps(dict(
        status='POSTFIT_CONTRACT_FROZEN_NO_LEARNED_OUTPUTS_READ',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        files={path.name: sha(path)}), indent=2) + '\n')
    print('POSTFIT_CONTRACT_FROZEN_NO_LEARNED_VALUES_READ')


if __name__ == '__main__':
    main()
