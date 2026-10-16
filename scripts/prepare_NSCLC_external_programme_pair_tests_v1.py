"""Freeze all source/context programme pair inputs before association effects.

Private inputs stay on the server. No RNA feature selection or correlations.
"""
import hashlib
import itertools
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_external_lineage_inputs_v2 import load
from NSCLC_native_MFA_interpretation_inputs_v1 import sha, VIEWS

CONTEXT = ['arm', 'dataset', 'cancer_label', 'sample_class', 'material_class', 'phase_key']
BLOCK = ['tissue', 'treatment_arm', 'capture_scope', 'ex_vivo']
AUDIT = 'd745f978950b0882f5ff30f3c975cd6c1e626296b2a55e4ab4464b081d44d936'


def phase_key(row):
    # Preserve W2/W4/M3, progression versus recurrence and vaccination wording.
    text = str(row.timepoint).strip().casefold().replace('-', '_').replace(' ', '_')
    if text in {'', 'unknown', 'not_reviewed', 'nan', 'none'}:
        return 'unknown_recorded_time'
    if text in {'pre', 'pre_treatment', 'pretreatment', 'baseline'}:
        return 'pre_treatment'
    if text in {'post', 'post_treatment', 'posttreatment'}:
        return 'post_treatment'
    if text in {'on_trt', 'on_treatment', 'ontreatment'}:
        return 'on_treatment'
    return text


def packed(values):
    return json.dumps([str(x) for x in values], ensure_ascii=False, separators=(',', ':'))


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 64 * 1024**2)
    src = root / 'NSCLC_external_fixed_gene_programme_assay_20261016_v1'
    proof = src / 'private_independent_fixed_programme_assay_audit.json'
    assert sha(proof) == AUDIT
    assert json.loads(proof.read_text())['producer_terminal_sha256'] == sha(src / 'private_terminal.json')
    terminal = json.loads((src / 'private_terminal.json').read_text())
    for name, item in terminal['files'].items():
        assert sha(src / name) == item['sha256']
    out = root / 'NSCLC_external_programme_pair_inputs_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__),
        assay_audit_sha256=AUDIT, budget=budget, all_factor_pair_slots=True,
        context=CONTEXT, nuisance_exchangeability_block=BLOCK,
        primary='All ten frozen factors x six unordered pairs x every source/cancer/sample/material/phase context. Same-library co-measurements only. Raw fixed programme indices; partial rank correlation controlling both own five QC vectors and all categorical block indicators.',
        analysis_unit='Reviewed patient averaged across the same paired libraries in this context; otherwise original library fallback. Mixed categorical block labels within a patient/context exclude that unit with explicit reason. No cross-library filling of a missing endpoint.',
        eligibility='Both pools n>=20, positive measured RNA and all five QC finite, both source gene-universe gates pass. >=20 analysis units in categorical blocks of >=2, n-rank(full nuisance span)>=8, nonconstant individual residuals. Full nuisance span, no automatic covariate deletion.',
        uncertainty='19999 Freedman-Lane permutations of reduced-model residuals within categorical blocks; reproject after permutation; approximate conditional exchangeability. 2000 within-block analysis-unit bootstrap with ranks and nuisance refit. All-unit and reviewed-only families separate; unknown samples never called independent patients.',
        sensitivity='Unadjusted ranks; block-only; precomputed source-QC index; add both within-library whole-annotated-scope pool fractions; all leave-one-analysis-unit effects. No new inferential primary selected from sensitivities.',
        multiplicity='BH across all estimable all-unit slots, reviewed-only BH separate; BY and simultaneous 99% binomial Monte-Carlo-upper BH sensitivity. Unestimable slots retained with missing P and reason, not treated as negative. Family frozen before effects.',
        min_units=20, min_residual_df=8, permutations=19999, bootstrap=2000,
        exact_original_factor_replication=False, direct_LR_or_clinical_validation=False,
        no_new_expression_download_or_outcome_loading=True,
        overlap='Different sources are consistency checks, not untouched independent validation; accession overlap and unreviewed patient duplication remain limitations. In-model source GSE243013, if present, labelled and excluded from external inferential family.',
        method_reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC4010955/')
    (out / 'plan_before_any_external_pair_association.json').write_text(json.dumps(plan, indent=2) + '\n')
    units = pd.read_csv(src / 'private_all6760_programme_units_metadata_quality_TMM.tsv.gz', sep='\t', keep_default_na=False)
    assert units.unit_index.tolist() == list(range(6760))
    assert pd.api.types.is_bool_dtype(units.reviewed_patient)
    units['phase_key'] = [phase_key(r) for r in units.itertuples()]
    units['block_key'] = [packed([getattr(r, c) for c in BLOCK]) for r in units.itertuples()]
    with np.load(src / 'private_all6760_fixed_programme_indices.npz', allow_pickle=False) as z:
        scores = {k: z[k] for k in z.files}
    lib = load()['lib']
    assert lib.library_id.tolist() == units.iloc[::4].library_id.tolist()
    nlib = len(lib)
    assert nlib == 1690
    denominator = lib.annotated_scope_cells.to_numpy(float)
    assert (denominator > 0).all()
    numeric = ['n_cells', 'full_source_RNA_sum', 'common_native_RNA_sum']
    for name in numeric:
        units[name] = pd.to_numeric(units[name], errors='raise')
    measurable = units.RNA_measurement_available.to_numpy(bool)
    eligibility = measurable[:, None] & scores['source_gene_gate']
    factors = scores['factors'].astype(str)
    definitions = []
    positions = []
    pair_rows = []
    rejected = []
    payload = []
    # Pre-effect selected raw audit slots use hashes, never effect ranking.
    for levels, context in units.groupby(CONTEXT, sort=True, dropna=False):
        context_id = 'PC_' + hashlib.sha256(packed(levels).encode()).hexdigest()[:16]
        meta = dict(zip(CONTEXT, [str(v) for v in levels]))
        libraries = np.unique(context.unit_index.to_numpy(int) // 4)
        assert len(context) == len(libraries) * 4
        for va, vb in itertools.combinations(range(4), 2):
            ia, ib = libraries * 4 + va, libraries * 4 + vb
            for k, factor in enumerate(factors):
                slot = '|'.join([context_id, str(factor), VIEWS[va], VIEWS[vb]])
                paired = eligibility[ia, k] & eligibility[ib, k]
                ids = libraries[paired]
                members = units.iloc[ids * 4].copy()
                start = len(payload)
                rejected_count = 0
                for block, g in members.groupby('resampling_block', sort=True, dropna=False):
                    li = g.unit_index.to_numpy(int) // 4
                    a, b = li * 4 + va, li * 4 + vb
                    assert g.reviewed_patient.nunique() == 1
                    if g.block_key.nunique() != 1:
                        rejected.append(dict(slot=slot, resampling_block=str(block), reason='PATIENT_CONTEXT_HAS_MULTIPLE_TISSUE_TREATMENT_CAPTURE_EXVIVO_BLOCKS', paired_libraries=len(li)))
                        rejected_count += 1
                        continue
                    row_id = len(payload)
                    payload.append(dict(
                        row_id=row_id, slot=slot, resampling_block=str(block),
                        reviewed_patient=bool(g.reviewed_patient.iloc[0]),
                        block_key=g.block_key.iloc[0], paired_libraries=len(li),
                        raw_A=float(scores['raw'][a, k].mean()), raw_B=float(scores['raw'][b, k].mean()),
                        source_QC_A=float(scores['source_QC'][a, k].mean()) if np.isfinite(scores['source_QC'][a, k]).all() else np.nan,
                        source_QC_B=float(scores['source_QC'][b, k].mean()) if np.isfinite(scores['source_QC'][b, k]).all() else np.nan,
                        fraction_A=float(np.mean(units.n_cells.to_numpy()[a] / denominator[li])),
                        fraction_B=float(np.mean(units.n_cells.to_numpy()[b] / denominator[li])),
                        **{'Q_A_' + str(j): float(scores['quality'][a, j].mean()) for j in range(5)},
                        **{'Q_B_' + str(j): float(scores['quality'][b, j].mean()) for j in range(5)}))
                    pair_rows.extend(dict(row_id=row_id, slot=slot, unit_A=int(x), unit_B=int(y), library_index=int(l)) for x, y, l in zip(a, b, li))
                stop = len(payload)
                definitions.append(dict(
                    slot=slot, context_id=context_id, **meta, factor=str(factor),
                    view_A=VIEWS[va], view_B=VIEWS[vb], context_libraries=len(libraries),
                    both_RNA_measurable_libraries=int((measurable[ia] & measurable[ib]).sum()),
                    both_gene_gate_libraries=int(paired.sum()), mixed_block_units_excluded=rejected_count,
                    all_units=stop-start, reviewed_units=sum(p['reviewed_patient'] for p in payload[start:stop]),
                    row_start=start, row_stop=stop, in_model_source=meta['dataset']=='GSE243013'))
                positions.append(slot)
    table = pd.DataFrame(definitions)
    assert table.slot.is_unique and table.groupby('context_id').size().eq(60).all()
    pd.DataFrame(payload).to_csv(out / 'private_all_pair_analysis_units.tsv.gz', sep='\t', index=False, na_rep='NA', compression='gzip')
    pd.DataFrame(pair_rows).to_csv(out / 'private_same_library_pair_membership.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(rejected, columns=['slot', 'resampling_block', 'reason', 'paired_libraries']).to_csv(out / 'private_excluded_mixed_patient_contexts.tsv.gz', sep='\t', index=False, compression='gzip')
    table.to_csv(out / 'private_all_frozen_pair_slots_and_coverage.tsv', sep='\t', index=False)
    terminal = dict(status='COMPLETE_ALL_EXTERNAL_PAIR_INPUTS_AND_FAMILY_FROZEN_NO_ASSOCIATIONS',
                    finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__),
                    context_count=int(table.context_id.nunique()), slots=len(table),
                    files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
                    scientific_effects_or_P_computed=False)
    assert sum(x['bytes'] for x in terminal['files'].values()) < budget['output_cap_bytes']
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'], finished_UTC=terminal['finished_UTC'], terminal_sha256=sha(out/'private_terminal.json'))))


if __name__ == '__main__':
    main()
