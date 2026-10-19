"""Source-only descriptive calibration tables/figures for every fixed slot.

Reference labels are not TCR truth. No new P value, clinical model, threshold
selection, cell-level figure or formal identity adoption is performed.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from gamma_delta_RNA_rules_v1 import RULES
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); guard = resource(root, 3 * 1024**2)
    out = root / 'NSCLC_fixed_gamma_scGate_pilot_review_20261019_v1'; assert not out.exists(); out.mkdir(mode=0o700)
    audit = root / 'NSCLC_gamma_scGate_fixed_pilot_audit_20261019_v1/private_terminal.json'
    t = wait(root, out, 'audit_NSCLC_gamma_scGate_pilot_v1', audit, {'ACCEPTED_ALL_UNIFORM_PILOT_SAMPLING_ANONYMOUS_R_LABEL_TRANSFER_AND_AUTHOR_AGGREGATES_STRUCTURAL_ONLY'})
    src = root / 'NSCLC_gamma_scGate_fixed_pilot_20261019_v1'
    assert t['producer_terminal_sha256'] == sha(src / 'private_terminal.json')
    z = dict(np.load(src / 'private_complete_pilot_candidate_classification_sidecar.npz', allow_pickle=False))
    rows = z['native_row']; assert len(rows) == 20000
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    unit = np.load(native / 'private_cell_unit_index.npy', allow_pickle=False)
    u = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False).sort_values('unit_index')
    ref = u.iloc[unit[rows]].reset_index(drop=True)
    group = np.full(len(rows), 4, np.int8)
    group[ref.major_cell_type.eq('T/NK cell')] = 3
    group[ref.sub_cell_type.str.startswith('NK_')] = 2
    group[ref.sub_cell_type.eq('T_gdT_TRDV2')] = 1
    group[ref.sub_cell_type.eq('T_gdT_TRDV1')] = 0
    categories = ['Author Vdelta1', 'Author Vdelta2', 'Author NK', 'Other T/NK', 'Non-T/NK']
    totals = np.bincount(group, minlength=5)
    values = pd.DataFrame(z['classes'], columns=z['models'].tolist())
    hard = dict(np.load(root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1/private_candidate_identity_sidecar.npz', allow_pickle=False))
    fullbits = np.zeros(len(unit), np.uint8); fullbits[hard['native_row']] = hard['rule_bits']; bits = fullbits[rows]
    masks = {RULES[j]: (bits & (1 << j)) > 0 for j in range(4)}
    masks.update({name: values[name].eq('Pure').to_numpy() for name in ['PBMC_exact', 'HiTME_native_subset', 'NK_HiTME_exact']})
    masks['HiTME_full'] = None
    both = masks['PBMC_exact'] & masks['HiTME_native_subset']; anygamma = masks['PBMC_exact'] | masks['HiTME_native_subset']; nk = masks['NK_HiTME_exact']
    tiers = dict(Both_no_NK=both & ~nk, HiTME_only_no_NK=masks['HiTME_native_subset'] & ~masks['PBMC_exact'] & ~nk, PBMC_only_no_NK=masks['PBMC_exact'] & ~masks['HiTME_native_subset'] & ~nk, Gamma_NK_conflict=anygamma & nk, No_gamma_call=~anygamma)
    assert np.array_equal(sum(v.astype(int) for v in tiers.values()), np.ones(len(rows), int))
    masks['Both_no_NK'] = tiers['Both_no_NK']
    order = RULES + ['PBMC_exact', 'HiTME_native_subset', 'Both_no_NK', 'NK_HiTME_exact', 'HiTME_full']
    labels = ['Raw core5 count>=1', 'Raw supported count>=1', 'Raw core5 count>=2', 'Raw supported count>=2', 'PBMC exact: gamma target', 'HiTME subset: gamma target', 'Both gamma models; NK negative', 'NK control: NK target', 'HiTME full: unmeasurable']
    records = []; count = np.full((len(order), 5), np.nan)
    for i, name in enumerate(order):
        chosen = masks[name]
        if chosen is not None: count[i] = np.bincount(group[chosen], minlength=5)
        for j, category in enumerate(categories):
            calls = int(chosen.sum()) if chosen is not None else None
            amount = int(count[i, j]) if chosen is not None else None
            records.append(dict(source='GSE243013', rule_or_model=name, reference_group=category, pilot_reference_cells=int(totals[j]), called_cells=amount, all_model_calls=calls, reference_group_coverage=amount / int(totals[j]) if amount is not None and totals[j] else None, fraction_among_calls=amount / calls if calls else None, status='NATIVE_GENE_UNMEASURABLE' if chosen is None else 'DESCRIPTIVE_REFERENCE_COMPARISON_NOT_TCR_ACCURACY', target='NK' if name == 'NK_HiTME_exact' else 'gamma_RNA_candidate'))
    table = pd.DataFrame(records); table.to_csv(out / 'private_complete_source_reference_calibration.tsv', sep='\t', index=False)
    joint = []
    for name, selected in tiers.items():
        n = np.bincount(group[selected], minlength=5)
        for j, category in enumerate(categories): joint.append(dict(source='GSE243013', tier=name, reference_group=category, pilot_reference_cells=int(totals[j]), cells=int(n[j]), fraction_of_reference_group=n[j] / totals[j] if totals[j] else None))
    joint = pd.DataFrame(joint); joint.to_csv(out / 'private_complete_source_joint_candidate_tiers.tsv', sep='\t', index=False)
    coverage = np.divide(count, totals[None, :], out=np.full_like(count, np.nan), where=totals[None, :] > 0)
    ncall = np.nansum(count, axis=1); ncall[-1] = np.nan
    composition = np.divide(count, ncall[:, None], out=np.full_like(count, np.nan), where=ncall[:, None] > 0)
    plt.rcParams.update({'font.size': 9, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, axs = plt.subplots(1, 2, figsize=(15, 7.7))
    cmap = plt.get_cmap('Blues').copy(); cmap.set_bad('#dddddd')
    for ax, a, title in zip(axs, [coverage, composition], ['Coverage of each author reference group', 'Composition of calls (row target)']):
        im = ax.imshow(100 * a, cmap=cmap, vmin=0, vmax=100, aspect='auto')
        ax.set_xticks(range(5), [c + '\n(n=%d)' % n for c, n in zip(categories, totals)], rotation=30, ha='right')
        ax.set_yticks(range(len(labels)), labels)
        ax.set_title(title)
        for i in range(len(order)):
            for j in range(5):
                if i == len(order) - 1: text = 'not run'
                elif not np.isfinite(a[i, j]): text = 'undefined'
                else: text = '%.1f%%\nn=%d' % (100 * a[i, j], count[i, j])
                ax.text(j, i, text, ha='center', va='center', fontsize=8, color='white' if np.isfinite(a[i, j]) and a[i, j] > .6 else 'black')
        fig.colorbar(im, ax=ax, shrink=.65, label='%')
    fig.suptitle('GSE243013 | fixed 20,000-cell pilot | all prespecified slots', y=.98)
    fig.text(.01, .015, 'Original author types are reference labels, not TCR ground truth. NK control predicts NK, not gamma cells.\nCell fractions are descriptive; no patient-level accuracy, clinical test, external validation or final identity is claimed.', fontsize=9)
    fig.tight_layout(rect=(0, .085, 1, .95))
    for ext in ['png', 'pdf']: fig.savefig(out / ('source_complete_reference_calibration.' + ext), dpi=180, bbox_inches='tight')
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10.5, 6.4)); bottom = np.zeros(5)
    colors = ['#2166ac', '#67a9cf', '#8073ac', '#ef8a62', '#d9d9d9']
    for (name, selected), color in zip(tiers.items(), colors):
        n = np.bincount(group[selected], minlength=5); height = np.divide(100 * n, totals, out=np.zeros(5), where=totals > 0)
        ax.bar(range(5), height, bottom=bottom, label=name, color=color); bottom += height
    ax.set_xticks(range(5), [c + '\n(n=%d)' % n for c, n in zip(categories, totals)]); ax.set_ylim(0, 100); ax.set_ylabel('% within author reference group')
    ax.set_title('GSE243013 | gamma-model agreement and NK conflict')
    ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)
    fig.text(.01, .02, 'Natural-prevalence pilot; not independent TCR truth or patient-level validation.\nModel-only and NK-conflict candidates remain separate; no expression or cell-level plot is exported.', fontsize=9)
    fig.tight_layout(rect=(0, .085, 1, .97))
    for ext in ['png', 'pdf']: fig.savefig(out / ('source_complete_joint_candidate_tiers.' + ext), dpi=180, bbox_inches='tight')
    plt.close(fig)
    receipt = dict(status='COMPLETE_ALL_SOURCE_REFERENCE_CALIBRATION_AND_JOINT_TIERS_FIGURES_VISUAL_AND_SCIENTIFIC_REVIEW_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), structural_audit_sha256=sha(audit), producer_terminal_sha256=sha(src / 'private_terminal.json'), resource_guard=guard, original_model_slots=4, compared_rules_and_candidate_tiers=len(order), no_P_or_clinical_or_identity_adoption=True, figures_actual_visual_reviewed=False, scientific_outputs_public=False, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in receipt['files'].values()) < 3 * 1024**2
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
