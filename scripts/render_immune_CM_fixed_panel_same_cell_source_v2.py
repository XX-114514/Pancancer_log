"""Render audited, private source-level CM02 functional-RNA descriptions."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget = resource(root, 8 * 1024**2)
    derived = root / 'immune_CM_fixed_panel_same_cell_descriptions_20261018_v1'
    audit = root / 'immune_CM_fixed_panel_same_cell_descriptions_audit_20261018_v1'
    proof = json.loads((audit / 'private_terminal.json').read_text())
    assert proof['status'].startswith('ACCEPTED_ALL_UNIT_PANEL_PAIR_SCOPE_NODE_')
    assert sha(derived / 'private_terminal.json') == proof['derived_terminal_sha256']
    for folder in [derived, audit]:
        for n, v in json.loads((folder / 'private_terminal.json').read_text())['files'].items():
            assert sha(folder / n) == v['sha256']
    d = pd.read_csv(derived / 'private_all_scope_source_node_supported_capture.tsv', sep='\t')
    p = pd.read_csv(derived / 'private_all_scope_source_node_panel_and_pair_descriptions.tsv', sep='\t')
    d = d[d.original_node.eq('ALL_ORIGINAL_T_NODES') & d.identity_gate.eq('pan_T_RNA_ge2') & d.valid_complete_both_panels_cells.gt(0)]
    p = p[p.original_node.eq('ALL_ORIGINAL_T_NODES') & p.identity_gate.eq('pan_T_RNA_ge2') & p.panelA.eq('regulatory_associated') & p.panelB.eq('cytotoxic')]
    d = d.merge(p[['dataset', 'scope', 'threshold', 'positive_cells', 'ratio_observed_to_depth_expected']], on=['dataset', 'scope', 'threshold'], validate='one_to_one')
    d['FOXP3_supported_fraction_of_original_R'] = d.FOXP3_supported / d.regulatory_ge2
    d['effector_supported_fraction_of_original_C'] = d.effector_supported / d.cytotoxic_ge2
    d['same_cell_effector_fraction_of_FOXP3_supported'] = d.FOXP3_and_effector_supported / d.FOXP3_supported
    d['CD8_pair_no_effector_fraction_of_original_C'] = d.CD8_pair_without_any_effector_gene / d.cytotoxic_ge2
    out = root / 'immune_CM_fixed_panel_same_cell_source_report_20261018_v2'
    assert not out.exists();out.mkdir(mode=0o700)
    d.to_csv(out / 'private_all_source_all_scope_summary.tsv', sep='\t', index=False)
    x = d[d.scope.eq('original_CM_scope_libraries')]
    sources = sorted(x.dataset.unique());positions = np.arange(len(sources))
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.6))
    fields = ['FOXP3_supported_fraction_of_original_R', 'effector_supported_fraction_of_original_C',
              'same_cell_effector_fraction_of_FOXP3_supported', 'ratio_observed_to_depth_expected']
    titles = ['FOXP3-supported capture within original R panel',
              'Effector-supported capture within original C panel',
              'Same-cell effector capture within FOXP3-supported R',
              'Original R-C overlap / within-unit RNA-depth reference']
    colours = ['#1674a3', '#d77b28']
    for i, (ax, field, title) in enumerate(zip(axes.flat, fields, titles)):
        for j, threshold in enumerate([1, 2]):
            values = x[x.threshold.eq(threshold)].set_index('dataset').reindex(sources)[field].to_numpy(float)
            displayed = values * 100 if i < 3 else values
            bars = ax.bar(positions + (j - .5) * .36, displayed, .34, color=colours[j], label='Per-gene count >=%d' % threshold)
            for bar, value in zip(bars, displayed):
                if np.isfinite(value):ax.annotate(('%.1f%%' if i < 3 else '%.2f') % value,
                    (bar.get_x() + bar.get_width() / 2, value), xytext=(0, 3), textcoords='offset points', ha='center', fontsize=8, bbox=dict(facecolor='white', edgecolor='none', pad=.1))
        ax.set_xticks(positions, sources, fontsize=9);ax.set_title(title, fontsize=10)
        ax.spines[['top', 'right']].set_visible(False)
        ax.set_ylim(0, 110 if i < 3 else 1.15)
        ax.set_ylabel('Captured-cell fraction (%)' if i < 3 else 'Descriptive observed / expected')
        if i == 3:ax.axhline(1, color='#656565', ls='--', lw=1)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.5, .93), ncol=2, frameon=False)
    fig.suptitle('2026-10-14 immune atlas CM02 (v1): fixed functional proxy', y=.985, fontsize=13)
    fig.text(.02, .105, 'R: >=2 original regulatory-associated genes. FOXP3 support: FOXP3 plus another original R gene.\n'
             'C: >=2 original cytotoxic panel genes. Effector support: >=2 of GZMB / PRF1 / CTSW / NKG7.\n'
             'All cells pass >=2 pan-T genes at count>=1. Source aggregates; stages include pre / post / unknown.', fontsize=9)
    fig.text(.02, .025, 'Descriptive RNA capture, not validated Treg identity, functional killing, patient replication or direct signalling.\n'
             'Non-detection is not biological absence. No new P values or clinical refitting.', fontsize=9)
    fig.tight_layout(rect=[0, .20, 1, .89])
    fig.savefig(out / 'private_CM02_source_capture_interpretation.png', dpi=180)
    fig.savefig(out / 'private_CM02_source_capture_interpretation.pdf')
    plt.close(fig)
    (out / 'private_plot_scope_and_limits.json').write_text(json.dumps(dict(
        recorded_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), budget=budget,
        derived_acceptance_sha256=sha(audit / 'private_terminal.json'),
        CM='2026-10-14 immune atlas CM v1 T_whole_TME_primary_tumor CM02 fixed functional proxy',
        displayed_sources=sources, all_sources_and_scopes_kept_in_table=True,
        display_scope='Original CM material, both fixed thresholds, pan-T RNA gate; no source chosen by effect.',
        interpretation='Fractions describe captured RNA support for existing panels. Neither lineage/function nor bulk score source is proven by RNA co-detection; original min-percentile clinical score is not decomposed by these fractions.',
        actual_visual_review='PENDING'), indent=2) + '\n')
    files = [p for p in out.iterdir() if p.is_file()]
    terminal = dict(status='COMPLETE_PRIVATE_AUDITED_SOURCE_FIGURE_VISUAL_REVIEW_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in files}, no_new_statistical_tests=True)
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps(dict(status=terminal['status'], completed_UTC=terminal['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
