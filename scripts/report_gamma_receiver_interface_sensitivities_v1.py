"""Complete private follow-up synthesis and source-only forest summary."""
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
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 512 * 1024)
    request = root / 'private_gamma_receiver_interface_priority_request_20261019_v1.json'
    assert sha(request) == 'fe8587f44faadc8d0d61f376070f1875b69daecc7b9ed2f6d7f4fa52fe0045e4'
    expert = json.loads(request.read_text())
    for info in expert['inputs'].values():
        path = root / info['relative_terminal']
        assert sha(path) == info['sha256']
        for name, item in json.loads(path.read_text()).get('files', {}).items():
            assert sha(path.parent / name) == item['sha256']
    identity = pd.read_csv(root / 'gamma_receiver_identity_sensitivity_20261019_v1/private_complete7_identity_sensitivity_effects.tsv', sep='\t')
    assembly = pd.read_csv(root / 'fixed_assembly_clinical_sensitivity_20261019_v1/private_complete2_fixed_assembly_followup_effects.tsv', sep='\t')
    assert len(identity) == 7 and len(assembly) == 2
    original = pd.read_csv(root / 'NSCLC_frozen92_axis_source_consistency_20261019_v1/private_all92_source_consistency_results.tsv', sep='\t').set_index('axis_id')
    assert len(original) == 92
    axis = identity.axis_id.unique().item()
    assert original.loc[axis, 'BH92'] < .05 <= original.loc[axis, 'BY92']
    out = root / 'gamma_receiver_interface_sensitivity_review_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    combined = pd.concat([identity.assign(followup='IDENTITY'), assembly.assign(followup='ASSEMBLY')], ignore_index=True, sort=False)
    assert len(combined) == 9
    combined.to_csv(out / 'private_complete9_followup_descriptions.tsv', sep='\t', index=False)
    save('private_scientific_review_and_next_scope.json', dict(expert, original_primary_BH=float(original.loc[axis, 'BH92']), original_primary_BY=float(original.loc[axis, 'BY92']), all_nine_descriptions_preserved=True, new_inferential_test=False))
    take = lambda mask, comparison: identity[(identity['mask'] == mask) & (identity.comparison == comparison)].iloc[0]
    display = [
        ('Original reference', take('AUTHOR_REFERENCE', 'REDEFINED_RECEIVER_OWN_ELIGIBLE_RECORDS'), '#346b96'),
        ('Original reference, matched records', take('REFERENCE_WITHOUT_COMPLETE_AB_ANNOTATION', 'ORIGINAL_RECEIVER_SAME_ELIGIBLE_RECORDS'), '#777777'),
        ('Exclude complete alpha/beta annotation', take('REFERENCE_WITHOUT_COMPLETE_AB_ANNOTATION', 'REDEFINED_RECEIVER_OWN_ELIGIBLE_RECORDS'), '#c57c27'),
        ('Specified ligand assembly extension', assembly[assembly.ligand_proxy == 'SPECIFIED_ASSEMBLY_EXTENSION'].iloc[0], '#43816a'),
    ]
    # The repeated assembly baseline is identical, not an independent fifth fit.
    np.testing.assert_allclose(assembly[assembly.ligand_proxy == 'ORIGINAL_COMPONENT'][['effect', 'CV3_se']].iloc[0].to_numpy(float), display[0][1][['effect', 'CV3_se']].to_numpy(float), atol=1e-10, rtol=1e-10)
    plt.rcParams.update({'font.size': 10, 'font.family': 'DejaVu Sans', 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(9.5, 4.7))
    fig.subplots_adjust(left=.39, right=.91, top=.79, bottom=.31)
    ax.axvline(0, color='#777777', lw=.8, ls='--')
    for i, (label, row, colour) in enumerate(display):
        assert row.status == 'ESTIMABLE'
        ax.errorbar(row.effect, i, xerr=[[row.effect - row.CV3_t_CI_low], [row.CV3_t_CI_high - row.effect]], fmt='o', color=colour, ms=6, capsize=3, lw=1.4)
        ax.text(1.025, i, 'n=%d' % row.author_records, transform=ax.get_yaxis_transform(), va='center', fontsize=9)
    ax.set_yticks(range(4), [label for label, _, _ in display])
    ax.set_ylim(3.6, -.6)
    ax.set_xlim(min(float(row.CV3_t_CI_low) for _, row, _ in display) - .08, .08)
    ax.grid(axis='x', alpha=.15)
    ax.set_xlabel('Adjusted RNA-proxy difference: pCR/MPR minus non-MPR\n(original frozen axis SD)', fontsize=10)
    fig.suptitle('GSE243013 | Identity and ligand-assembly sensitivities', fontsize=13, y=.96)
    fig.text(.5, .865, 'Within-cohort, post-treatment follow-up; four nonduplicated estimable comparisons', ha='center', fontsize=10)
    fig.text(.03, .13, 'Conditional CV3 t-reference 95% intervals; post-result descriptions, no new P/FDR.\nRNA-rule-only subsets: n=17 and15, not estimable. All9 descriptions retained in the source table.\nIdentity and assembly were separate fits. No independent-patient, receptor-identity or binding validation.', fontsize=8.5, va='center')
    for extension in ['png', 'svg']:
        fig.savefig(out / ('source_identity_and_assembly_sensitivity.' + extension), dpi=160, facecolor='white')
    plt.close(fig)
    receipt = dict(status='COMPLETE_NINE_PRIVATE_DESCRIPTIONS_SCIENTIFIC_SCOPE_AND_SOURCE_FIGURE_ACTUAL_VISUAL_REVIEW_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), private_priority_request_sha256=sha(request), all_followup_positions=9, displayed_nonduplicated_estimable_positions=4, full_tables_include_unestimable_positions=True, resource_guard=guard, new_P_or_FDR=False, original_inference_or_formal_identity_changed=False, source_aggregate_only=True, scientific_outputs_public=False, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in receipt['files'].values()) < 512 * 1024
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
