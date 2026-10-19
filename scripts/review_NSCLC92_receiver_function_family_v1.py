"""Source-only complete function plots, no molecular endpoints or individual rows."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC92_receiver_function_math_v2 import QUESTIONS, PARTITIONS


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 4 * 1024**2)
    out = root / 'NSCLC92_receiver_function_review_20261019_v1'
    assert not out.exists(); out.mkdir(mode=0o700)
    src = root / 'NSCLC92_receiver_function_family_20261019_v1'
    proof = root / 'NSCLC92_receiver_function_family_audit_20261019_v1/private_terminal.json'
    at = json.loads(proof.read_text())
    assert at['status'].startswith('ACCEPTED_COMPLETE_FUNCTION_FAMILY_ALL_MAIN')
    assert at['producer_terminal_sha256'] == sha(src / 'private_terminal.json')
    t = json.loads((src / 'private_terminal.json').read_text())
    for name, d in t['files'].items(): assert sha(src / name) == d['sha256']
    a = dict(np.load(root / 'NSCLC92_receiver_function_inputs_20261019_v2/private_complete_function_inputs.npz'))
    r = pd.read_csv(src / 'private_complete_function_family_results.tsv', sep='\t')
    labels = [str(v).removeprefix('HALLMARK_').replace('_', ' ') for v in a['references']]
    titles = dict(clinical_response='Function difference: pCR/MPR minus non-MPR', axis_given_response='Function association with axis score, response controlled', sender_given_receiver_and_response='Function association with sender RNA, receiver RNA and response controlled')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    for scope in a['scope_ids']:
        z = r[r.scope_id.eq(scope)]
        for question in QUESTIONS:
            frame = z[z.question.eq(question)]
            fig, axes = plt.subplots(1, 2, figsize=(16, 15), sharey=True)
            fig.subplots_adjust(left=.29, right=.97, top=.91, bottom=.085, wspace=.18)
            fig.suptitle('GSE243013: complete exploratory receiver-function family\n' + titles[question] + '\nOriginal partitions; all references retained; pointwise WCU-B intervals; joint family BH/BY.', fontsize=12, y=.985)
            finite = frame[['CI_low', 'CI_high']].to_numpy(); finite = finite[np.isfinite(finite)]
            limit = max(1., np.max(np.abs(finite)) * 1.08) if finite.size else 1.
            for ax, partition, label in zip(axes, PARTITIONS, ['Original 129-record partition', 'Original 87-record partition']):
                sub = frame[frame.partition.eq(partition)].set_index('reference').reindex(a['references'])
                ax.set_title(label); ax.set_ylim(49.7, -.8); ax.set_xlim(-limit, limit); ax.axvline(0, color='black', lw=.7); ax.grid(axis='x', alpha=.15)
                ax.set_yticks(np.arange(50)); ax.set_yticklabels(labels); ax.tick_params(axis='y', labelsize=8)
                for y, row in enumerate(sub.itertuples()):
                    if np.isfinite(row.effect):
                        colour = '#a12b37' if row.BY_complete_family < .05 else '#266ab2' if row.BH_complete_family < .05 else '#777777'
                        ax.errorbar(row.effect, y, xerr=[[row.effect - row.CI_low], [row.CI_high - row.effect]], fmt='o', color=colour, ms=3, lw=.7, capsize=1)
                    else: ax.text(.04, y, 'Not estimable', transform=ax.get_yaxis_transform(), va='center', fontsize=7, color='#777777')
                ax.set_xlabel('Function score SD per specified target contrast')
            fig.legend(handles=[Line2D([0], [0], color='#a12b37', marker='o', label='Joint BH and BY < .05'), Line2D([0], [0], color='#266ab2', marker='o', label='Joint BH only < .05'), Line2D([0], [0], color='#777777', marker='o', label='No joint corrected support')], loc='lower center', ncol=3, frameon=False)
            name = '%s_%s_complete50' % (scope, question)
            fig.savefig(out / (name + '.png'), dpi=140); fig.savefig(out / (name + '.pdf')); plt.close(fig)
        fig, ax = plt.subplots(figsize=(13, 15)); fig.subplots_adjust(left=.35, right=.95, top=.92, bottom=.09)
        columns = [(p, q) for p in PARTITIONS for q in QUESTIONS]
        counts = np.zeros((50, 6)); eligible = np.zeros((50, 6), bool)
        for c, (p, q) in enumerate(columns):
            sub = z[z.partition.eq(p) & z.question.eq(q)].set_index('reference').reindex(a['references'])
            counts[:, c] = sub.author_records.to_numpy(); eligible[:, c] = sub.status.eq('ESTIMABLE').to_numpy()
        cmap = plt.get_cmap('Blues').copy(); cmap.set_bad('#eeeeee')
        ax.imshow(np.where(eligible, counts, np.nan), aspect='auto', cmap=cmap, vmin=0, vmax=129)
        ax.set_yticks(np.arange(50)); ax.set_yticklabels(labels, fontsize=8); ax.set_xticks(np.arange(6)); ax.set_xticklabels(['129 clinical', '129 axis', '129 sender', '87 clinical', '87 axis', '87 sender'], rotation=30, ha='right')
        for y in range(50):
            for c in range(6): ax.text(c, y, str(int(counts[y, c])), ha='center', va='center', fontsize=7, color='white' if eligible[y, c] and counts[y, c] > 85 else 'black')
        fig.suptitle('GSE243013: all 50 references and complete design coverage\nCounts are author records, not verified independent patients. Gray means ineligible function or design.', fontsize=11)
        fig.savefig(out / (str(scope) + '_complete50_coverage.png'), dpi=140); fig.savefig(out / (str(scope) + '_complete50_coverage.pdf')); plt.close(fig)
    files = {p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}
    assert sum(d['bytes'] for d in files.values()) < 4 * 1024**2
    receipt = dict(status='COMPLETE_FULL_SOURCE_AGGREGATE_FUNCTION_PLOTS_ACTUAL_VISUAL_INSPECTION_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), result_audit_sha256=sha(proof), resource_guard=guard, files=files, learned_molecular_endpoint_names_exported=False, gene_weights_exported=False, individual_records_exported=False)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
