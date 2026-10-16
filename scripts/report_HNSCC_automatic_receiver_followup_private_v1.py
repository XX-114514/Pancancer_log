"""Combine audited prior, signed attribution and all-sender evidence privately."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    directories = {
        'prior': root/'HNSCC_receiver_prior_ligand_ranking_20261016_v1',
        'receptor': root/'HNSCC_prior_receptor_qualification_20261016_v1',
        'components': root/'HNSCC_frozen_receiver_component_attribution_20261016_v1',
        'sender': root/'HNSCC_all_sender_ligand_measurements_20261016_v1'}
    audits = {}
    for key, path in directories.items():
        a = json.loads((path/'private_independent_audit.json').read_text())
        assert a['status'].startswith('PASS_')
        audits[key] = a
    assert sha(directories['prior']/'private_all_ligand_prior_rankings.tsv') == audits['prior']['prior_ranking_sha256']
    assert sha(directories['receptor']/'private_all_frozen_ligand_qualification.tsv') == audits['receptor']['qualified_table_sha256']
    assert sha(directories['components']/'private_patient_signed_and_quality_components.tsv') == audits['components']['component_table_sha256']
    assert sha(directories['sender']/'all400_private_sender_ligand_coverage.tsv') == audits['sender']['coverage_sha256']
    out = root/'HNSCC_automatic_receiver_private_evidence_20261016_v1'
    assert not out.exists()
    out.mkdir()
    text = '''# Automatic receiver-prior and native sender follow-up (private)

## Question and provenance

The receiver follows an already observed HNSCC clinical programme lead. This
is exploratory follow-up, not independent clinical validation. The source is
the 2026-10-14 immune atlas v1, not a relabeling of the original v20261010
nine CM or an accepted mapping between their nodes. All native data and the
original programme weights, means, SDs and quality correction remain fixed.

## What the ranking measures

Every prior ligand was ranked against positive and negative original gene
sets separately, using a native expressed prior-target background. The score
is custom average precision above foreground prevalence, with AUROC/Pearson
alongside; it is not official NicheNet activity, a P value, receptor-specific
causality or the predicted sign of a ligand effect. The prior is unsigned.
The forty prospective sign/ligand slots were frozen before sender expression
effects. Failed receptor entries were retained, never replaced by rank21+.

Reference constituent support and same-cell receptor RNA do not establish
protein-complex assembly, direct binding, or current-tumour signaling. Patient
coverage is included because pooled detection can hide sparse support.

## Signed and technical attribution

The full programme equals positive-weight centred gene contribution plus
negative-weight centred gene contribution plus the frozen quality correction.
These components are not mutually independent and can each take either sign.
Exact reconstruction and every clinical/functional descriptor and leave-one
patient calculation were independently replayed. No additional P/BH family
was computed, and a favorable component does not replace the original full
programme endpoint. The gene-disjoint Hallmarks exclude both learned endpoint
signatures from the signature and ranking background. Descriptive correlations
do not establish IFN causality or formal cell states.

## Sender evidence

All frozen ligands were measured in every prospectively capture-eligible
major sender using all native matched baseline cells, full-RNA normalization,
and unique measured source aliases. C1QC receivers are excluded from Myeloid
senders. Missing genes and absent cell groups are missing expression, measured
zeros remain zeros, and no sample is counted as an independent extra patient.
Reported any/1%/10% detection levels describe coverage; no favorable threshold
was selected and no sender-response effect was tested in this report.

'''
    sections = [
        ('Complete prior ligand rankings', 'prior', 'private_all_ligand_prior_rankings.tsv'),
        ('All frozen ligand receptor qualifications', 'receptor', 'private_all_frozen_ligand_qualification.tsv'),
        ('All literal receptor reference and joint RNA edges', 'receptor', 'private_all_shortlist_receptor_edges.tsv'),
        ('Receptor per-patient coverage aggregate', 'receptor', 'private_all_edge_patient_RNA_coverage.tsv'),
        ('All sender ligand coverage slots', 'sender', 'all400_private_sender_ligand_coverage.tsv'),
        ('Complete signed component clinical descriptors', 'components', 'all45_private_component_clinical_descriptors.tsv'),
        ('Complete signed component independent functional descriptors', 'components', 'all1500_private_component_function_descriptors.tsv'),
        ('Component covariance decomposition', 'components', 'private_component_covariance_decomposition.tsv')]
    manifest = []
    for title, key, filename in sections:
        path = directories[key]/filename
        table = pd.read_csv(path, sep='\t')
        text += '\n## '+title+'\n\n```text\n'+table.to_string(index=False)+'\n```\n'
        manifest.append(dict(path=str(path.relative_to(root)), sha256=sha(path), rows=len(table)))
    for key, path in directories.items():
        for filename in ['private_terminal.json', 'private_independent_audit.json']:
            p = path/filename
            text += '\n## '+key+' '+filename+'\n\n```json\n'+p.read_text()+'\n```\n'
            manifest.append(dict(path=str(p.relative_to(root)), sha256=sha(p)))
    frame = pd.read_csv(directories['sender']/'all400_private_sender_ligand_coverage.tsv', sep='\t')
    prior = pd.read_csv(directories['prior']/'private_all_ligand_prior_rankings.tsv', sep='\t')
    selected = prior[prior.prospective_top20].sort_values(['programme_sign', 'all_prior_rank'])
    sender_order = sorted(frame.sender_group.unique())
    ligand_order = selected.ligand.tolist()
    pct = frame.pivot(index='ligand', columns='sender_group', values='pooled_sender_positive_fraction').reindex(index=ligand_order, columns=sender_order)
    n = frame.pivot(index='ligand', columns='sender_group', values='patients_with_any_positive_ligand_cells').reindex(index=ligand_order, columns=sender_order)
    qualified = frame.groupby('ligand').receiver_prior_qualified.first()
    labels = [gene+' ['+sign+']'+(' *' if qualified[gene] else '') for gene, sign in zip(selected.ligand, selected.programme_sign)]
    fig, axes = plt.subplots(1, 2, figsize=(15, 15), layout='constrained')
    for ax, data, title, vmax in zip(axes, [pct, n], ['Fraction of native sender cells detecting ligand', 'Matched patients with any ligand-positive sender cell'], [1, 26]):
        im = ax.imshow(data.to_numpy(float), vmin=0, vmax=vmax, aspect='auto', cmap='viridis')
        ax.set_yticks(np.arange(len(labels)), labels, fontsize=7)
        ax.set_xticks(np.arange(len(sender_order)), sender_order, rotation=60, ha='right', fontsize=8)
        ax.set_title(title, fontsize=11)
        fig.colorbar(im, ax=ax, shrink=.65)
    fig.suptitle('Frozen HNSCC ligand coverage; * = constituent reference + pooled receptor RNA only', fontsize=12)
    figure = out/'PRIVATE_ALL400_SENDER_LIGAND_COVERAGE.png'
    fig.savefig(figure, dpi=160)
    plt.close(fig)
    manifest.append(dict(path=str(figure.relative_to(root)), sha256=sha(figure)))
    text += '''
## Next decision and non-claims

The next stage must use actual sender measurements and a prospectively frozen
patient-level model with required abundance, native quality, own-ligand and
receptor controls. Check sample coverage, variable ranks and residual degrees
of freedom first. Do not automatically drop required covariates or recode
missing to negative to obtain an estimable association. Keep the full planned
family, negative results and source restriction. Interface-disjoint programme
sensitivity and independent cohorts remain separate requirements.

This report does not verify a direct interaction, clinical prediction,
resistance programme, complete receptor assembly or formal cell state. The
original clinical and negative communication families remain intact. The full
biological/clinical interaction goal is active.
'''
    path = out/'HNSCC_AUTOMATIC_RECEIVER_PRIVATE_EVIDENCE.md'
    path.write_text(text)
    terminal = dict(status='AUDITED_PRIOR_RECEPTOR_COMPONENT_AND_NATIVE_SENDER_EVIDENCE_COMBINED_PRIVATE',
        finished_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(Path(__file__)), inputs=manifest,
        report_sha256=sha(path), new_effect_tests=False, direct_interaction_verified=False, goal_completed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')
    print(json.dumps(terminal))


if __name__ == '__main__':
    main()
