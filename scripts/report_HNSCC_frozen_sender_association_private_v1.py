"""Render the complete audited association family privately, retaining nulls."""
import hashlib
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


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    result = root/'HNSCC_frozen_sender_association_results_20261016_v2'
    preflight = root/'HNSCC_frozen_sender_association_preflight_20261016_v1'
    deadline = time.monotonic()+3600
    while not (result/'private_independent_audit.json').exists():
        receipt = json.loads((root/'HNSCC_receiver_followup_job_control_20261016_v2/audit_HNSCC_frozen_sender_association_v2/launch_receipt.json').read_text())
        process = Path('/proc')/str(receipt['pid'])/'stat'
        assert process.exists(), 'Original numerical audit handle missing before acceptance'
        state = process.read_text().rsplit(')', 1)[1].split()
        assert state[19] == receipt['proc_start_ticks'] and state[0] != 'Z', 'Numerical audit stopped before acceptance'
        assert time.monotonic() < deadline, 'Audit observation window expired; preserve original audit'
        time.sleep(10)
    audit = json.loads((result/'private_independent_audit.json').read_text())
    assert audit['status'].startswith('PASS_ALL400_V2_FAMILY_FULL_QR')
    stats_path = result/'all400_private_primary_association_statistics.tsv'
    desc_path = result/'all800_private_control_and_interface_descriptors.tsv'
    assert sha(stats_path) == audit['result_sha256'] and sha(desc_path) == audit['descriptor_sha256']
    stats = pd.read_csv(stats_path, sep='\t'); desc = pd.read_csv(desc_path, sep='\t')
    out = root/'HNSCC_frozen_sender_private_evidence_20261016_v2'; assert not out.exists(); out.mkdir()
    text = '''# Frozen HNSCC sender / receiver programme association (private)

All forty frozen ligand/sign candidates and ten capture-eligible sender groups
remain in the four-hundred-slot primary family. The receiver and programme
follow a prior internal clinical lead; this is exploratory and source-limited.
It is not an accepted old172-node identity bridge, independent clinical
validation, ligand-specific causal perturbation, or protein binding evidence.

Mean native sender ligand expression is related to the original full frozen
receiver programme score. All original full controls remain: both RNA and
gene-complexity measures, both native counts, sender retained-TME fraction,
receiver within-Myeloid fraction, receiver own ligand and same-cell receptor
RNA readout, treatment and HPV. Continuous variables use standardized ranks.
Beta is a conditional standardized-rank regression slope, not a correlation;
it may exceed one in magnitude. Required controls were never auto-deleted.

Inference uses twenty-thousand treatment-by-HPV-restricted residual
permutations and two-thousand patient bootstraps. Every rank transform and
complete model is recomputed for bootstrap, with at least eighteen-hundred
valid draws required for CI. This conditional exchangeability model is an
assumption, not proof of causal identification. Unestimable planned slots
remain missing and retain internal P1 in the full primary BH family.

Every eligible slot includes complete-control and interface-disjoint
complete-control point/LOO descriptors. Interface removal subtracts the exact
original frozen weighted gene contribution of ligand and all qualified
receptor components from the original sampled receiver score; remaining genes
are not refit. Descriptors receive no new P family and cannot replace the
original full-score primary result. Missing LOO due to insufficient residual
df is not silently repaired. A missing CI is not counted as a null or positive.

The original v1 primary slopes/t/P/BH are preserved exactly within roundtrip
precision. Coefficient-bootstrap and descriptor eligibility are separated from
the nonzero error variance needed only for a studentized t test; original
coverage/rank/condition/exposure guards, all bootstrap index draws and the
eighteen-hundred-valid-draw CI gate are unchanged. The failed first audit and
original results remain. Independent full QR uses orthogonal residuals to
avoid amplifying roundoff by subtracting an inverse-derived fitted value.

All native gene contributions, every primary statistic, each randomization,
bootstrap validity and CI, BH, and all descriptor/LOO records were independently
recomputed using manual ranks and full QR, distinct from production FWL.

'''
    manifest = []
    sections = [('All400 primary statistics', stats_path),
        ('All800 complete-control and interface-disjoint descriptors', desc_path),
        ('All400 prespecified design coverage and exclusion reasons', preflight/'all400_private_design_identifiability.tsv'),
        ('All400 native ligand detection coverage', root/'HNSCC_all_sender_ligand_measurements_20261016_v1/all400_private_sender_ligand_coverage.tsv')]
    for title, path in sections:
        table = pd.read_csv(path, sep='\t')
        text += '\n## '+title+'\n\n```text\n'+table.to_string(index=False)+'\n```\n'
        manifest.append(dict(path=str(path.relative_to(root)), sha256=sha(path), rows=len(table)))
    for path in [result/'private_terminal.json', result/'private_independent_audit.json', result/'contract_before_bootstrap_revision.json',
                 root/'HNSCC_frozen_sender_association_results_20261016_v1/contract_before_effect_tests.json',
                 preflight/'contract_before_any_association_effects.json', preflight/'private_independent_audit.json']:
        text += '\n## '+str(path.relative_to(root))+'\n\n```json\n'+path.read_text()+'\n```\n'
        manifest.append(dict(path=str(path.relative_to(root)), sha256=sha(path)))
    estimated = stats[stats.p.notna()].sort_values(['sender_group', 'ligand']).reset_index(drop=True)
    if len(estimated):
        fig, ax = plt.subplots(figsize=(12, max(8, .19*len(estimated))), layout='constrained')
        y = np.arange(len(estimated)); valid = np.isfinite(estimated.lower) & np.isfinite(estimated.upper)
        for i in np.flatnonzero(valid):
            r = estimated.iloc[i]
            ax.plot([r.lower, r.upper], [i, i], color='0.6', lw=.7)
        ax.scatter(estimated.beta, y, c=np.where(estimated.BH.lt(.05), '#b2182b', '#2166ac'), s=12)
        if (~valid).any():
            ax.scatter(estimated.loc[~valid, 'beta'], y[~valid], marker='x', color='black', s=16)
        ax.axvline(0, color='0.2', ls='--', lw=.7)
        ax.set_yticks(y, estimated.sender_group+' / '+estimated.ligand, fontsize=6)
        ax.invert_yaxis()
        ax.set_xlabel('Conditional standardized-rank slope; bars = patient-bootstrap 95% CI')
        ax.set_title('Complete estimable primary family; red = BH < .05, x = unavailable CI\nExploratory source-specific associations; no direct binding or causal verification', fontsize=10)
        figure = out/'PRIVATE_ALL_ESTIMABLE_PRIMARY_FOREST.png'; fig.savefig(figure, dpi=140); plt.close(fig)
        manifest.append(dict(path=str(figure.relative_to(root)), sha256=sha(figure)))
    text += '''
## Decision boundary

Keep the entire family. An uncorrected trend, unsigned prior-target rank or
RNA-compatible receptor does not establish an axis. If primary BH does not
retain a candidate, stop tuning this family and continue the existing original
CM DIALOGUE / independent function / eligible spatial directions. If a primary
candidate remains, assess CI, patient influence and interface-disjoint score
before planning independent receiver-function or cohort validation. Neither
path replaces the biological/clinical interaction goal with software success.
'''
    report = out/'HNSCC_FROZEN_SENDER_PRIVATE_ASSOCIATION_EVIDENCE.md'; report.write_text(text)
    terminal = dict(status='ALL400_PRIMARY_AND800_INTERFACE_DESCRIPTORS_AUDITED_PRIVATE_REPORT_COMPLETE',
        finished_UTC=datetime.now(timezone.utc).isoformat(), planned_slots=len(stats), estimated=len(estimated),
        primary_BH_below05=int(stats.BH.lt(.05).sum()), report_sha256=sha(report), inputs=manifest,
        direct_interaction_verified=False, independent_validation=False, goal_completed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')
    print(json.dumps(terminal))


if __name__ == '__main__':
    main()
