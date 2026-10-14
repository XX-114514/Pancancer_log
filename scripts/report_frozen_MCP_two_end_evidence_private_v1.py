"""Combine audited evidence on server; actual scientific values are runtime inputs.

No patient rows or learned signatures appear in the report/stdout. Existing
statistical families stay unchanged; this does not turn diagnostics into tests.
"""
import hashlib, json, os
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'frozen_MCP_two_end_private_evidence_20261014_v1'
    assert not out.exists()
    audit_specs = [
        ('frozen_MCP_two_end_function_review_20261014_v1/private_independent_audit.json', 'PASS_ALL1500'),
        ('HNSCC_all_T_group_sender_associations_20261014_v1/private_independent_audit.json', 'PASS_ALL576'),
        ('HNSCC_receiver_adopted_ambient_sensitivity_20261014_v1/private_independent_audit.json', 'PASS_ALL18'),
        ('full_parent_denominator_independent_audit_20261014_v1/private_independent_audit.json', 'PASS_ALL_RETAINED'),
    ]
    audits = []
    for relative, prefix in audit_specs:
        value = json.loads((root / relative).read_text())
        assert value['status'].startswith(prefix)
        audits.append((relative, value))
    tables = [
        ('Original frozen clinical family; unchanged', 'HNSCC_fixed_MCP_clinical_results_20261014_v2/all24_private_clinical_statistics.tsv'),
        ('All source/end/function diagnostics; no new significance', 'frozen_MCP_two_end_function_review_20261014_v1/all1500_private_function_diagnostics.tsv'),
        ('All T-group sender associations including diagnostic-only signatures', 'HNSCC_all_T_group_sender_associations_20261014_v1/all576_private_T_group_associations.tsv'),
        ('All adopted ambient-removal descriptors; no new clinical P/BH', 'HNSCC_receiver_adopted_ambient_sensitivity_20261014_v1/all18_private_adopted_ambient_diagnostics.tsv'),
    ]
    text = '''# Fixed two-end programme evidence (private server report)

## Questions and scope

Separate the original CM composition, frozen programme covariance, one-end
clinical association, functional interpretation and direct molecular interaction.
They are different evidence requirements. The old172 nodes/9CM and original
weights remain frozen. New-atlas CD8/Treg/naive-memory and macrophage labels
are multigene subtype candidates; this is not an accepted old-node mapping or
formal state annotation.

## Statistical and measurement contracts

All50 external functional references and all available sources/endpoints are
retained. Learned programme genes are removed from both functional signatures
and rank backgrounds. RNA quality/count and leave-one-patient diagnostics are
descriptive only. Missing coverage/df and unavailable patient review are not
negative findings.

All existing T candidate groups were included with the previously frozen T
dictionary. Exact duplicate signatures are grouped before interpreting effects.
The primary family is all groups times admissible physical signatures; every
model is corrected across that complete family. Unadjusted, demographic and
overlap-removal sensitivities are separate families. These exploratory sender
associations do not supply independent clinical or direct-interaction validation.

Ambient sensitivity uses only removal already adopted by the atlas. It does
not assume all estimated contamination has been corrected. Raw/corrected
rank agreement cannot establish absence of other confounders.

Parent counts are compared for the exact retained libraries, not unmatched
source-wide totals. Historical uncertainty is resolved by the independent
count audit below; no old scientific result has been replaced.

## Audits
'''
    manifest = []
    for relative, value in audits:
        text += '\n### ' + relative + '\n\n```json\n' + json.dumps(value, indent=2) + '\n```\n'
        manifest.append(dict(path=relative, sha256=hashlib.sha256((root / relative).read_bytes()).hexdigest()))
    for title, relative in tables:
        frame = pd.read_csv(root / relative, sep='\t')
        text += '\n## ' + title + '\n\nFull authoritative table: `' + relative + '`\n\n'
        text += '```text\n' + frame.to_string(index=False) + '\n```\n'
        manifest.append(dict(path=relative, rows=len(frame), sha256=hashlib.sha256((root / relative).read_bytes()).hexdigest()))
    readiness = root / 'Bassez_fixed_function_readiness_20261014_v2/private_terminal.json'
    if readiness.exists():
        value = json.loads(readiness.read_text())
        assert value['status'].startswith('BASSEZ_FIXED_PUBLIC_FUNCTION_READINESS')
        text += '\n## Independent-cohort input readiness; no effect test\n\n```json\n' + json.dumps(value, indent=2) + '\n```\n'
        manifest.append(dict(path=str(readiness.relative_to(root)), sha256=hashlib.sha256(readiness.read_bytes()).hexdigest()))
    text += '''

## Next decision rules

Retain the unchanged original clinical family and all negative sender families.
A stable one-end clinical descriptor is a follow-up lead, not verified CM
clinical value or a direct T–myeloid axis. Source functional meanings must be
checked separately; do not assign a universal biological name by one source.

For an external cohort, first accept endpoint provenance, independent patient
keys, lineage capture and feature measurement. Bassez E/NE denotes T clonal
expansion, not pathological regression. A new fixed public partial functional
readout has a different estimand from old MCP transfer; the old95% weight
failure remains. Global columns do not prove per-source measurement, and
positive counts provide only a conservative measurement lower bound.

Do not rerun the current T screen with more thresholds to seek significance.
Plan any new independent-cohort readout family before looking at its effects,
retain all outputs and require molecular/spatial/functional validation before
claiming a cell interaction. The full research goal remains active.
'''
    out.mkdir()
    (out / 'FROZEN_TWO_END_PRIVATE_EVIDENCE.md').write_text(text)
    report = dict(status='AUDITED_TWO_END_FUNCTION_SENDER_AMBIENT_AND_DENOMINATOR_EVIDENCE_COMBINED_PRIVATE',
        finished_UTC=datetime.now(timezone.utc).isoformat(), inputs=manifest,
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        original_scientific_results_changed=False, formal_state_or_direct_interaction_verified=False,
        new_expression_download=False, goal_completed=False)
    (out / 'private_terminal.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'inputs'}))


if __name__ == '__main__':
    main()
