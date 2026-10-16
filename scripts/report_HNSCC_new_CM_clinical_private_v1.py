"""Combine the fully audited internal CM clinical family; private runtime values."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    result=root/'HNSCC_new_CM_clinical_results_20261016_v1'
    inputs=root/'HNSCC_new_CM_baseline_preparation_20261016_v1'
    model=root/'HNSCC_new_CM_frozen_model_export_20261016_v2'
    capture=root/'HNSCC_new_CM_capture_diagnostics_20261016_v1'
    out=root/'HNSCC_new_CM_clinical_private_evidence_20261016_v1';assert not out.exists()
    audit=json.loads((result/'private_independent_audit.json').read_text());assert audit['status']=='PASS_ALL30_CLINICAL_TESTS_FULL_RESAMPLING_BH_AND_ALL120_CONTROL_LOO'
    stats=result/'all30_private_clinical_statistics.tsv'
    assert hashlib.sha256(stats.read_bytes()).hexdigest()==audit['result_sha256']
    input_audit=json.loads((inputs/'private_independent_input_audit.json').read_text());assert input_audit['status'].startswith('PASS_ALL30')
    text='''# New atlas CM baseline clinical evidence (private)

## Identity, source and clinical scope

These are the new October2026 immune-atlas scope/component discoveries, not
the original October2026 v20261010 nine CM. All four eligible HNSCC scopes and
all30 component slots were fixed before clinical effects. Scope/component
overlap is retained, never counted as independent studies or30 distinct states.
All original history, labels and earlier single-end programme results remain.

Exact original author library keys and baseline timepoint were matched before
scoring. Each patient contributes one baseline composition per scope. The
clinical outcome is continuous pathological tumour regression, adjusted for
treatment and HPV. Author binary-label discrepancies were not repaired.

## Original model and baseline projection

The original NsNMF model has fitted matrix W S H. Omitting fixed smoothing S
would change coefficient meaning. Original R fitted values, saved W/H and
Python W S H were checked before clinical effects. Original minmax scaling is
fixed from the full discovery matrix. Independent BVLS reconstructs all
baseline NNLS coefficients and exact patient/clinical keys from raw frequency.
Original fitted H is replayed on its reconstructed matrix only; observed
baseline scores are new fixed-basis NNLS projections, not copied mixed-time H.

No W/S/scaling parameters were relearned. However, original unsupervised
discovery included post-treatment cells. The present analysis is internal
clinical association, not independent prediction or validated causality.

## Complete primary family and numerical audit

All30 primary slots share BH. Patient residual permutation uses20000 draws;
patient bootstrap2000 with ranks recomputed; fewer than1800 valid draws leaves
CI unavailable. Unestimable slots retain internal P1 for BH and are displayed
as missing. Main and all120 descriptive/LOO results were independently replayed
using manual ranks and full QR. Raw-W-only and age/sex remain descriptive.

'''
    manifest=[]
    sections=[('All30 primary clinical tests',stats),
        ('All120 control and LOO descriptors',result/'all120_private_descriptive_and_LOO.tsv'),
        ('All30 component coverage and structural interpretation',capture/'all30_private_capture_and_activity_diagnostics.tsv'),
        ('Cohort projection diagnostics',capture/'private_scope_projection_diagnostics.tsv')]
    for title,path in sections:
        table=pd.read_csv(path,sep='\t')
        text+='\n## '+title+'\n\n```text\n'+table.to_string(index=False)+'\n```\n'
        manifest.append(dict(path=str(path.relative_to(root)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),rows=len(table)))
    for title,path in [('Independent clinical audit',result/'private_independent_audit.json'),
        ('Independent baseline input audit',inputs/'private_independent_input_audit.json'),
        ('Original fixed model export',model/'private_terminal.json'),
        ('Capture diagnostic terminal',capture/'private_terminal.json')]:
        text+='\n## '+title+'\n\n```json\n'+path.read_text()+'\n```\n'
        manifest.append(dict(path=str(path.relative_to(root)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    text+='''

## Scientific decision boundary

Do not treat an uncorrected trend, rank label or isolated sensitivity as proof
of a clinically valuable CM. Limited active patients, repeated rank patterns
and high relative L2 projection residual identify applicability limits; they
are not an after-result justification to remove slots and refit FDR. Relative
L2 residual is not a centered variance-explained statistic or clinical R2.
Nor does lack of a robust association refute all CM biology or the original
nine modules in other cancers.

Preserve this complete clinical family and the previous programme evidence.
Do not keep changing normalization, rank or outcome cutoffs to seek a positive.
Continue the scheduled source-specific expression MCP / response / molecular
interface and qualified spatial directions. New discovery should keep a
patient/source separation and may require cohort-specific structure before
clinical generalisation. No direct interaction, resistance mechanism or
clinical prediction is verified by this report. The full goal remains active.
'''
    out.mkdir();path=out/'HNSCC_NEW_CM_PRIVATE_CLINICAL_EVIDENCE.md';path.write_text(text)
    report=dict(status='AUDITED_ALL30_NEW_CM_CLINICAL_AND_CAPTURE_EVIDENCE_COMBINED_PRIVATE',finished_UTC=datetime.now(timezone.utc).isoformat(),
        inputs=manifest,result_sha256=audit['result_sha256'],report_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        no_new_statistical_tests=True,original_20261010_9CM_changed=False,direct_interaction_verified=False,goal_completed=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='inputs'}))


if __name__=='__main__':main()
