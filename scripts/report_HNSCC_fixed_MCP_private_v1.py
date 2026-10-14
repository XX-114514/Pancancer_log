"""Build the full private HNSCC evidence card, including negative results."""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    clinical=root/'HNSCC_fixed_MCP_clinical_results_20261014_v2'
    paired=root/'HNSCC_fixed_MCP_paired_results_20261014_v1'
    quality=root/'HNSCC_native_quality_descriptive_20261014_v1'
    composition=root/'HNSCC_clinical_descriptive_composition_20261014_v1'
    out=root/'HNSCC_fixed_MCP_private_evidence_20261014_v1';assert not out.exists()
    c=pd.read_csv(clinical/'all24_private_clinical_statistics.tsv',sep='\t')
    p=pd.read_csv(paired/'all12_private_paired_statistics.tsv',sep='\t')
    d=pd.read_csv(composition/'all15_private_descriptive_diagnostics.tsv',sep='\t')
    q=pd.read_csv(quality/'all9_private_native_quality_diagnostics.tsv',sep='\t')
    ca=json.loads((clinical/'private_independent_audit.json').read_text())
    pa=json.loads((paired/'private_independent_audit.json').read_text())
    assert ca['status'].startswith('PASS_REAL_PRIMARY') and pa['status'].startswith('PASS_ALL6')
    main=c[(c.program=='MCP1')&(c.method=='primary_adjusted_continuous')]
    assert len(main)==3
    supported=main[main.BH_six_slot_family.lt(.05)&main.lower.mul(main.upper).gt(0)]
    paired_supported=p[p.BH_twelve_slot_family.lt(.05)]
    out.mkdir()
    text='''# HNSCC frozen CRC programme transfer: private evidence card

## Evidence identity

This is fixed-gene programme transfer into new-atlas multigene CD8 and
macrophage subtype candidates, not an accepted mapping to the old 172 nodes
and not a validation of an unchanged original CoVarNet CM. The original CM
numbers, cell labels, gene weights and quality coefficients remain unchanged.
Stable-ID aliases restored measurable features without fuzzy/expression-based
matching. Full-programme measurement still fails; >=95% fixed absolute-weight
partial-measurement sensitivity passes. Unmeasured centered contributions are
omitted, never RNA zeros or renormalized weights. CRC full/partial calibration
was recorded before clinical effects.

## Clinical definition

All author sample keys and treatment groups align exactly. Author binary
labels contain visit/threshold discrepancies and are preserved unchanged.
Continuous Tumor.regression is consistent across patient visits and is the
pre-score declared primary outcome. The >=50% secondary label is explicitly
derived from the paper definition, not a silent correction of the original.
The concordant-author-label sensitivity subset fails its minimum responder
coverage gate. No RECIST or survival conclusions.

## Primary clinical family

Baseline patient-level rank association adjusted for treatment arm and HPV.
Three endpoints times two original programme slots form six tests; the absent
MCP2 slot remains absent and participates internally as P=1 in correction.
20000 residual permutations and 2000 patient bootstraps, and all primary
coefficients/P/CI/BH, were independently replayed. No post-treatment scores
are described as prospective response predictors.

'''
    text+=main.to_string(index=False)+'\n\n## All clinical hypotheses\n\n'+c.to_string(index=False)
    text+='\n\n## Predeclared paired/coherence family\n\n'+p.to_string(index=False)
    text+='\n\n## All descriptive composition/LOO variants (no new P/BH)\n\n'+d.to_string(index=False)
    text+='\n\n## All descriptive native RNA/age/LOO variants (no new P/BH)\n\n'+q.to_string(index=False)
    text+='''

## Interpretation and remaining requirements

A macrophage-end programme clinical association is an observational lead,
not proof of a T--macrophage interaction, immunosuppression, drug resistance,
causal therapy effect or clinical prediction performance. A paired result
that does not pass the frozen family is retained; neither biopsy/resection
differences nor state composition are claimed to have been eliminated.
Post-result diagnostic adjustments give descriptive direction/LOO robustness,
not additional significant models. Smaller concordant-label subsets cannot
replace the original planned primary population.

Prior breast/glioma fixed-programme covariance and gastric counterexample
remain source-specific context; the HNSCC clinical signal does not turn them
into universal multicellular clinical effects. Review all frozen Hallmark
directions before a biological name, and retain old holdout and missingness.

Next: audit source/sample duplication and baseline clinical covariates;
separate sender/receiver and within-subtype state readouts under a frozen
complete candidate family; seek an eligible independent clinical or full-RNA
spatial cohort. Targeted panels that miss most programme weight cannot serve
as complete-programme validation. Do not tune endpoints, weight coverage,
programme signs, cell count thresholds or statistics to make a joint result.
The overall cell-interaction biological/clinical objective remains active.
'''
    (out/'HNSCC_FIXED_MCP_PRIVATE_EVIDENCE.md').write_text(text)
    audit=dict(status='HNSCC_PRIVATE_EVIDENCE_CARD_COMPLETE',finished_UTC=datetime.now(timezone.utc).isoformat(),
      supported_primary_endpoints=supported.endpoint.tolist(),supported_paired_questions=paired_supported.question.tolist(),
      direct_interaction_verified=False,joint_program_clinical_value_verified=False,
      original_author_binary_overwritten=False,original_CM_labels_changed=False)
    (out/'private_terminal.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit))


if __name__=='__main__':main()
