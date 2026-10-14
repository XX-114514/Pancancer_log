"""Private evidence extension; no gene lists or patient rows in stdout."""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    audit=root/'HNSCC_subtype_function_independent_audit_20261014_v1/private_audit.json'
    a=json.loads(audit.read_text());assert a['status'].startswith('PASS_ALL471')
    sub=root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1'
    funcs=root/'HNSCC_full50_orthogonal_function_diagnostics_20261014_v1'
    sex=root/'HNSCC_sex_function_descriptive_20261014_v1'
    external=root/'ICB_clinical_programme_scope_20261014_v2'
    out=root/'HNSCC_subtype_function_private_evidence_20261014_v2';assert not out.exists();out.mkdir()
    text='''# HNSCC subtype/function evidence extension (private)

This follows the previous frozen baseline clinical family, not a new set of
confirmatory clinical tests. All macrophage candidate subtypes, all50 Hallmarks,
all nuisance variants, and unestimable cases remain in scope. Original cell and
CM labels are unchanged; no formal cell-state assignment is made.

Within-subtype programme variation, subtype abundance, coverage and quality
sensitivity are separate questions. All subtype rows below determine their
direction and stability. A lineage candidate is not automatically a true state.

Function scores exclude both frozen T and Mac programme genes from signatures
AND rank background; fixed subtype cell samples and full measured native axes
are retained. All50 are measurable. Average-tie maxrank1500 UCell scalar checks
agree with independent count/order calculations. Functional associations and
all471 scalar/LOO diagnostics were independently replayed.

Functional associations and demographic sensitivity are reported in full below.
Gene-set collection names do not prove a cell identity or tissue phenotype.
Demographic category coverage must be checked before a claim of adequate
adjustment. These post-result diagnostics do not supply new significance.

## All subtype clinical descriptors

'''
    for directory,filename,title in [(sub,'all_private_subtype_descriptive_diagnostics.tsv','all subtype variants'),
      (funcs,'all300_private_descriptive_functions.tsv','all orthogonal functions'),
      (sex,'all12_private_sex_clinical_diagnostics.tsv','all sex/age clinical variants'),
      (sex,'all150_private_sex_function_diagnostics.tsv','all sex/age functional variants')]:
        text+='\n## '+title+'\n\n'+pd.read_csv(directory/filename,sep='\t').to_string(index=False)+'\n'
    text+='\n## Existing external source gates\n\n'+(external/'aggregate_terminal.json').read_text()
    text+='''

All8 supplemental ICB source opportunities were checked. Exact native/stable-ID
reconciliation and the fixed95% weight gate are reported in the native scope
record. Missing weights are neither declared zero nor compensated by changing
that threshold.
Bassez E/NE is T-cell clonal expansion, not HNSCC pathological regression.
Broad unsplit myeloid/T annotations are not an accepted old fine-node mapping;
donor-celltype minimum-cell filtering invalidates a population-census claim.
Overlap-quarantined studies stay quarantined. No external effect test ran.

## New public validation entrances and method assessment

Author site: https://lambrechtslab.sites.vib.be/en/single-cell
Author HNSCC site: https://lambrechtslab.sites.vib.be/en/hnscc
Bassez primary paper: https://www.nature.com/articles/s41591-021-01323-8
Franken primary DOI: 10.1016/j.immuni.2024.02.007
Franken spatial repository: https://doi.org/10.48804/992X8C
scIVCCC primary method: https://pmc.ncbi.nlm.nih.gov/articles/PMC13092271/

Current public author HTML links expose no processed-count download; no new
expression matrix was fetched. Franken spatial entry is controlled raw FASTQ,
with publicly readable README/inventory/access information. No access-control
bypass or application/message submitted. Published downloadable-count statements
are not counted as a currently accessible gene-complete matrix.

scIVCCC's treatment-instrument and conditional independence design is worth
considering after paired source/target/measurement review. For current HNSCC,
we have not established that treatment affects the receiver only through the
chosen sender, and pre biopsy/post resection differences remain. Non-rejection
of conditional-independence tests is not proof of mediation. Do not label a
new causal graph as validated based on method nomenclature alone.

## Next executable directions

Retain the original whole-programme result and this descriptive subtype lead.
Continue full source-specific programme interpretation, including the T end,
and seek independent biological/clinical evidence under explicit endpoint
definitions. A new transferable-gene-universe programme would need a NEW
phenotype-blind discovery/held-patient model and complete frozen transfer family;
it cannot replace the original95% failure by changing old weights. For CCC,
require sender programme, molecular interface, receiver response and appropriate
patient/source/spatial support. No new direct interaction is verified yet.
'''
    (out/'HNSCC_SUBTYPE_FUNCTION_PRIVATE_EVIDENCE.md').write_text(text)
    terminal=dict(status='FULL_SUBTYPE_FUNCTION_AND_EXTERNAL_GATE_EVIDENCE_RECORDED_PRIVATE',
      finished_UTC=datetime.now(timezone.utc).isoformat(),formal_state_or_interaction_verified=False,
      new_significance_tests=False,new_expression_download=False,goal_completed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
