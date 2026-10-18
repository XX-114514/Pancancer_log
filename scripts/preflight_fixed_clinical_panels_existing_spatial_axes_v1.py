"""All unchanged clinical functional panels on all existing targeted spatial axes.

Availability only; no expression, phenotype, geometry or clinical tests.
"""
import hashlib
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget=resource(root,2*1024**2)
    state=Path(json.loads((root/'state_discovery_evidence_register_20261014_v1.json').read_text())['model_output_root'])
    folder=state/'state_spatial_proxy_calibration_inputs_20261016_v1';inventory=json.loads((folder/'inventory.json').read_text())
    path=folder/'fixed_all14_panel_gene_axes.json';assert sha(path)==inventory['files'][path.name]
    assays=json.loads(path.read_text());assert len(assays)==14
    source=Path(os.environ['PANCANCER_PROJECT_ROOT'])/'Project_v5/runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    fixed=json.loads(source.read_text());out=root/'fixed_clinical_panels_existing_spatial_axes_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),budget=budget,all14_accepted_axes_sha256=sha(path),original_clinical_panel_plan_sha256=sha(source),scope='All7 unchanged panels and both unchanged functional axes on all14 existing assays. Complete gene set required; shared panels not independent technical replications.',limits='Metadata gene availability only; no new expression downloads, scores, state labels, spatial effects or patient independence. Missing is unmeasurable, never negative.')
    (out/'plan_before_complete_fixed_panel_coverage.json').write_text(json.dumps(plan,indent=2)+'\n')
    rows=[];axisrows=[]
    for assay in assays:
        genes=set(assay['unambiguous_symbols']);digest=hashlib.sha256(('\n'.join(sorted(genes))+'\n').encode()).hexdigest();assert digest==assay['gene_axis_sha256']
        common=dict(dataset_id=assay['dataset_id'],sample_id=assay['sample_id'],gene_axis_sha256=digest,input_sha256=assay['input_sha256'])
        coverage={}
        for panel,required in fixed['panels'].items():
            complete=set(required)<=genes;coverage[panel]=complete
            rows.append(dict(**common,panel=panel,required_genes=len(required),unique_measured_genes=len(set(required)&genes),complete_fixed_panel=complete,missing_genes='|'.join(sorted(set(required)-genes))))
        for axis,(a,b) in fixed['axes'].items():axisrows.append(dict(**common,axis=axis,both_arm_sets_complete=coverage[a] and coverage[b],pan_T_complete=coverage['pan_T'],interferon_complete=coverage['interferon'],proliferation_complete=coverage['proliferation'],all_seven_panels_complete=all(coverage.values()),same_clinical_assay_or_model_validated=False))
    pd.DataFrame(rows).to_csv(out/'private_all98_fixed_panel_assay_coverage.tsv',sep='\t',index=False)
    pd.DataFrame(axisrows).to_csv(out/'private_all28_functional_axis_assay_coverage.tsv',sep='\t',index=False)
    files=[p for p in out.iterdir() if p.is_file()]
    terminal=dict(status='COMPLETE_ALL98_PANEL_AND28_AXIS_MEASUREMENT_AVAILABILITY_ONLY',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),assays=len(assays),unique_gene_axes=len({r['gene_axis_sha256'] for r in assays}),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in files},new_expression_or_spatial_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(dict(status=terminal['status'],completed_UTC=terminal['completed_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))

if __name__=='__main__':main()
