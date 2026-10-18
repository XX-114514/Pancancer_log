"""Outcome-free complete projection scope and previously frozen auxiliary genes."""
import json
import os
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from measure_immune_RNA_MCP_fixed_lineage_diagnostic_v1 import PANELS


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    opp=root/'immune_RNA_MCP_all_atlas_opportunities_20261019_v1';audit=root/'immune_RNA_MCP_all_atlas_opportunities_audit_20261019_v1/private_terminal.json';assert sha(audit)=='7f181b10ce21de7b84b0ebb0f34ca88a34946afcaabba38a0c48db2b24157e0e';assert sha(opp/'private_terminal.json')==json.loads(audit.read_text())['producer_terminal_sha256']
    for n,v in json.loads((opp/'private_terminal.json').read_text())['files'].items():assert sha(opp/n)==v['sha256']
    old=root/'immune_RNA_MCP_disjoint_receiver_preflight_20261018_v1';oa=root/'immune_RNA_MCP_disjoint_receiver_preflight_audit_20261018_v1/private_terminal.json';assert sha(oa)=='001f82b90faf11e67e8a198653d90444df5e1881f754aee99e34dc4c3a62ab6e';assert sha(old/'private_terminal.json')==json.loads(oa.read_text())['producer_terminal_sha256']
    filename='private_all_unique_directional_endpoint_slots.tsv';ot=json.loads((old/'private_terminal.json').read_text());assert sha(old/filename)==ot['files'][filename]['sha256']
    registry=pd.read_csv(old/filename,sep='\t');assert len(registry)==394
    fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';capture_panels=json.loads(fixed.read_text())['panels']
    wp=root/'immune_RNA_MCP_frozen_projection_20261018_v1/private_R_projection/private_frozen_gene_weights.tsv';assert sha(wp)=='fd7316c67beb687c403a752c36b5989644e432f59a83119b21eade5efb0d89a7';weights=pd.read_csv(wp,sep='\t');weights=weights[weights.fit_id.eq('K2_seed20261018_source_adjusted')]
    roles=defaultdict(set)
    for g in set(weights.gene):roles[g].add('full_frozen_programme_gene')
    for group in capture_panels.values():
        for g in group:roles[g].add('original_fixed34_capture_gene')
    for group in PANELS.values():
        for g in group:roles[g].add('original_fixed43_diagnostic_gene')
    for r in registry.itertuples(index=False):
        for g in str(r.ligand).split('_'):roles[g].add('frozen394_ligand_component')
        for g in str(r.receptor).split('_'):roles[g].add('frozen394_receptor_component')
    lib=pd.read_csv(opp/'private_all400_material_RNA_capture_scope.tsv',sep='\t',keep_default_na=False);pairs=pd.read_csv(opp/'private_all2400_material_programme_pair_qualification.tsv.gz',sep='\t');axes=pd.read_csv(opp/'private_all54_source_axis_qualification.tsv',sep='\t').set_index(['dataset','view','program'])
    needed=set()
    for r in pairs[pairs.eligible].itertuples(index=False):
        for view in r.pair.split('|'):needed.add((r.sampleID,view,r.program))
    scope=[];views=['R_FOXP3_capture','E_effector_capture','Myeloid']
    for m in lib.itertuples(index=False):
        for view in views:
            for program in ['MCP1','MCP2']:
                complete=bool(axes.loc[(m.dataset,view,program),'complete_axis']);selected=(m.sampleID,view,program) in needed
                assert not selected or complete
                scope.append(dict(sampleID=m.sampleID,dataset=m.dataset,view=view,program=program,complete_gene_axis=complete,selected_for_projection=selected,status='PROJECT_FULL_FROZEN_AXIS' if selected else 'GENE_AXIS_UNAVAILABLE' if not complete else 'NO_QUALIFIED_PAIR_CAPTURE',reviewed_patient=m.reviewed_patient,case_block=m.case_block,phase_class=m.phase_class,material_class=m.material_class,capture_scope=m.capture_scope,previously_RNA_measured_material=m.previously_RNA_measured_material,original_discovery_material=m.original_discovery_material,reviewed_discovery_case_overlap=m.reviewed_discovery_case_overlap))
    table=pd.DataFrame(scope);assert len(table)==2400
    out=root/'immune_RNA_MCP_all_atlas_projection_contract_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,opportunity_acceptance_sha256=sha(audit),frozen394_acceptance_sha256=sha(oa),frozen_registry_sha256=sha(old/filename),full_weight_sha256=sha(wp),fixed_capture_panel_sha256=sha(fixed),fixed43_panel_source_sha256=sha(Path(__file__).with_name('measure_immune_RNA_MCP_fixed_lineage_diagnostic_v1.py')),
        scope='All400 materials/all2400 single-axis positions retained. Project every full axis participating in at least one qualified same-programme pair, including original materials as technical anchors; no clinical score/outcome selection.',
        primary='Fixed global-row hash2026101802, up100 cells/material/view. Full RNA CP10K, original gene/discovery quality scales and original nuisance coefficients; no programme/PCA/gene refit.',
        source_constants='Known HNSCC sourceCC0,ccRCC1 as original. Other source-specific additive constants unknown/omitted, no absolute cross-source score comparison.',
        sensitivities='All captured-cell score AND own quality means; fixed15 stress-complement scores on every full axis; optional raw gene component for measurement diagnostics, not replacement of main score.',
        auxiliary='All existing394 LR component genes plus original43 diagnostics and34capture genes measured in same raw pass; no new LR selection or tests. Per source unique measurement/complex missingness saved, never fill unmeasurable negative.',
        original_anchors='All396 official full programme material means plus all624 previously accepted HNSCC full/stress and all-capture means. Entire native capture/selection/slots/scores/means/axes independently replayed before new inference.',
        identity='RNA capture views remain not formal Treg/CD8/states. CAR product separate contextual/control domain, not tumour ICI replication; inherited ex_vivo false does not mean uncultured. Patients/contexts and all prior overlaps retained.',
        clinical='No outcome data read for this contract. CAR35legacy R/NR labels currently definition not_frozen; independently qualify exact endpoint/time/missingness/alias first. New inference source/pair/programme/clinical families frozen separately after assay acceptance, old families stay closed.',
        no_new_scores_or_effects_or_expression_download=True,prospective_assay_output_cap_bytes=256*1024**2,minimum_free_bytes=500000000000,additional_free_buffer_bytes=1073741824)
    (out/'plan_before_any_new_atlas_programme_expression.json').write_text(json.dumps(plan,indent=2)+'\n')
    table.to_csv(out/'private_all2400_frozen_projection_axis_scope.tsv',sep='\t',index=False)
    pd.DataFrame([dict(gene=g,roles='|'.join(sorted(tags))) for g,tags in sorted(roles.items())]).to_csv(out/'private_frozen_all_programme_LR_and_diagnostic_gene_roles.tsv',sep='\t',index=False)
    summary=dict(material_axis_positions=2400,selected_complete_axes=int(table.selected_for_projection.sum()),selected_materials=int(table.loc[table.selected_for_projection,'sampleID'].nunique()),sources=int(lib.dataset.nunique()),auxiliary_gene_role_records=len(roles),prior_fixed_LR_endpoints=394,clinical_outcome_or_new_programme_values_read=False)
    (out/'private_technical_projection_scope_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()}
    t=dict(status='FROZEN_COMPLETE_ALL_ATLAS_NATIVE_PROJECTION_AND_FIXED_AUXILIARY_MEASUREMENT_NO_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'),technical=summary)))


if __name__=='__main__':main()
