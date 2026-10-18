"""Complete existing function/scope synthesis and source-only CAR role diagnostics."""
import json
import os
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    functional=root/'immune_RNA_MCP_functional_description_20261018_v1';fa=root/'immune_RNA_MCP_functional_description_audit_20261018_v1/private_terminal.json';assert sha(fa)=='e96be1d1a2630885557cec7443b931723a8f162f76dbff2b920f92b638f37b42'
    assert sha(functional/'private_terminal.json')==json.loads(fa.read_text())['producer_terminal_sha256']
    for n,h in json.loads((functional/'private_terminal.json').read_text())['files'].items():assert sha(functional/n)==h
    opp=root/'immune_RNA_MCP_all_atlas_opportunities_20261019_v1';oa=root/'immune_RNA_MCP_all_atlas_opportunities_audit_20261019_v1/private_terminal.json';assert sha(oa)=='7f181b10ce21de7b84b0ebb0f34ca88a34946afcaabba38a0c48db2b24157e0e';assert sha(opp/'private_terminal.json')==json.loads(oa.read_text())['producer_terminal_sha256']
    for n,v in json.loads((opp/'private_terminal.json').read_text())['files'].items():assert sha(opp/n)==v['sha256']
    out=root/'immune_RNA_MCP_function_and_atlas_scope_review_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_complete_function_scope_and_role_review.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,scope='All600 existing functional annotations/90source partitions and full400-material coverage. Fixed public author-label categories for CAR capture roles; clinical fields only availability/concordance, no programme scores or outcome associations.',
        interpretation='Hallmark gene-set enrichment is descriptive after learning. UV_RESPONSE_DN names a UV-down reference gene set, not evidence of UV exposure or exhaustion. Original paper4-1BB/Tpex/IL2 observations are comparison hypotheses, not current programme validation.',
        no_new_P_or_programme_fit_or_identity_change=True),indent=2)+'\n')
    annotation=pd.read_csv(functional/'private_all600_Hallmark_signed_annotations.tsv',sep='\t');components=pd.read_csv(functional/'private_all90_source_sentinel_RNA_component_descriptions.tsv',sep='\t');assert len(annotation)==600 and len(components)==90
    assert pd.api.types.is_integer_dtype(annotation.signed_active_genes.dtype)
    source_counts=pd.read_csv(opp/'private_all54_source_programme_pair_opportunities.tsv',sep='\t');contexts=pd.read_csv(opp/'private_all_context_programme_pair_opportunities.tsv',sep='\t',keep_default_na=False);lib=pd.read_csv(opp/'private_all400_material_RNA_capture_scope.tsv',sep='\t',keep_default_na=False)
    atlas=root.parents[1]/'runs/20261013_immunotherapy_atlas_v1';z=zarr.open_group(str(atlas/'outputs/immune_validation_authority_v1.zarr'),mode='r')
    src=np.asarray(read_elem(z['obs/dataset'])).astype(str);libraries=np.asarray(read_elem(z['obs/library_id'])).astype(str);labels=np.asarray(read_elem(z['obs/source_author_label'])).astype(str);target='CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea'
    pair=pd.read_csv(opp/'private_all2400_material_programme_pair_qualification.tsv.gz',sep='\t');valid=pair[pair.dataset.eq(target)&pair.pair.eq('R_FOXP3_capture|E_effector_capture')&pair.program.eq('MCP2')&pair.eligible]
    eligible=set(valid.sampleID);assert len(eligible)==58
    old=root/'immune_CM_fixed_panel_same_cell_RNA_20261018_v1';cache=np.load(old/'private_fixed34_gene_T_counts_and_depth.npz');raw=cache['raw_gene_counts'];genes=cache['genes'].tolist();rows=cache['global_original_T_row_index']
    fixed=root.parents[1]/'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json';panels=json.loads(fixed.read_text())['panels']
    def gate(names,n):return np.count_nonzero(raw[:,[genes.index(g) for g in names]]>=1,axis=1)>=n
    core=gate(panels['pan_T'],2);reg=gate(panels['regulatory_associated'],2)&gate(['FOXP3'],1);eff=gate(['GZMB','PRF1','CTSW','NKG7'],2)
    where=(src[rows]==target)&np.array([target+'::'+key in eligible for key in libraries[rows]])
    def category(label):
        label=str(label).lower()
        if 'regulatory' in label or 'treg' in label:return 'author_regulatory_T_annotation'
        if 'cd4' in label:return 'author_CD4_annotation'
        if 'cd8' in label:return 'author_CD8_annotation'
        if 'natural killer' in label or label.startswith('nk'):return 'author_NK_annotation'
        if 't cell' in label or 't-cell' in label:return 'author_T_annotation_unspecified'
        return 'other_or_unresolved_author_annotation'
    cats=np.array([category(x) for x in labels[rows]],str);roles=[]
    masks={'parent_core_T':where&core,'R_FOXP3_capture':where&core&reg&~eff,'E_effector_capture':where&core&eff&~reg}
    for view,mask in masks.items():
        counts=Counter(cats[mask]);n=int(mask.sum())
        for label,count in sorted(counts.items()):roles.append(dict(dataset=target,view=view,author_category=label,cells=count,view_cells=n,fraction=count/n,materials_with_pair_coverage=58,not_new_formal_identity=True))
    roles=pd.DataFrame(roles);roles.to_csv(out/'private_CAR_product_author_role_capture_summary.tsv',sep='\t',index=False)
    mask=src==target;clinical_fields=['library_id','patient_block','source_response','native_response_endpoint','clinical_endpoint_definition_status','clinical_validation_eligible','response_source']
    md=pd.DataFrame({f:np.asarray(read_elem(z['obs'][f])).astype(str)[mask] for f in clinical_fields});assert md.library_id.nunique()==59
    clinicalrows=[]
    for material,g in md.groupby('library_id',sort=True):
        assert g.patient_block.nunique()==1 and g.source_response.nunique()==1 and g.native_response_endpoint.nunique()==1
        response=g.source_response.iloc[0].strip().upper();clinicalrows.append(dict(library_id=material,patient_block=g.patient_block.iloc[0],source_response_class=response if response in ['R','NR'] else 'OTHER_OR_UNREVIEWED',native_response_endpoint=g.native_response_endpoint.iloc[0],definition_status=g.clinical_endpoint_definition_status.iloc[0],all_cells_legacy_clinical_eligible=bool(g.clinical_validation_eligible.str.lower().eq('true').all()),R_E_capture_eligible=target+'::'+material in eligible))
    clinical=pd.DataFrame(clinicalrows);clinical.to_csv(out/'private_all59_CAR_metadata_concordance_no_score_read.tsv',sep='\t',index=False)
    cs=dict(dataset=target,materials=len(clinical),reviewed_case_blocks=int(clinical.patient_block.nunique()),source_response_material_counts=clinical.source_response_class.value_counts().to_dict(),legacy_clinical_eligible_materials=int(clinical.all_cells_legacy_clinical_eligible.sum()),response_material_counts_in_R_E_coverage=clinical[clinical.R_E_capture_eligible].source_response_class.value_counts().to_dict(),definition_status_not_frozen_materials=int(clinical.definition_status.eq('not_frozen').sum()),new_clinical_endpoint_accepted=False,new_programme_scores_or_effects=False,
        context='CAR-T infusion product, not a tumour ICI sample. Endpoint timing/mapping and material/patient bridge must be independently reviewed before clinical models. Do not substitute24 unreviewed outcomes with negative response.')
    (out/'private_CAR_product_clinical_metadata_coverage.json').write_text(json.dumps(cs,indent=2)+'\n')
    refs=[dict(url='https://pubmed.ncbi.nlm.nih.gov/41045934/',role='HNSCC original article abstract/figure captions;4-1BB and Tpex hypotheses only; PMC full text captcha, not claimed read'),dict(url='https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_UV_RESPONSE_DN.html?ex=1',role='Reference definition: UV-down genes; no inferred UV exposure'),dict(url='https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_ALLOGRAFT_REJECTION.html',role='Reference describes genes up during transplant rejection, not a tumour-specific exhausted state'),dict(url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169246',role='CD45-enriched TNBC tumour/blood, paclitaxel alone or plus atezolizumab; separate contexts'),dict(url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176021',role='NSCLC source includes many libraries per patient and distinct tissue/functional branches; no110-independent-patient claim'),dict(url='https://pubmed.ncbi.nlm.nih.gov/37738975/',role='CAR infusion products, original resource endpoint described as3month PET/CT; current source labels still require exact bridge/definition review'),dict(url='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241783',role='Related59-patient CAR resource public primary record, alias to current CXG requires provenance qualification')]
    review=dict(reviewed_UTC=datetime.now(timezone.utc).isoformat(),functional_acceptance_sha256=sha(fa),all400_acceptance_sha256=sha(oa),all600_functional_annotations=annotation.to_dict(orient='records'),all90_source_background_partitions=components.to_dict(orient='records'),all54_source_pair_coverage=source_counts.to_dict(orient='records'),all_context_coverage=contexts.to_dict(orient='records'),author_role_summary=roles.to_dict(orient='records'),clinical_metadata_coverage=cs,primary_references=refs,
        interpretation='Current MCP2 shared R/E temporal clue has broad mixed functional annotation (E positive KRAS/allograft/UV-down, R positive UV-down and negative IFN-alpha). No basis for a specific exhausted/Treg-suppressive/resistance name. Fixed15 stress complement preserves prior direction but does not remove all processing effects.',
        next='Full frozen programme projection to all newly pair-eligible atlas materials, strict full axes and contexts/known patient blocks; all original programme/gene/quality parameters unchanged, original scores serve anchors. Include CAR products as separate domain/role diagnostic, not ICI TME replication. Exact clinical provenance first; freeze any new complete source-level test family before effects.',no_new_association_P_or_gene_fit=True)
    (out/'private_complete_function_scope_and_next_decision.json').write_text(json.dumps(review,indent=2)+'\n')
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(v['bytes'] for v in files.values())<8*1024**2
    t=dict(status='COMPLETE_ALL_FUNCTION_SCOPE_REVIEW_AND_CAR_ROLE_CLINICAL_METADATA_DIAGNOSTICS_NO_NEW_EFFECTS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),files=files)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
