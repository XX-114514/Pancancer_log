"""Explicit unit-level analysis gates; unknown is distinct from rejected.

Inputs must come from source/asset and expression-unit review, never inferred
from tumor labels, the presence of immune cells, model success or sample names.
This computes eligibility only. It does not assign patient identities or run CM.
"""
ACCEPTED='accepted'
REJECTED='rejected'
UNKNOWN='unknown'
STATUSES={ACCEPTED,REJECTED,UNKNOWN,'not_applicable'}
CAPTURES={'unselected_tissue','whole_tissue_sampling_candidate','immune_enriched','T_enriched_nuclei','selected_lineage','selected_nonmalignant_subset','cell_product','mixed_sorted_fractions','unknown'}


def _requirements(review, fields, restrictions=()):
    missing=[field for field in fields if review.get(field,UNKNOWN)!=ACCEPTED]
    rejected=[field for field in missing if review.get(field)==REJECTED]
    reasons=list(restrictions)+[f'{field}:{review.get(field,UNKNOWN)}' for field in missing]
    status=REJECTED if restrictions or rejected else UNKNOWN if missing else ACCEPTED
    return {'status':status,'eligible':status==ACCEPTED,'reasons':reasons,'required_review_fields':list(fields)}


def qualify_unit_views(review):
    """Return explicit archive/composition/reference/product/clinical views.

    Shared source reviews and per-unit reviews must first be joined explicitly.
    Unknown provenance stays unknown even when observed values look plausible.
    The caller retains source, native unit, library, sample and patient keys in
    separate private columns; no native key is promoted to a patient here.
    """
    if review.get('archive_membership')!=ACCEPTED:
        raise ValueError('This view operates on already accepted archive members')
    capture=review.get('capture_class','unknown')
    if capture not in CAPTURES:raise ValueError('Capture class requires explicit review')
    for field,value in review.items():
        if field!='capture_class' and value not in STATUSES:
            raise ValueError('Review fields must use explicit accepted/rejected/unknown/not_applicable status')
    # Native evidence can be referenced without a finalized identity. Labeled
    # composition and clinical analyses require the separate annotation gate.
    expression=('source_asset_join','native_measurement','prior_QC_scope','expression_unit_join','assay_scope')
    composition=(*expression,'annotation_release','source_alias_scope','material_and_anatomy','capture_design')
    non_tissue=() if capture in ['unselected_tissue','whole_tissue_sampling_candidate','unknown'] else ('capture_not_unselected_whole_tissue:'+capture,)
    tissue_fields=(*composition,'unselected_tissue_sampling')
    capture_unknown=('capture_class:'+capture,) if capture in ['unknown','whole_tissue_sampling_candidate'] else ()
    tissue=_requirements(review,tissue_fields,non_tissue)
    if capture_unknown and tissue['status']==ACCEPTED:
        tissue.update(status=UNKNOWN,eligible=False,reasons=list(capture_unknown))
    within=_requirements(review,(*composition,'lineage_selection_scope','within_lineage_capture_comparability'),('cell_product_requires_product_specific_view',) if capture=='cell_product' else ())
    if capture=='unknown' and within['status']==ACCEPTED:
        within.update(status=UNKNOWN,eligible=False,reasons=['capture_class:unknown'])
    product=_requirements(review,(*expression,'annotation_release','product_design'),() if capture=='cell_product' else ('not_a_qualified_cell_product_capture',))
    clinical=_requirements(review,(*composition,'patient_identity','patient_sample_library_join','treatment_phase','clinical_endpoint','endpoint_time_origin','repeated_sampling_plan','clinical_cohort_scope'))
    return {'all_archive':{'status':ACCEPTED,'eligible':True,'reasons':[]},
            'expression_reference':_requirements(review,expression),
            'whole_tissue_composition':tissue,
            'within_selected_lineage_composition':within,
            'cell_product':product,
            'treatment_clinical':clinical,
            'source_alias_acceptance_does_not_prove_independent_patients':True,
            'new_patient_or_identity_assignment':False,
            'CM_or_communication_analysis_performed':False}
