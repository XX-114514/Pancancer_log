"""Record the exact aggregate figures actually inspected in the authorized chat.

This is a visual review receipt, not a numerical or biological validation.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    resource(root, 1024**2)
    review = root / 'NSCLC_native_interface_association_review_20261018_v1'
    proof = review / 'private_independent_full_review_audit.json'
    assert sha(proof) == '738da26502c67aab026e8233ba4427be529a0a56191a4361638f7b446f3fd4a4'
    files = {
        'private_all120_record_count_overview.png': 'f855b396a7570370d97d350c738cf142836f6190d1634d984819c2a42f5145a3',
        'private_all1338_existing_sensitivity_comparisons.png': 'b901c85cd3fc0972c2bc2406f4ade3a5164b06c8ec2f441947ad7a1a7e614447',
        'private_all_internal_candidates_page_001.png': 'e2faa01b9e0d517544046cb68f40aca21fe5577810697a6d13e6c976467ae4db',
    }
    for name, digest in files.items():
        assert sha(review / name) == digest
    external = root / 'NSCLC_external_coordination_forests_20261018_v2'
    extra = 'private_NSCLC_MFA_20261016_F03_reviewed_patients_all_context_forests_v2.png'
    extra_sha = '1397485d9c602a9ee1123266c5e04122f8c2ee77426865b5d0d97c489b3c17b6'
    assert sha(external / extra) == extra_sha
    report = dict(
        status='ACTUAL_PIXEL_REVIEW_THREE_NATIVE_FIGURES_AND_ONE_ADDITIONAL_EXTERNAL_FOREST',
        recorded_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        native_derivation_audit_sha256=sha(proof), native_figures=files,
        additional_external_figure=dict(name=extra, sha256=extra_sha),
        visual_findings=['Titles, directional labels, axes, confidence intervals and annotations are readable.',
                         'Main family overview and all sensitivity panels retain their full original scope.',
                         'Existing external hollow points indicate unavailable intervals, not zero effects.',
                         'External material and phase labels remain visible; source contexts are not pooled.'],
        limits=['All three native images were actually viewed after full derivation acceptance.',
                'This adds one external forest to the previously recorded representative review; no claim that all20 external images were viewed.',
                'Numerical audit and visual inspection do not establish direct communication or clinical validity.'],
        new_models_P_or_endpoint_changes=False)
    target = review / 'private_actual_visual_review_20261018_v1.json'
    assert not target.exists()
    target.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(status=report['status'], recorded_UTC=report['recorded_UTC'], receipt_sha256=sha(target))))


if __name__ == '__main__':
    main()
