"""Complete private scientific review, preserving author/candidate distinction."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 256 * 1024)
    audit = root / 'gamma_group_lineage_orthogonal_description_audit_20261019_v1/private_terminal.json'
    assert sha(audit) == '32d0dba881c9029c61839ccef8b449a27ddb1142453c3d889551fcff3a1e3ec9'
    src = root / 'gamma_group_lineage_orthogonal_description_20261019_v1'
    assert json.loads(audit.read_text())['producer_terminal_sha256'] == sha(src / 'private_terminal.json')
    out = root / 'gamma_group_lineage_scientific_review_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    hroot = root / 'public_GSE243013_TCR_header_qualification_20261019_v1'
    ht = wait(root, out, 'qualify_public_GSE243013_TCR_header_v1', hroot / 'private_terminal.json', {'QUALIFIED_PUBLIC_PROCESSED_TCR_HEADER_ONLY_EXPLICIT_ALPHA_BETA_NOT_GAMMA_CONTENT_VERIFIED'})
    assert sha(hroot / 'private_header_qualification.json') == ht['header_qualification_sha256']
    header = json.loads((hroot / 'private_header_qualification.json').read_text())
    assert not header['explicit_TRG_TRD_columns_present'] and not header['complete_file_downloaded']
    plan = json.loads((src / 'plan_before_lineage_statistic_or_TCR_schema.json').read_text())
    table = pd.read_csv(src / 'private_all100_within_author_record_marker_contrasts.tsv', sep='\t')
    descriptors = []
    for (a, b), g in table.groupby(['reference_gamma_group', 'comparison_group']):
        assert set(g.marker) == set(sum(plan['panels'].values(), []))
        for panel, genes in plan['panels'].items():
            v = g[g.marker.isin(genes)].median_detection_difference
            descriptors.append(dict(source='GSE243013', reference=a, comparison=b, panel=panel, fixed_markers=len(v), positive_paired_median_differences=int(v.gt(0).sum()), negative_paired_median_differences=int(v.lt(0).sum()), equal_paired_median_differences=int(v.eq(0).sum()), unknown=int(v.isna().sum()), comparable_author_records=int(g.comparable_author_records.iloc[0]), reviewed_independent_patients=0, inferential_test=False))
    assert len(descriptors) == 28
    pd.DataFrame(descriptors).to_csv(out / 'private_all_four_comparisons_seven_panel_direction_descriptions.tsv', sep='\t', index=False)
    request = root / 'private_gamma_group_scientific_review_request_20261019_v1.json'
    assert sha(request) == 'c560c553827a02d35acb14d0f83c3a7c99662a4b688dc66753aa1b96ff93168b'
    expert = json.loads(request.read_text())
    decision = dict(reviewed_UTC=datetime.now(timezone.utc).isoformat(), source='GSE243013', full_group_profile_audit_sha256=sha(audit), public_TCR_header_terminal_sha256=sha(hroot / 'private_terminal.json'), scientific_status=expert['scientific_status'], findings=expert['findings'], next_priorities=expert['next_priorities'], preserve_all_author_references=True, preserve_all_existing_raw_and_pilot_candidate_flags=True, no_final_label_or_state_write=True, no_new_communication_P=True, scientific_outputs_public=False)
    (out / 'private_scientific_review_and_next_scope.json').write_text(json.dumps(decision, indent=2) + '\n')
    receipt = dict(status='COMPLETE_FULL_GROUP_LINEAGE_SCIENTIFIC_REVIEW_AND_PUBLIC_HEADER_SCOPE_NO_FORMAL_GAMMA_ADOPTION', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), private_scientific_decision_sha256=sha(out / 'private_scientific_review_and_next_scope.json'), profile_audit_sha256=sha(audit), public_header_terminal_sha256=sha(hroot / 'private_terminal.json'), public_scientific_results=False, resource_guard=guard)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
