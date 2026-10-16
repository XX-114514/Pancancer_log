"""Join audited full-family statistics with all native molecular diagnostics."""
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    native = root / 'NSCLC_candidate_lineage_cofactor_RNA_20261016_v1'
    job = root / 'HNSCC_receiver_followup_job_control_20261016_v2/audit_NSCLC_candidate_lineage_cofactor_RNA_v1/launch_receipt.json'
    launch = json.loads(job.read_text())
    deadline = time.monotonic() + 14400
    while not (native / 'private_independent_audit.json').exists():
        proc = Path('/proc') / str(launch['pid']) / 'stat'
        assert proc.exists()
        state = proc.read_text().rsplit(')', 1)[1].split()
        assert state[19] == str(launch['proc_start_ticks']) and state[0] != 'Z'
        assert time.monotonic() < deadline
        time.sleep(15)
    audit = json.loads((native / 'private_independent_audit.json').read_text())
    assert audit['status'].startswith('PASS_ALL_NATIVE_UNIT_MARGINALS_')
    original = root / 'NSCLC_automatic_ligand_response_review_20261016_v1'
    original_audit = json.loads((original / 'private_independent_audit.json').read_text())
    assert original_audit['status'].startswith('PASS_ALL44528_VALUES_')
    evidence = original / 'private_all_primary_BH_supported_rows.tsv'
    terminal = json.loads((original / 'private_terminal.json').read_text())
    assert sha(evidence) == terminal['files'][evidence.name]
    table = pd.read_csv(evidence, sep='\t')
    summary_path = native / 'private_all_record_lineage_cofactor_summaries.tsv'
    summary = pd.read_csv(summary_path, sep='\t')
    native_terminal = json.loads((native / 'private_terminal.json').read_text())
    assert sha(summary_path) == native_terminal['files'][summary_path.name]
    assert set(summary.record) == set(table.record)
    rules_path = root / 'NSCLC_automatic_ligand_molecular_review_20261016_v1/private_public_literature_component_rules.json'
    rules = json.loads(rules_path.read_text())['ligand_rules']
    table['molecular_reference_status'] = table.ligand.map(lambda g: rules[g]['status'])
    joined = summary.merge(table, on='record', how='left', validate='many_to_one')
    assert len(joined) == len(summary) and joined.molecular_reference_status.notna().all()
    out = root / 'NSCLC_molecular_and_native_evidence_review_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    joined.to_csv(out / 'private_all_native_diagnostics_with_original_evidence.tsv.gz', sep='\t', index=False, na_rep='NA', compression='gzip')
    report = ['# 全部原主筛记录的分子及原生诊断', '', '原完整家族统计和各类临床/协调注释不变；下面是额外RNA诊断，不创造新的P、FDR或独立验证。', '', '所有>=1/2阈值及阴性/低覆盖项均保留。谱系标志缺失可能源于低测序深度，异类标志可有多种解释，不据此自动删细胞。', '', '每一受体/配体复合物需要同细胞组分；RNA不能代替表面蛋白、释放/剪切、氧化还原或功能传导。', '']
    columns = ['role', 'author_subtype', 'group', 'threshold', 'measured', 'cases', 'native_cells', 'positive_cases', 'GOI_positive_cells', 'GOI_positive_fraction', 'core_lineage_among_GOI_positive', 'mesenchymal_markers_co_positive_cells', 'endothelial_markers_co_positive_cells', 'epithelial_markers_co_positive_cells', 'platelet_markers_co_positive_cells']
    for row in table.itertuples():
        report.extend([f'## {row.record}', '', f'{row.ligand}: {row.sender_author_subtype} -> {row.receiver_author_subtype}; original possible receptor(s): {row.qualified_receptors}.', '', f'Literal reference interpretation: {row.molecular_reference_status}', '', f'Original rho={row.partial_rank_rho:.6g}, BH={row.BH_q_complete_unique_family:.6g}, BY={row.BY_q_dependence_sensitivity:.6g}; sender-programme additional-control rho={row.sender_programme_adjusted_rho:.6g}. Original clinical annotations concern the FULL signature, not this masked outcome.', '', summary[summary.record.eq(row.record)][[c for c in columns if c in summary]].to_string(index=False), ''])
    report.extend(['## 下一步', '', '基因不重叠的接收TF读数仅是同来源探索性功能检查；外部样本、蛋白实验/阻断、明确患者随访及跨来源仍是不同证据。不能把必要组分补齐后的替代轴继承原q值，不能直接命名免疫抑制或耐药。'])
    (out / 'PRIVATE_MOLECULAR_NATIVE_REVIEW.md').write_text('\n'.join(report) + '\n')
    result = dict(status='COMPLETE_PRIVATE_JOIN_OF_ALL_AUDITED_NATIVE_DIAGNOSTICS_AND_ORIGINAL_STATISTICAL_EVIDENCE_NO_NEW_TESTS', finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), original_evidence_audit_sha256=sha(original / 'private_independent_audit.json'), native_audit_sha256=sha(native / 'private_independent_audit.json'), molecular_rules_sha256=sha(rules_path), parents=len(table), joined_rows=len(joined), no_new_P_BH_or_identity_changes=True, files={p.name: sha(p) for p in out.iterdir() if p.is_file()}, receptor_specific_or_causal_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
