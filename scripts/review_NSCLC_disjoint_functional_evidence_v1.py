"""Private complete functional evidence, shared-response diagnostics and figures."""
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    d = root / 'NSCLC_gene_disjoint_TF_response_descriptions_20261016_v1'
    s = root / 'NSCLC_TF_specificity_shared_response_20261016_v1'
    e = root / 'existing_immune_cohort_candidate_extension_preflight_20261016_v2'
    receipt = json.loads((root / 'HNSCC_receiver_followup_job_control_20261016_v2/audit_existing_immune_candidate_extension_inputs_v3/launch_receipt.json').read_text())
    deadline = time.monotonic() + 3600
    while not (e / 'private_independent_audit.json').exists():
        p = Path('/proc') / str(receipt['pid']) / 'stat'
        assert p.exists()
        state = p.read_text().rsplit(')', 1)[1].split()
        assert state[19] == str(receipt['proc_start_ticks']) and state[0] != 'Z' and time.monotonic() < deadline
        time.sleep(10)
    for folder in [d, s, e]:
        audit = json.loads((folder / 'private_independent_audit.json').read_text())
        assert audit['status'].startswith('PASS_')
        assert sha(folder / 'private_terminal.json') == audit['producer_terminal_sha256']
        terminal = json.loads((folder / 'private_terminal.json').read_text())
        for name, digest in terminal['files'].items():
            assert sha(folder / name) == digest
    points = pd.read_csv(d / 'private_all2568_TF_functional_descriptions.tsv', sep='\t')
    specific = pd.read_csv(s / 'private_all_regulator_specificity_and_shared_response.tsv', sep='\t')
    # Overlapping fields are checked, then retained once with original names.
    assert specific.TF_record.is_unique
    for row in specific.itertuples():
        old = points[points.TF_record.eq(row.TF_record)].iloc[0]
        assert old.record == row.record and old.TF == row.TF and old.ligand == row.ligand
        assert np.isclose(old.partial_rank_rho, row.original_partial_rank_rho)
    extra = specific.drop(columns=['record', 'TF', 'ligand', 'original_partial_rank_rho', 'original_sender_programme_adjusted_rho'])
    merged = points.merge(extra, on='TF_record', how='left', validate='one_to_one')
    assert len(merged) == len(points)
    merged['regulator_interpretation'] = np.where(merged.partial_rank_rho.isna(), 'ORIGINAL_UNAVAILABLE_SLOT', np.where(merged.literal_regulator_gene_measured.eq(False), 'NONLITERAL_COMPLEX_NOT_ABSENT', np.where(merged.receiver_regulator_RNA_positive_cells.eq(0), 'NO_NATIVE_REGULATOR_RNA_DETECTED_SIGNATURE_ONLY', 'TARGET_SIGNATURE_SOME_REGULATOR_RNA_NOT_PROTEIN_ACTIVITY')))
    out = root / 'NSCLC_disjoint_functional_evidence_review_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    merged.to_csv(out / 'private_all2568_functional_evidence_and_limits.tsv.gz', sep='\t', index=False, na_rep='NA', compression='gzip')
    external = pd.read_csv(e / 'private_all_source_parent_coverage.tsv', sep='\t')
    external.to_csv(out / 'private_all108_external_input_opportunities.tsv', sep='\t', index=False)
    pc = pd.read_csv(s / 'private_all_parent_residual_PC_descriptions.tsv', sep='\t')
    report = ['# 完整功能RNA及外部输入机会审阅', '', f'UTC {datetime.now(timezone.utc).isoformat()}', '', '全部2568原槽、1570可测调控签名保留；原LR家族P/BH不变。本轮只是同来源筛选后的方向和敏感性描述，无新P/BH或独立验证。', '', '## 名称与证据边界', '', '这里测量有符号调控靶RNA模式，不是TF蛋白、DNA结合或磷酸化活性。原生调控因子mRNA为零也不证明蛋白不存在；非字面复合物不是缺表达。接收签名可能有共享靶基因/共同表达程序。', '', '加入发送者同一调控签名是共同程序诊断，也可能控制中介，不能据衰减直接证明混杂或排除通讯。路径符号是受体到调控因子参考方向，不能自动说明配体激活受体。', '', '## 全部原端点的描述', '']
    for parent, f in merged[merged.partial_rank_rho.notna()].groupby('record', sort=True):
        first = f.iloc[0]
        columns = ['TF', 'partial_rank_rho', 'sender_programme_adjusted_rho', 'sender_same_regulator_adjusted_rho', 'receiver_regulator_RNA_positive_fraction', 'maximum_absolute_residual_score_correlation', 'regulator_interpretation']
        strongest = f.assign(absrho=f.partial_rank_rho.abs()).sort_values(['absrho', 'TF'], ascending=[False, True]).head(8)
        report += [f'### {parent}: {first.ligand}, {first.sender_author_subtype} -> {first.receiver_author_subtype}', '', '下面仅显示全部描述中绝对系数最高8项用于审阅，不另定义显著家族或筛选后确认。完整表保留所有项。', '', strongest[columns].to_string(index=False), '', pc[pc.record.eq(parent)].to_string(index=False), '']
    report += ['## 外部来源', '', '以下只有测量和细胞池机会，不是外部效应。原作者细亚群未一一映射。CD8_Effector/Treg等是已有图谱候选、未正式晋级；不会替换原身份。CAR产品/离体刺激、不同组织/时点单列，不能将多个文库当独立患者。原临床8项控制并非在所有外部来源可恢复，不能静默删控制拟合。', '', external.drop_duplicates(['dataset', 'ligand', 'sender_diagnostic_pool', 'receiver_diagnostic_pool']).to_string(index=False), '', 'v2按已有稳定ID的pipe别名精确匹配，排除歧义和物理列碰撞，保留95%门槛；v1只有初步字面匹配。后续仍需正式病例/取样和临床sidecar及精确端点/新谱系扩展范围预冻结。', '', '## 下一步', '', '对有生物结合依据、原生来源支持及可解释功能读数的候选，规划独立来源中的固定验证。保留全部原分子家族及否定/不确定项；不新增同来源显著性筛选来宣称机制。']
    (out / 'PRIVATE_COMPLETE_FUNCTIONAL_REVIEW.md').write_text('\n'.join(report) + '\n')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none', 'font.size': 9})
    colors = {'SDC2': '#2e6f9e', 'ADAM17': '#ba7c29', 'HMGB1': '#8b5ca6', 'LTB': '#548856'}
    fig, ax = plt.subplots(figsize=(7.5, 6))
    for ligand, f in merged[merged.partial_rank_rho.notna()].groupby('ligand'):
        ax.scatter(f.partial_rank_rho, f.sender_same_regulator_adjusted_rho, s=11, alpha=.45, label=ligand, color=colors[ligand])
    ax.plot([-1, 1], [-1, 1], color='#777777', linestyle='--', lw=.8)
    ax.axhline(0, color='#aaaaaa', lw=.6)
    ax.axvline(0, color='#aaaaaa', lw=.6)
    ax.set(xlim=(-1, 1), ylim=(-1, 1), xlabel='Original full-control partial rank correlation', ylabel='After adding sender same-regulator signature', title='All 1,570 fixed functional descriptions; one NSCLC source')
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out / 'private_all_functional_shared_response.png', dpi=160)
    fig.savefig(out / 'private_all_functional_shared_response.svg')
    plt.close(fig)
    # A named display follow-up, not a selected inferential family.
    selected = merged[merged.ligand.eq('SDC2') & merged.TF.isin(['FOXO3', 'NR3C1'])].sort_values(['TF', 'record'])
    selected.to_csv(out / 'private_display_FOXO3_NR3C1_all_six_endpoints.tsv', sep='\t', index=False)
    fig, ax = plt.subplots(figsize=(10, 6))
    y = np.arange(len(selected))
    for k, row in enumerate(selected.itertuples()):
        ax.plot([row.partial_rank_rho, row.sender_same_regulator_adjusted_rho], [k, k], color='#888888', lw=1)
    ax.scatter(selected.partial_rank_rho, y, label='Original controls', color='#2e6f9e', s=30)
    ax.scatter(selected.sender_same_regulator_adjusted_rho, y, label='+ sender same-regulator signature', color='#ba7c29', marker='^', s=32)
    ax.set_yticks(y, [f'{r.TF} | {r.receiver_author_subtype}' for r in selected.itertuples()], fontsize=8)
    ax.invert_yaxis()
    ax.set(xlim=(-.02, .5), xlabel='Descriptive partial rank correlation; no new P or BH', title='SDC2 follow-up display selected after full exploration; not independent validation')
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    fig.tight_layout()
    fig.savefig(out / 'private_FOXO3_NR3C1_display.png', dpi=160)
    fig.savefig(out / 'private_FOXO3_NR3C1_display.svg')
    plt.close(fig)
    result = dict(status='COMPLETE_PRIVATE_ALL_FIXED_FUNCTIONAL_EVIDENCE_EXTERNAL_INPUTS_AND_TWO_DESCRIPTIVE_FIGURES', finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), functional_audit_sha256=sha(d / 'private_independent_audit.json'), specificity_audit_sha256=sha(s / 'private_independent_audit.json'), external_input_audit_sha256=sha(e / 'private_independent_audit.json'), all_fixed_rows=len(merged), external_rows=len(external), figures=2, files={p.name: sha(p) for p in out.iterdir() if p.is_file()}, visual_QA_pending=True, no_new_P_BH_or_independent_validation=True)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
