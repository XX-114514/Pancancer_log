"""Review complete source aggregates and plot all frozen opportunities privately."""
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import argparse
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import GSE288199_state_pair_math_v1 as math
import NSCLC_TCR_reported_expansion_statistics_v1 as family


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def short(name):
    return name.replace('CD4T_Treg_', 'Treg ').replace('CD8T_', '')


def cell(result):
    if result.get('rho') is None:
        return 'NA: '+result['status']
    ci = result.get('CI95')
    prefix = 'coordination40' if result['kind'] == 'coordination' else 'conditional80'
    return f"{result['rho']:+.3f}; CI "+(f'[{ci[0]:+.3f}, {ci[1]:+.3f}]' if ci else 'NA')+f"; BH {result[prefix+'_BH']:.4g} / BY {result[prefix+'_BY']:.4g}"


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--version', type=int, choices=[1, 2], default=1)
    a = p.parse_args()
    t = a.task
    read = lambda name: json.loads((t/name).read_text())
    plan = read('private_frozen_local_statistics_v1.json')
    terminal = read('local_statistics_v1/private_terminal.json')
    transport = read('private_local_RNA_transport_terminal_v1.json')
    measurement = read('private_complete_local_RNA_bridge_v1.json')
    results = read('local_statistics_v1/private_all120_components.json')
    primary = read('local_statistics_v1/private_primary40.json')
    points = read('local_statistics_v1/private_all1080_point_diagnostics.json')
    audits = read('local_statistics_v1/private_all_real_random_stream_R_audits.json')
    calibrations = read('local_statistics_v1/private_all_fixed_working_null_diagnostics.json')
    assert transport['local_payload_released'] and terminal['all120_slots_and40_pairs_retained']
    assert sha(t/'private_local_RNA_bridge_terminal_v1.json') == plan['measurement_terminal_SHA256']
    assert sha(t/'private_complete_local_RNA_bridge_v1.json') == plan['measurement_summary_SHA256']
    assert terminal['started_UTC'] > plan['frozen_UTC']
    for name, digest in plan['code_SHA256'].items():
        assert sha(Path(__file__).parent/name) == digest
    pairs = list(map(tuple, plan['all40_pairs']))
    assert len(results) == 120 and len(points) == 1080 and len(primary) == 40
    lookup = {(v['R_group'], v['E_group'], v['kind']): v for v in results}
    assert len(lookup) == 120
    assert set(lookup) == {(r, e, k) for r, e in pairs for k in family.KINDS}
    fitted = [v for v in results if v['status'] == 'ESTIMABLE']
    assert len(calibrations) == len(fitted)
    by_label = {v['label']: v for v in audits}
    assert len(by_label) == len(audits)
    for v in fitted:
        key = '|'.join([family.NAMESPACE, v['kind'], v['R_group'], v['E_group']])
        assert by_label[key]['full_wild_statistics_checked'] == 19999 and by_label[key]['full_bootstrap_statistics_checked'] == 2000
    fitted_labels = {'|'.join([family.NAMESPACE, v['kind'], v['R_group'], v['E_group']]) for v in fitted}
    assert all(v['full_wild_statistics_checked'] == 0 and v['full_bootstrap_statistics_checked'] == 0 for v in audits if v['label'] not in fitted_labels)
    assert max(v['max_numeric_error'] for v in audits) < 1e-8
    for v in fitted:
        assert np.isfinite(v['rho']) and abs(v['rho']) <= 1
        assert 1/20000 <= v['raw_P'] <= 1
        assert (v['CI95'] is not None) == (v['bootstrap_valid'] >= 1900)
        assert sum(v['bootstrap_status_counts'].values()) == 2000
        assert v['residual_df'] >= 10
        if v['CI95'] is not None:
            assert -1 <= v['CI95'][0] <= v['CI95'][1] <= 1
    for name, table, key in family.families(results, pairs):
        original = read('local_statistics_v1/private_'+name+'.json')
        assert len(original) == len(table)
        for u, v in zip(original, table):
            assert u == v
    full = [v for v in points if v['setting'] == 'full']
    for v in full:
        original = lookup[v['R_group'], v['E_group'], v['kind']]
        assert v['n'] == original['n'] and v['status'] == original['status'] and v['rho'] == original.get('rho')
    null_ok = all(s['passed_finite_alarm'] for v in calibrations for s in v['scenarios'])
    assert ('WORKING_NULL_ACCEPTED' in terminal['status']) == null_ok
    out = t/('source_figures_v'+str(a.version))
    out.mkdir(exist_ok=False)
    rs = list(dict.fromkeys(r for r, e in pairs))
    es = list(dict.fromkeys(e for r, e in pairs))
    coverage = {(v['R_group'], v['E_group']): v for v in measurement['all40_clinical_capture_positions']}
    fig, axes = plt.subplots(1, 3, figsize=(19, 4.8), constrained_layout=True)
    for ax, key, title in zip(axes, ['qualified_author_materials', 'response1', 'response0'],
                             ['Qualified author material units', 'MPR / pCR units', 'Non-MPR units']):
        values = np.asarray([[coverage[r, e][key] for e in es] for r in rs])
        im = ax.imshow(values, vmin=0, vmax=max(1, values.max()), cmap='Blues', aspect='auto')
        for i, r in enumerate(rs):
            for j, e in enumerate(es):
                ax.text(j, i, str(values[i, j]), ha='center', va='center', fontsize=8,
                        color='white' if values[i, j] > values.max()*.65 else 'black')
        ax.set_xticks(range(len(es)), [short(e) for e in es], rotation=60, ha='right', fontsize=8)
        ax.set_yticks(range(len(rs)), [short(r) for r in rs], fontsize=9)
        ax.set_title(title, fontsize=11)
        fig.colorbar(im, ax=ax, shrink=.7)
    fig.suptitle('GSE243013: all 40 fixed TCR / RNA type opportunities\nAuthor material units; patient independence unreviewed', fontsize=13)
    fig.savefig(out/'all40_coverage.png', dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(20, 16), sharey=True)
    colors = ['#3b6f94', '#3d887a', '#8c6391', '#858585']
    for ax, kind, title in zip(axes, family.KINDS, ['R–E expansion coordination', 'R–pathology conditional on E', 'E–pathology conditional on R']):
        for j, (r, e) in enumerate(pairs):
            v = lookup[r, e, kind]
            color = colors[rs.index(r)]
            if v.get('rho') is None:
                ax.text(.04, j, 'NA', fontsize=8, color='#777777', va='center')
                continue
            ci = v.get('CI95')
            if ci:
                ax.plot(ci, [j, j], color=color, lw=1.4)
            ax.scatter([v['rho']], [j], c=color, marker='o' if ci else 'x', s=20, zorder=3)
        ax.axvline(0, color='#999999', lw=1, ls='--')
        ax.set_xlim(-1, 1)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel('Partial rank correlation; 95% material bootstrap interval')
        ax.grid(axis='x', alpha=.2)
        for cut in [9.5, 19.5, 29.5]:
            ax.axhline(cut, color='#bbbbbb', lw=.8)
    axes[0].set_yticks(range(40), [f'{short(r)} / {short(e)}   n={coverage[r,e]["qualified_author_materials"]}' for r, e in pairs], fontsize=8)
    axes[0].set_ylim(39.8, -.8)
    fig.suptitle('GSE243013: complete adjusted reported-expansion family\nPostoperative source exploration; neither direct communication nor predictive clinical gain', fontsize=13)
    fig.subplots_adjust(left=.28, right=.98, top=.94, bottom=.075, wspace=.13)
    fig.text(.63, .015, 'o: qualified interval   |   x: point only, interval unqualified   |   NA: coverage / design unqualified', ha='center', fontsize=10)
    fig.savefig(out/'all40_adjusted_effects.png', dpi=160)
    plt.close(fig)
    status = dict(Counter(v['status'] for v in results))
    supported = [v for v in primary if v['primary40_BH'] < .05 and v['all_CI_qualified']]
    def stable_ci(r, e, kind):
        ci = lookup[r, e, kind].get('CI95')
        return ci is not None and ci[0]*ci[1] > 0
    supported_ci = [v for v in supported if all(stable_ci(v['R_group'], v['E_group'], k) for k in family.KINDS)]
    captured = {v['author_group']: v for v in measurement['all_author_group_capture']}
    target_calls = sum(captured[name]['known_TCR_expansion_calls'] for name in ['Treg_POOL', 'CD8_POOL', 'other_CD4_POOL'])
    auxiliary = read('local_statistics_v1/private_conditional80.json')
    aux_bh = [v for v in auxiliary if v['conditional80_BH'] < .05]
    aux_by = [v for v in auxiliary if v['conditional80_BY'] < .05]
    lines = ['# CM02支线：报告克隆扩增的完整来源结果', '',
        '## 本轮回答什么', '',
        '保持2026-10-14免疫图谱CM02 v1（T_whole_TME_primary_tumor）优先。本轮是GSE243013术后来源的作者Treg/CD8细类型内扩增状态解释；原作者已报告相关克隆现象，不当作新发现。源材料未审阅为独立患者，不能声称独立验证、治疗前预测、临床增量、直接通讯或原CM投影。', '',
        '## 测量与可检验分母', '',
        f"- 全部{measurement['original_TCR_rows']:,}条TCR旁表记录按原双键精确连回RNA；未匹配{measurement['unmatched_TCR_cells']}。",
        '- 作者flag与报告大小的实际边界为≥3；≥2差异全部来自大小2的记录，保留为固定敏感性。没有按病理效果改flag。',
        '- 两张表的CXCL13 CD8命名差异已按精确细胞键解决，不能把字符串缺失当该细类不存在。',
        f'- 原冻结三个池合计{target_calls:,}个已知调用；另{measurement["original_TCR_rows"]-target_calls:,}条属于本轮未选的RNA类型，未强行纳入。其他CD4背景仅原冻结六类，不能声称已控制全部剩余T状态。',
        '- 已知TCR状态的细胞作为比例分母，未知不补未扩增。克隆大小必要检查与字段语义通过不等于受体序列重建。',
        f'- 完整120分量状态：`{json.dumps(status, ensure_ascii=False)}`。增殖Treg的材料覆盖不足保留NA，不下调资格门槛。', '',
        '![完整覆盖](all40_coverage.png)', '', '## 完整调整后的结果', '',
        f'- 全部实际可拟合设计{len(fitted)}项已重放完整19999次wild与2000次bootstrap，独立R最大数值误差{max(v["max_numeric_error"] for v in audits):.3g}。工作null资格：{null_ok}。',
        '- 工作null通过只接受预定工作模型下的数值资格，不证明患者独立性。若未通过，以下P/区间仅存计算结果，不能据其支持推断。',
        f'- 主40项IUT BH<.05且三项区间资格通过：{len(supported)}；三项区间都不跨0：{len(supported_ci)}。这些仍为本来源多轮探索，不是全研究FDR或确证。',
        '- 正病理rho表示MPR/pCR方向较高；协调正负表示两类内部扩增比例的同向/反向变化。所有点值均按同一原门槛和完整背景计算。',
        '- 区间依赖作者材料独立假设；未知患者重复可能使区间偏窄。调节后减弱不证明因果中介，也不否认原整体共变的价值。', '',
        f'- 辅助80项条件病理关联：BH<.05共{len(aux_bh)}项，BY<.05共{len(aux_by)}项；按R/E分量分布为`{json.dumps(dict(Counter(v["kind"] for v in aux_bh)))}`。这是辅助状态关联，不能用它替代未成立的主双状态组合。',
        '- 多个CD8背景下重复的同一Treg信号使用重叠材料，不能算多次独立复现。全部主和辅助校正值见下表，未选择新的检验家族。', '',
        '![全部调整效应](all40_adjusted_effects.png)', '',
        '|R群|CD8群|材料n|协调rho / CI|R条件病理rho / CI|E条件病理rho / CI|主BH / BY|',
        '|---|---|---:|---|---|---|---|']
    pri = {(v['R_group'], v['E_group']): v for v in primary}
    for r, e in pairs:
        j = pri[r, e]
        lines.append('|'+ '|'.join([r, e, str(coverage[r,e]['qualified_author_materials']),
            *[cell(lookup[r,e,k]) for k in family.KINDS], f"{j['primary40_BH']:.4g} / {j['primary40_BY']:.4g}"])+ '|')
    lines += ['', '## 如何继续', '',
        '先按全部主/辅助校正、区间和留一材料结果解释，不另挑小P重开家族。若支持仅为单臂临床或共同扩增，按该范围保留，不能升格双状态特异互作。来源限定组合可以有价值；需要明确同癌种与治疗阶段、患者身份和同终点的另一来源后，才能检验重复性。', '',
        '现有作者类型空间线索与该RNA/TCR来源不是同一组患者，尚不能拼成同状态的空间临床链。TCR旁表只有alpha/beta字段，gamma-delta候选仍需TRG/TRD或独立RNA/受体证据；本轮不确认也不排除。GSE288199 TableS1按用户安排后补。', '',
        '## 技术接受范围', '',
        '完整源码/测量/合同SHA、固定位置、分母/校正/随机流已数值核。图像另需实际视查；本文件生成不自动表示视觉或科学发布接受。所有科学表图仅私有，不输出个体记录、序列、矩阵或学得权重。']
    (out/'REPORT_zh.md').write_text('\n'.join(lines)+'\n')
    receipt = dict(status='PASS_COMPLETE_SOURCE_FAMILY_NUMERICAL_RECONCILIATION_VISUAL_REVIEW_PENDING',
        observed_UTC=datetime.now(timezone.utc).isoformat(),
        result_SHA256=sha(t/'local_statistics_v1/private_all120_components.json'),
        primary_SHA256=sha(t/'local_statistics_v1/private_primary40.json'),
        terminal_SHA256=sha(t/'local_statistics_v1/private_terminal.json'),
        actual_RAM_release_SHA256=sha(t/'private_local_RNA_transport_terminal_v1.json'),
        all_frozen_slots_coverage_and_complete_adjustments_checked=True,
        actual_complete_random_stream_R_replay_checked=True, full_points_exactly_match_main=True,
        complete_source_report_SHA256=sha(out/'REPORT_zh.md'), source_scientific_values_private=True)
    (t/('private_full_source_numeric_acceptance_v'+str(a.version)+'.json')).write_text(json.dumps(receipt, indent=2)+'\n')
    print('COMPLETE_SOURCE_NUMERICAL_REVIEW_AND_FIGURES_VISUAL_REVIEW_PENDING')


if __name__ == '__main__':
    main()
