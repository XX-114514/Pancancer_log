"""Render complete source aggregate evidence; no individual inputs or arrays."""
from pathlib import Path
import argparse
import json
import os


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--task', type=Path, required=True)
    args = p.parse_args()
    source = args.task / 'state_statistics_v1'
    read = lambda name: json.loads((source / name).read_text())
    terminal = read('private_terminal.json')
    plan = json.loads((args.task / 'private_frozen_NSCLC_state_statistics_v1.json').read_text())
    synthetic = bool(plan.get('synthetic_only'))
    clinical = read('private_all6_clinical_state_associations.json')
    coordination = read('private_pair_state_coordination.json')
    joint = read('private_joint_coordination_clinical_IUT.json')
    coverage = read('private_main_same_case_coverage.json')
    reference = read('private_matched_other_T_differences.json')
    sensitivity = read('private_point_sensitivities.json')
    audit = read('private_all_real_design_R_random_stream_audit.json')
    calibration = read('private_all_fixed_working_null_diagnostics.json')
    assert len(clinical) == 6 and len(coordination) == len(joint) == len(reference) == 9
    assert len(sensitivity) == 129 and all(x['primary9_BH'] >= 0 for x in joint)
    out = args.task / 'summary_figures_v1'
    out.mkdir(exist_ok=False)
    os.environ.setdefault('MPLCONFIGDIR', str(out / 'matplotlib_config'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    def forest(rows, title, filename):
        fig, ax = plt.subplots(figsize=(10.5, max(4, .52 * len(rows) + 1.6)))
        for i, row in enumerate(rows):
            value = row.get('rho')
            interval = row.get('CI95')
            if value is not None:
                ax.plot(value, i, 'o', color='#176c87')
                if interval:
                    ax.plot(interval, [i, i], color='#176c87', linewidth=1.7)
            else:
                ax.text(0, i, 'Not estimable', color='#666666', va='center')
        labels = [x.get('state', x.get('R_state', '') + ' / ' + x.get('E_state', '')) for x in rows]
        ax.set_yticks(range(len(rows)), labels)
        ax.invert_yaxis()
        ax.axvline(0, color='#888888', linestyle='--', linewidth=.8)
        ax.set_xlim(-1.02, 1.02)
        ax.set_xlabel('Adjusted partial rank correlation; 95% author-case bootstrap CI')
        ax.set_title(('SYNTHETIC SOFTWARE CHECK ONLY\n' if synthetic else '') + title)
        ax.spines[['right', 'top']].set_visible(False)
        fig.tight_layout()
        fig.savefig(out / filename, dpi=180)
        plt.close(fig)

    forest(coordination, 'All nine observed R/E combinations — native NSCLC source', 'all9_coordination_forest.png')
    forest(clinical, 'All six observed patterns — MPR/pCR versus non-MPR', 'all6_clinical_forest.png')
    fmt = lambda x: 'NA' if x is None else f'{x:.4g}'
    interval = lambda x: 'NA' if not x else '[' + ', '.join(fmt(y) for y in x) + ']'
    lines = ['# NSCLC：固定T表达模式的完整来源级结果', '',
             '**仅合成软件检查，不是真实数据结果。**' if synthetic else '以下为真实来源级汇总。', '',
             '范围：2026-10-14免疫图谱CM02 v1（T_whole_TME_primary_tumor）的状态落实支线。'
             '这里测量固定marker的RNA观测模式，不改原CM或正式身份。', '',
             '## 执行与分母', '', f"- 统计终态：`{terminal['status']}`。",
             f"- 严格临床并双父群捕获：{coverage['clinical_and_both_parent_capture']}个作者材料单位；"
             f"具体药物类别至少10单位后保留{coverage['retained_after_agent_group10']}；"
             f"MPR/pCR={coverage['response1']}，non-MPR={coverage['response0']}。",
             '- 不能将作者材料单位当成已审独立患者。病理二元标签与HNSCC连续退缩或生存不是同终点。',
             f"- 真实设计R复算块数={len(audit)}；工作null检查设计数={len(calibration)}。", '',
             '## 九个组合：协调与双方临床关联的交并检验', '',
             '每行主要P为三个分量P的最大值；统一BH9/BY9。协调辅助BH不能替代主要交并结论。', '',
             '| R模式 | E模式 | 协调rho [95%CI] | 协调P / BH9 | 三分量交并P / BH9 / BY9 |',
             '|---|---|---|---|---|']
    for c, j in zip(coordination, joint):
        assert (c['R_state'], c['E_state']) == (j['R_state'], j['E_state'])
        lines.append(f"| {c['R_state']} | {c['E_state']} | {fmt(c.get('rho'))} {interval(c.get('CI95'))} | "
                     f"{fmt(c.get('raw_P'))} / {fmt(c['coordination9_BH'])} | "
                     f"{fmt(j['triple_IUT_P'])} / {fmt(j['primary9_BH'])} / {fmt(j['primary9_BY'])} |")
    lines += ['', '## 六种模式与病理标签', '',
              '正rho表示该模式的父群内比例在MPR/pCR方向更高；不是获益概率或预测准确率。', '',
              '| RNA模式 | rho | 95%区间 | P | 辅助BH6 / BY6 |', '|---|---|---|---|---|']
    for r in clinical:
        lines.append(f"| {r['state']} | {fmt(r.get('rho'))} | {interval(r.get('CI95'))} | "
                     f"{fmt(r.get('raw_P'))} | {fmt(r['auxiliary6_BH'])} / {fmt(r['auxiliary6_BY'])} |")
    lines += ['', '## 匹配的其他T参照', '',
              '两条关系使用同一病例与全部三父群控制，区间是共同重采样的差值区间；不由两个P之差推特异性。', '',
              '这些是逐组合的描述区间，没有提供九组合同时覆盖或特异性显著性的保证。', '',
              '| R模式 | E模式 | 共同单位数 | R-E减其他T-E | 差值95%区间 |', '|---|---|---|---|---|']
    for r in reference:
        lines.append(f"| {r['R_state']} | {r['E_state']} | {r['n']} | "
                     f"{fmt(r.get('delta_R_minus_other'))} | {interval(r.get('delta_CI95'))} |")
    lines += ['', '## 解释边界与下一步', '',
              '所有九组合和六临床边际都报告；没有按本轮显著性换marker、药物、终点或检验。'
              '129项固定敏感性只提供点值，完整表另存于本目录的上级统计目录，不形成替代P池。', '',
              '来源/癌种特有限定信号可以保留，但与外部来源相反不能叫同向复现。'
              '两端协调加病理关联仍不能证明直接通讯、调节相关T特异性或独有临床增量。'
              '术后材料不能叫治疗前预测；本来源已被项目探索过，不是未触及验证。', '',
              '后续按完整家族结果选择明确状态和同终点队列；组织定位需原生可测且不按目标marker选区的空间来源。'
              'GSE288199 Table S1依用户安排待补。', '',
              '![九组合协调](all9_coordination_forest.png)', '',
              '![六模式临床](all6_clinical_forest.png)', '']
    (out / 'REPORT_zh.md').write_text('\n'.join(lines))
    print(json.dumps({'status': 'COMPLETE_SOURCE_AGGREGATE_REPORT_AND_TWO_FIGURES',
                      'scientific_rows_published': False}))


if __name__ == '__main__':
    main()
