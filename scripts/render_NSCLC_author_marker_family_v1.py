"""Complete source-level conditional-family report, no individual data."""
from pathlib import Path
import argparse
import json
import os


def main():
    p=argparse.ArgumentParser();p.add_argument('--task',type=Path,required=True);a=p.parse_args();task=a.task
    folder=task/'state_family_v1';read=lambda n:json.loads((folder/n).read_text())
    terminal=read('private_terminal.json');coverage=read('private_main_coverage.json')
    primary=read('private_primary30.json');coord=read('private_coordination30.json')
    conditional=read('private_conditional60.json');marginal=read('private_marginal13.json')
    references=read('private_all30_matched_other_CD4_references.json')
    context=read('private_all33_author_subtype_context_points.json');points=read('private_all751_point_diagnostics.json')
    plan=json.loads((task/'private_frozen_author_state_family_statistics_v1.json').read_text())
    assert [len(v) for v in [primary,coord,conditional,marginal,references,context,points]]==[30,30,60,13,30,33,751]
    co={(r['R_state'],r['E_state']):r for r in coord}
    cond={(r['kind'],r['R_state'],r['E_state']):r for r in conditional}
    prim={(r['R_state'],r['E_state']):r for r in primary}
    ref={(r['R_state'],r['E_state']):r for r in references}
    assert set(co)==set(prim)==set(ref)=={tuple(p) for p in plan['all30_pairs']}
    out=task/'summary_figures_v1';out.mkdir(exist_ok=False)
    os.environ.setdefault('MPLCONFIGDIR',str(out/'matplotlib_config'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    accepted=terminal['status']=='COMPLETE_AUTHOR_MARKER_FAMILY_WORKING_NULL_ACCEPTED'
    def color(pair):
        row=prim[pair]
        return '#147c89' if accepted and row['primary30_BH']<.05 and row['coordination_CI_qualified'] and row['regulatory_addition_CI_qualified'] else '#888888'
    def dot(ax,row,i,tone,value='rho',interval='CI95'):
        if row.get(value) is None:
            ax.text(-.95,i,'NA',va='center',fontsize=8,color='#777777');return
        ax.plot(row[value],i,'o',color=tone,ms=4)
        if row.get(interval):ax.plot(row[interval],[i,i],lw=1.4,color=tone)
    fig,axes=plt.subplots(3,2,figsize=(15,13),sharex=True)
    for i,r in enumerate(plan['R_states']):
        for j,e in enumerate(plan['E_states']):
            pair=(r,e);dot(axes[i,0],co[pair],j,color(pair));dot(axes[i,1],cond[('conditional_R',r,e)],j,color(pair))
        for k in [0,1]:
            ax=axes[i,k];ax.set_xlim(-1,1);ax.axvline(0,color='#999999',lw=.7,ls='--')
            ax.set_yticks(range(10),plan['E_states'] if k==0 else []);ax.invert_yaxis();ax.spines[['top','right']].set_visible(False)
            ax.set_title(r+'\n'+('R/E coordination' if k==0 else 'R pathology association conditional on E'),fontsize=10)
            ax.tick_params(axis='y',labelsize=8)
            ax.set_xlabel('Partial rank correlation; 95% author-case bootstrap CI',fontsize=8)
    fig.suptitle('Complete exploratory source family — GSE243013\nTeal: primary IUT BH30<0.05 and CI qualified; gray otherwise'
                 +('' if accepted else '\nWORKING NULL GATE NOT ACCEPTED — no inferential priority'),fontsize=12)
    fig.tight_layout(rect=(0,0,1,.94));fig.savefig(out/'all30_coordination_and_regulatory_addition.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(16,5.5),sharey=True)
    for i,r in enumerate(plan['R_states']):
        for j,e in enumerate(plan['E_states']):dot(axes[i],ref[(r,e)],j,'#356c8a',value='rho_difference')
        ax=axes[i];ax.axvline(0,color='#999999',lw=.7,ls='--');ax.set_xlim(-1,1)
        ax.set_yticks(range(10),plan['E_states']);ax.set_title(r,fontsize=10)
        ax.set_xlabel('Matched R minus other-CD4 rho\nPointwise author-case bootstrap 95% CI',fontsize=8)
        ax.spines[['top','right']].set_visible(False);ax.tick_params(axis='y',labelsize=8)
    axes[0].invert_yaxis();fig.suptitle('Same-case, same-marker non-Treg CD4 reference — no simultaneous specificity claim',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.92));fig.savefig(out/'all30_matched_other_CD4_differences.png',dpi=180);plt.close(fig)
    fmt=lambda x:'NA' if x is None else f'{x:.4g}'
    ci=lambda x:'NA' if x is None else '['+', '.join(fmt(v) for v in x)+']'
    lines=['# 作者Treg内marker与CD8亚型：完整来源结果', '',
        '范围：2026-10-14免疫图谱CM02 v1（T_whole_TME_primary_tumor）的来源状态落实支线。'
        '这是同一已探索来源的新测量问题，原RNA父群家族保持，不是独立验证或全癌图谱身份更新。', '',
        '## 主分析对象与问题', '',f"- 实际终态：`{terminal['status']}`。",
        f"- 同QC双作者父群捕获：{coverage['clinical_and_both_author_parent_capture']}个材料；药物身份覆盖后保留"
        f"{coverage['retained_after_agent_group10']}，MPR/pCR={coverage['response1']}、non-MPR={coverage['response0']}。",
        '- 单位是作者材料，不能当已审独立患者。',
        '- R：先采用作者Treg身份，再测TNFRSF9、LAG3、MKI67各自≥1原始计数的比例。MKI67_pos不是旧多基因联合定义。',
        '- E：全部作者CD8T细亚型的父群内比例，保留作者MAIT/NK-like等细类。分母是该作者CD8T谱系，细类比例有闭合约束。',
        '- 问题：在CD8状态相近的材料中，Treg marker差异是否还随病理标签变化？这是条件rank关联，不是ΔAUC/Brier或治疗前预测。', '',
        '主P=max（协调P，R病理关联控制配对E后的P），完整BH30/BY30。正临床rho表示MPR/pCR方向更高。'
        '103主分量、751点值和所有33作者亚型背景完整保留，缺资格不补阴性。', '',
        '## 完整30配对', '',
        '| Treg marker状态 | 作者CD8T亚型 | 协调rho [CI] | R给定E的临床rho [CI] | 主交并P | BH30 / BY30 |',
        '|---|---|---|---|---|---|']
    for r in primary:
        pair=(r['R_state'],r['E_state']);c=co[pair];x=cond[('conditional_R',*pair)]
        lines.append(f"| {pair[0]} | {pair[1]} | {fmt(c.get('rho'))} {ci(c.get('CI95'))} | "
                     f"{fmt(x.get('rho'))} {ci(x.get('CI95'))} | {fmt(r['regulatory_addition_IUT_P'])} | "
                     f"{fmt(r['primary30_BH'])} / {fmt(r['primary30_BY'])} |")
    lines += ['', '## 全13临床边际（辅助家族）', '',
              '| 状态 | rho | 95%材料bootstrap区间 | P | BH13 / BY13 |','|---|---|---|---|---|']
    for r in marginal:
        lines.append(f"| {r['state']} | {fmt(r.get('rho'))} | {ci(r.get('CI95'))} | {fmt(r.get('raw_P'))} | "
                     f"{fmt(r['marginal13_BH'])} / {fmt(r['marginal13_BY'])} |")
    lines += ['', '## 同材料、同marker的其他CD4参照', '',
              '这里的其他CD4是作者非Treg CD4T集合，不把所有未知T混入。匹配设计同时调整三父群质量/比例。'
              '差值CI逐对描述，没有全30同时覆盖或特异性P保证，不能相减两条非匹配关系推特异性。', '',
              '| R状态 | E亚型 | 匹配单位数 | R rho | 非Treg CD4 rho | 差值 [95%CI] |','|---|---|---|---|---|---|']
    for r in references:
        lines.append(f"| {r['R_state']} | {r['E_state']} | {r['n']} | {fmt(r.get('R_rho'))} | "
                     f"{fmt(r.get('other_CD4T_rho'))} | {fmt(r.get('rho_difference'))} {ci(r.get('CI95'))} |")
    lines += ['', '## 其余完整证据与边界', '',
        '[全部双向条件临床辅助结果](../state_family_v1/private_conditional60.json)、'
        '[751点值诊断](../state_family_v1/private_all751_point_diagnostics.json)、'
        '[33作者亚型背景点值](../state_family_v1/private_all33_author_subtype_context_points.json)、'
        '[全部实际随机流R接受](../state_family_v1/private_all_real_random_stream_R_audits.json)。', '',
        '所有主要关系共同材料、完整临床/药物与双方质量列空间；未调整/仅临床/全保留、父群10/30和肺癌亚型点值并列。'
        '条件分量在“无背景”设置中仍控制配对状态。加术后病理的协调仅诊断，不能证明因果。'
        '全作者保留与QC相同来自实测集合相同，不是额外独立复现。', '',
        '本来源已多轮探索，作者亚型也不是功能真值。当前定义与旧FOXP3/IL2RA双阳性父群不同，不能把相同marker名当同一定义的复现。'
        '原生空间测量亦需核父群/状态一致、组织/取样/病例资格，才可连接成证据链。'
        '本轮未证明直接LR、耐药机制、预测增量、来源间同终点复现或原CM02活动本身。', '',
        '![完整30配对](all30_coordination_and_regulatory_addition.png)', '',
        '![全部匹配参照](all30_matched_other_CD4_differences.png)', '']
    (out/'REPORT_zh.md').write_text('\n'.join(lines))
    print(json.dumps(dict(status='COMPLETE_AUTHOR_LINEAGE_MARKER_SOURCE_REPORT_AND_TWO_FIGURES')))


if __name__=='__main__':main()
