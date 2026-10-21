"""Private source-level forest and report, never individual rows or arrays."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def stable_interval(row):
    ci=row.get('case_bootstrap_95_CI')
    loo=row.get('leave_one_case_effect_range')
    return bool(ci and loo and ((ci[0]>0 and loo[0]>0) or (ci[1]<0 and loo[1]<0)))


def grade(main,sensitivity):
    if main.get('P') is None:return '覆盖不足，不可估'
    if main['source_case_units']<10:return '病例较少，区间未报告'
    if main['BH']>.05:return '主条件检验未支持'
    if not stable_interval(main):return '条件重排显著，病例间方向未稳定'
    if sensitivity.get('P') is None or sensitivity['source_case_units']<10:return '主结果候选，区域背景稳健性未确定'
    if sensitivity['BH']<=.05 and stable_interval(sensitivity) and main['effect']*sensitivity['effect']>0:
        return '本来源内稳健组织候选'
    return '主结果候选，对区域背景敏感'


def identity(row):
    return tuple(row[key] for key in ['regulatory_context','left_author_label','right_author_label','left_marker','right_marker'])


def number(value):return 'NA' if value is None else f'{value:.6g}'


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--result',type=Path,required=True)
    parser.add_argument('--wrapper',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data=json.loads(args.result.read_text())
    assert data['status']=='COMPLETE_CONDITIONAL_CONTINUOUS_STATE_ORGANIZATION_NOT_CLINICAL_VALIDATION'
    assert data['wrapper_SHA256']==hashlib.sha256(args.wrapper.read_bytes()).hexdigest()
    assert data['independent_numeric_replay']['all_checks_passed']
    assert data['family_size']==18 and data['permutations']==data['case_bootstraps']==2000
    first=data['family']['ROI_parent_shuffle'];second=data['family']['ROI_parent_grid4_shuffle']
    assert len(first)==len(second)==18 and len(set(map(identity,first)))==18
    assert list(map(identity,first))==list(map(identity,second))
    for row in first+second:
        assert row['source_case_units']>=0
        for name in ['effect','P','BH']:
            assert row.get(name) is None or np.isfinite(row[name])
    args.output.mkdir(mode=0o700,parents=False,exist_ok=False)
    labels=[('RegT' if row['regulatory_context']=='author_regulatory' else 'CD4 control')+
        ': '+row['left_marker']+' / CD8 '+row['right_marker'] for row in first]
    fig,ax=plt.subplots(figsize=(11,9.6))
    for rows,offset,color,marker,label in [(first,-.12,'#232323','o','Within-ROI shuffle (main)'),
        (second,.12,'#287a9b','s','Within-grid shuffle (sensitivity)')]:
        legend_done=False
        for index,row in enumerate(rows):
            y=index+offset
            if row.get('effect') is None:continue
            value=row['effect'];ci=row.get('case_bootstrap_95_CI')
            caption=label if not legend_done else None
            if ci:
                # Draw the CI directly; the point need not lie inside a
                # bootstrap interval in every possible finite sample.
                ax.hlines(y,ci[0],ci[1],color=color,lw=1.25)
            ax.plot(value,y,marker=marker,color=color,markersize=4.8,
                markerfacecolor=color if ci else 'white',label=caption)
            legend_done=True
    ax.axvline(0,color='#888888',lw=.8)
    ax.axhline(8.5,color='#c9c9c9',lw=.8,ls='--')
    ax.set_yticks(range(18),labels);ax.invert_yaxis()
    for index,(a,b) in enumerate(zip(first,second)):
        ax.text(1.02,index,f"n={a['source_case_units']}/{b['source_case_units']}",
            transform=ax.get_yaxis_transform(),va='center',fontsize=8)
    ax.set_xlabel('Conditional rank covariance (observed minus null expectation)')
    ax.set_title('2026-10-14 immune-atlas CM02 v1: protein-state exploration\nPretreatment melanoma, invasive margin; native k10',fontsize=11)
    handles,_=ax.get_legend_handles_labels()
    if handles:ax.legend(loc='best',fontsize=8,frameon=False)
    ax.spines[['top','right']].set_visible(False)
    fig.subplots_adjust(left=.35,right=.85,bottom=.13,top=.91)
    fig.text(.35,.035,'95% case-bootstrap intervals; open points: interval suppressed.\nConditional tests do not establish clinical value or direct ligand-receptor signaling.',fontsize=8)
    fig.savefig(args.output/'source_state_forest.png',dpi=180)
    fig.savefig(args.output/'source_state_forest.svg')
    plt.close(fig)
    grades=[grade(a,b) for a,b in zip(first,second)]
    counts={key:grades.count(key) for key in dict.fromkeys(grades)}
    text=['# 2026-10-14免疫图谱CM02（v1，T_whole_TME_primary_tumor）启发的蛋白状态探索','',
        '当前结果为治疗前黑色素瘤浸润边缘的蛋白探索；没有投影原CM02权重，不是原CM02临床验证。',
        '图定义：'+data['graph']+'。该图与未完整读取的作者距离图不同；局部几何前缀检查不证明全图或微米校准。',
        '所有18组合均保留；主BH18与区域网格敏感性BH18分别报告。区间按病例单位bootstrap；P/BH来自条件细胞状态重排，不是总体临床关联检验。','',
        '|左端|左marker|CD8 marker|主n|主效应|病例95%区间|主BH|网格n|网格效应|网格BH|证据分组|',
        '|---|---|---|---:|---:|---|---:|---:|---:|---:|---|']
    for a,b,g in zip(first,second,grades):
        ci=a.get('case_bootstrap_95_CI');interval='未报告' if not ci else '['+', '.join(map(number,ci))+']'
        text.append('|'+ '|'.join(['作者调节T' if a['regulatory_context']=='author_regulatory' else '常规CD4对照',
            a['left_marker'],a['right_marker'],str(a['source_case_units']),number(a.get('effect')),interval,
            number(a.get('BH')),str(b['source_case_units']),number(b.get('effect')),number(b.get('BH')),g])+'|')
    text+=['','## 解读边界','',
        '- 正协方差同时包含高—高和低—低；负协方差表示邻接端状态反向，不等于两类细胞空间排斥。',
        '- 固定细胞身份、数量和邻接，只在同ROI/身份内重排联合marker向量；网格分析增加局部区域条件。',
        '- 网格敏感性减弱可能涉及区域梯度或网格内可交换信息减少，不能仅凭减弱断言混杂解释已成立。',
        '- 一类通过检验、另一类未通过，不是调节T与常规CD4差异的直接检验。',
        '- marker为相对蛋白强度，不是经校准的阳性细胞比例或新正式细胞身份；marker组合不是配体—受体轴。',
        '- “本来源内稳健组织候选”要求主/敏感性均BH≤0.05、病例区间及留一范围不跨零、方向相同；不是独立来源复现。',
        '- 几何前缀QC、完整输入与临床覆盖应另审。此计算没有读取疗效；临床最小病例与完整设计门槛保持。','',
        '![来源级森林图](source_state_forest.png)','',
        '结果SHA256：'+hashlib.sha256(args.result.read_bytes()).hexdigest()]
    (args.output/'source_report.md').write_text('\n'.join(text)+'\n')
    print(json.dumps(dict(status='PRIVATE_SOURCE_REPORT_AND_FOREST_WRITTEN',group_counts=counts,all18_preserved=True,no_individual_records=True),ensure_ascii=False))


if __name__=='__main__':main()
