"""Private data-driven source-summary report; no embedded scientific results."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

METRICS=['continuous_rank_covariance','binary_standardized_covariance','binary_double_positive_edge_fraction']
NAMES=['Continuous Ki67','Binary Ki67: covariance','Binary Ki67: double-positive edges']
RELATIONS=['regulatory_CD8','conventional_CD4_CD8']
COLORS=['#ba5b25','#2878a5']


def main():
    p=argparse.ArgumentParser();p.add_argument('--summary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--mode',choices=['matched','nuclear'],default='matched');a=p.parse_args()
    nuclear=a.mode=='nuclear'
    r=json.loads(a.summary.read_text());assert r['status'].startswith('COMPLETE_DECLARED_NUCLEAR_' if nuclear else 'COMPLETE_MATCHED_')
    rows=r['results']['family'];assert len(rows)==12
    metrics=['raw_Ki67_rank_covariance','nuclear_residual_Ki67_covariance','DNA_proxy_covariance'] if nuclear else METRICS
    names=['Raw Ki67 ranks','Ki67 nuclear residual','DNA proxy'] if nuclear else NAMES
    a.output.mkdir(mode=0o700,exist_ok=False)
    fig,axes=plt.subplots(1,3,figsize=(13,4.6))
    for ax,metric,title in zip(axes,metrics,names):
        scale=100 if metric=='binary_double_positive_edge_fraction' else 1
        for relation,color in zip(RELATIONS,COLORS):
            selected=[x for x in rows if x['parent_relation']==relation and x['metric']==metric]
            assert len(selected)==2
            points=[(i,x['effect']*scale,x) for i,x in enumerate(selected) if x['effect'] is not None]
            if points:
                ax.plot([x[0] for x in points],[x[1] for x in points],'-o',color=color,
                        label='Regulatory-CD8' if relation==RELATIONS[0] else 'Conventional CD4-CD8')
                for x,y,row in points:
                    ax.annotate(f"q={row['BH']:.3g}",(x,y),xytext=(7,5),textcoords='offset points',fontsize=8,color=color)
        ax.axhline(0,color='#999999',linewidth=.8,zorder=-1)
        ax.set_xticks([0,1],['Grid','Grid + tumor\nneighbor'])
        ax.set_xlim(-.15,1.7);ax.margins(y=.32)
        ax.set_title(title,fontsize=10)
        ax.set_ylabel('Observed minus null: percentage points' if scale==100 else 'Observed minus null: standardized covariance',fontsize=8)
        ax.spines[['top','right']].set_visible(False)
    handles,labels=axes[0].get_legend_handles_labels()
    if handles:fig.legend(handles,labels,loc='upper center',ncol=2,frameon=False,bbox_to_anchor=(.5,.97))
    fig.suptitle('Matched-material nuclear signal sensitivity' if nuclear else 'Matched-material representation and spatial-null sensitivity',fontsize=13,y=1.02)
    fig.text(.5,-.025,'Material-unit exploration, patient independence unreviewed; no patient CIs. BH over all 12 positions.\nDifferent metric units are not pooled; lines describe changes, not formal differences or causal effects.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.86))
    for suffix in ['png','svg']:fig.savefig(a.output/f'matched_representation_null.{suffix}',dpi=180,bbox_inches='tight')
    plt.close(fig)
    description=('Ki67秩对预设Histone H3与DNA代理回归，残差按单位RMS标准化；DNA代理为两个Ir通道单位RMS秩的等均值再标准化。所有声明通道必需、设计不满秩/残差无变异时不可估；不自动删项。' if nuclear else
                 '连续与二值指标先在父类内平均秩/单位RMS。二值HH边比例另作指标衔接，单位不同。')
    lines=['# '+('同材料核信号敏感性检查' if nuclear else '同材料表示与空间条件交叉检查'),'','主线：2026-10-14免疫图谱CM02（v1，T_whole_TME_primary_tumor）启发的独立蛋白状态问题。',
           f"完成UTC：{r['completed_UTC']}。",'',
           '每条父类关系的六项检查使用完全相同材料、父类和边集；'+description+'BH同时覆盖12项。材料独立患者身份未核，不提供患者置信区间或临床模型。',
           '', '## 来源级结果','','|关系|指标|空间条件|共同材料|效果|BH|正/负/零材料|','|---|---|---|---:|---:|---:|---:|']
    for row in rows:
        effect='不可估' if row['effect'] is None else f"{row['effect']:+.5f}"
        q='未检验' if row['BH'] is None else f"{row['BH']:.5g}"
        signs='—' if row['effect'] is None else '/'.join(str(row[x]) for x in ['positive_materials','negative_materials','zero_materials'])
        lines.append(f"|{row['parent_relation']}|{row['metric']}|{row['null']}|{row['common_source_materials']}|{effect}|{q}|{signs}|")
    lines+=['','效果为观测减该空间条件下零假设期望；所有协方差及来源置换均值有独立标量复算。跨位置P/效果变化仅作描述，不是效应差检验，也不能证明某种条件造成了变化。','',
            '## 核/背景资格（在效果前保存）','']
    if nuclear:
        lines+=['|父类|设计资格原因计数|','|---|---|']
        for role in r['preeffect_nuclear_design_qualification']['roles']:
            lines.append(f"|{role['role']}|{json.dumps(role['reason_counts'],ensure_ascii=False)}|")
    else:
        lines+=['|父类|通道|有变异材料|恒定材料|非有限材料|','|---|---|---:|---:|---:|']
        for role in r['preeffect_background_qualification']['roles']:
            for channel in role['per_channel_source_material_counts']:
                count=channel['counts'];lines.append(f"|{role['role']}|{channel['field']}|{count.get('variable',0)}|{count.get('constant',0)}|{count.get('nonfinite',0)}|")
    lines+=['',('本轮另行固定核信号模型；旧恒定控制不可估家族保留，拟合系数仅RAM。' if nuclear else '本轮只核通道可用性，未拟合新背景效果。')+'不自动删通道。核/背景信号可能含真实周期生物学，也不是完全技术性负对照；衰减不证明伪影或无增殖。','',
            '## 解释范围','','- 原三、八、十八位置家族保留，不以新BH覆盖旧结果，不选择最小P当主线。',
            ('- 核调整是另一观察量，不是独立病例验证；不以新家族最小P取代原关系。' if nuclear else '- 阳性概率阈值固定≥0.5；严格>0.5仅核输入变异/阈值等值数量，没有根据效果改阈值。'),
            ('- 连续同向关系含高高和低低；核残差关系不直接对应Ki67阳性组合。' if nuclear else '- 连续同向关系含高高和低低，二值协方差同样不能证明仅高高组合主导；HH指标检验双阳性边富集。'),
            '- 原生k10一跳肿瘤近邻不等于微米接触、病理侵袭边缘、TLS或免疫突触；局部自相关和分割/共同通道背景仍有限制。',
            '- 本次没有重新核实患者—材料映射、临床阶段或终点，不作同终点独立临床复现。原CM02权重未投影。',
            '- Ki67较低并不能排除周期活动；[多指标质谱原始研究](https://www.nature.com/articles/s41467-025-63883-4)支持多指标解释，但不是本队列细胞周期状态的验证。','',
            '[公开输入和作者方法](https://zenodo.org/records/5903179)。本报告数值来自本轮计算，不是论文原结果。']
    (a.output/'REPORT_zh.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':main()
