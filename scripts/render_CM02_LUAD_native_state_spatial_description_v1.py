"""Render source aggregates only; image spread is never a confidence interval."""
from pathlib import Path
import argparse
import json
import os


def main():
    p=argparse.ArgumentParser();p.add_argument('--task',type=Path,required=True)
    p.add_argument('--output-name',default='summary_figures_v1');args=p.parse_args()
    assert Path(args.output_name).name==args.output_name and args.output_name not in {'.','..'}
    source=args.task/'native_description_v1'
    terminal=json.loads((source/'private_terminal.json').read_text())
    assert terminal['status']=='COMPLETE_FROZEN_ALL11_DESCRIPTIVE_SCOPE'
    result=json.loads((source/'private_complete_source_descriptive_summary.json').read_text())
    rows=result['all27_fixed_neighborhood_pair_slots'];assert len(rows)==27
    keys={(r['k'],r['R_state'],r['E_state']) for r in rows};assert len(keys)==27
    assert {r['k'] for r in rows}=={5,10,20}
    statuses=result['source_status_counts'];assert sum(statuses.values())==11
    accepted=statuses.get('COMPLETE_ALL_FIXED_DESCRIPTIVE_SLOTS',0)
    mainrows=[r for r in rows if r['k']==10];assert len(mainrows)==9
    for r in rows:
        for key in ['R_summary','other_T','matched_R_minus_other']:
            v=r[key]
            assert 0<=v['n']<=accepted
            assert v['positive']+v['negative']<=v['n']
            assert (v['median'] is None)==(v['n']==0)
        assert r['matched_R_minus_other']['n']<=min(r['R_summary']['n'],r['other_T']['n'])
    out=args.task/args.output_name;out.mkdir(exist_ok=False)
    os.environ.setdefault('MPLCONFIGDIR',str(out/'matplotlib_config'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    counts=result['source_total_pattern_counts']
    capture_keys=['T_support','regulatory_only','CD8_only','R_TNFRSF9','R_LAG3',
                  'R_cycling','E_dysfunction','E_cycling','E_memory']
    fig,ax=plt.subplots(figsize=(10,5.5))
    values=[counts[k] for k in capture_keys]
    ax.barh(range(len(values)),values,color='#217a88',height=.65)
    ax.set_xscale('symlog',linthresh=1);ax.set_xlim(0,max(1,max(values))*3)
    for i,value in enumerate(values):ax.annotate(f'{value:,}',(value,i),xytext=(5,0),
                                               textcoords='offset points',va='center',fontsize=9)
    ax.set_yticks(range(len(values)),capture_keys);ax.invert_yaxis()
    ax.set_xlabel('Positive author cell-bin observations (symlog; zero retained)')
    ax.set_title('GSE328481: source RNA capture in all 11 measured files\n'
                 'Overlapping expression patterns; these are not validated cell identities')
    ax.spines[['top','right']].set_visible(False);fig.tight_layout()
    fig.savefig(out/'source_capture.png',dpi=180);plt.close(fig)

    def plot(key,title,filename):
        fig,ax=plt.subplots(figsize=(11,6.2))
        labels=[];limits=[0.]
        for i,row in enumerate(mainrows):
            value=row[key];labels.append(row['R_state']+' / '+row['E_state']+f"   n={value['n']}")
            if value['n']:
                lo,hi=value['IQR'];limits.extend([lo,hi,value['median']])
                ax.plot([lo,hi],[i,i],color='#217a88',lw=2)
                ax.plot(value['median'],i,'o',color='#217a88',ms=5)
            else:
                ax.text(.02,i,'Insufficient capture / no eligible images',transform=ax.get_yaxis_transform(),
                        va='center',fontsize=8,color='#777777')
        extent=max(1,max(abs(v) for v in limits))*1.15
        ax.set_xlim(-extent,extent);ax.axvline(0,color='#888888',lw=.8,ls='--')
        ax.set_yticks(range(9),labels);ax.set_ylim(-.65,8.65);ax.invert_yaxis()
        label=('Median matched difference in log edge odds ratio' if key=='matched_R_minus_other'
               else 'Median log edge odds ratio')
        ax.set_xlabel(label+' and between-image IQR (descriptive; no P or CI)')
        ax.set_title(title+'\nGSE328481 author cell bins; k=10; no independent-patient claim')
        ax.spines[['top','right']].set_visible(False);fig.tight_layout()
        fig.savefig(out/filename,dpi=180);plt.close(fig)
    plot('R_summary','All nine fixed R/E RNA patterns','all9_primary_k10_description.png')
    plot('matched_R_minus_other','R coupling minus other-T coupling on the same eligible images','all9_matched_other_T_description.png')
    fmt=lambda v:'NA' if v is None else f'{v:.4g}'
    estimable=sum(r['R_summary']['n']>0 for r in mainrows)
    matched=sum(r['matched_R_minus_other']['n']>0 for r in mainrows)
    lines=['# LUAD原生空间：完整测量与可估计性报告', '',
           '对象：2026-10-14免疫图谱CM02 v1（T_whole_TME_primary_tumor）的状态落实支线。'
           '不改原CM与正式身份，不提供新的P、FDR、患者区间或临床结果。', '',
           '## 输入与解释', '',
           f'- 全部11份冻结公开文件均完成处理尝试；完整测量{accepted}份。',
           f'- 主k10可估计R/E组合{estimable}/9；可估计匹配其他T差值的组合{matched}/9。',
           '- 测量完成不等于效应可估计；无法估计不能解释为没有细胞或没有互作。',
           '- 作者cell-bin、原始UMI、物理微米尺度、独立患者及治疗阶段均未由本轮确认。',
           '- 只使用deposited正整数的阳性模式；未套用另一scRNA来源的RNA/MT阈值。',
           '- 同一xy的全部重复观察排除。第k距离并列邻居全部纳入，超过5倍最近邻中位距离的边排除。',
           '- 正log OR表示R状态阳性锚点的邻接E父群中，E目标状态相对更多；不是通信强度、因果或预测。',
           '- 图中横线为合格图像之间的IQR，**不是置信区间**。没有用细胞/边数充当独立患者样本数。', '',
           '## 固定表达模式', '',
           'T证据：CD3D/CD3E/TRAC任一阳性。R另需FOXP3及IL2RA，E另需CD8A；R/E双阳性从两者排除。', '',
           '- R_TNFRSF9、R_LAG3：R父群内相应marker阳性。',
           '- R_cycling、E_cycling：对应父群内MKI67与TOP2A/UBE2C至少一个阳性。',
           '- E_dysfunction：E父群内PDCD1与TOX/HAVCR2至少一个阳性。',
           '- E_memory：E父群内TCF7/IL7R至少一个阳性。', '',
           '同一cell-bin可同时满足不同状态；其他T为pan-T且非R/E/双父群，不能直接称纯常规CD4。', '',
           '## 来源级捕获', '']
    if accepted:
        counts=result['source_total_pattern_counts']
        other=counts['T_support']-counts['regulatory_only']-counts['CD8_only']-counts['dual_ambiguous']
        assert other>=0
        lines += [f"成功处理部分共{result['native_author_rows']}个作者观察；重复坐标排除{result['duplicate_coordinate_rows_excluded']}个。",
                  '以下状态数均在唯一坐标观察内汇总，不能当独立患者或互斥亚群总和。', '',
                  '| 类别/模式 | 总阳性观察数 |','|---|---|']
        selected=['T_support','regulatory_only','CD8_only','dual_ambiguous',
                  'R_TNFRSF9','R_LAG3','R_cycling','E_dysfunction','E_cycling','E_memory']
        for key in selected:lines.append(f'| {key} | {counts[key]} |')
        lines.append(f'| 其他T | {other} |')
    else:
        lines += ['没有图像完成测量，不能把空汇总中的零初始化解释成生物学零。']
    lines += ['', '![来源RNA捕获](source_capture.png)', '',
              '## 全部27个预定描述位置', '',
              'n为符合父群、阳/阴状态和实际接收邻居机会门槛的图像数。'
              'R和其他T单独汇总可有不同分母；比较特异性只看最后两列的同图像差值，不能相减两个非匹配中位数。', '',
              '| k | R模式 | E模式 | R n | R中位log OR | 其他T n | 其他T中位log OR | 匹配n | 匹配差值中位 |',
              '|---|---|---|---|---|---|---|---|---|']
    for row in rows:
        r,o,d=[row[k] for k in ['R_summary','other_T','matched_R_minus_other']]
        lines.append(f"| {row['k']} | {row['R_state']} | {row['E_state']} | {r['n']} | {fmt(r['median'])} | "
                     f"{o['n']} | {fmt(o['median'])} | {d['n']} | {fmt(d['median'])} |")
    lines += ['', '## 主k10的缺资格原因', '',
              '| R模式 | E模式 | R侧原因计数 | 其他T侧原因计数 |','|---|---|---|---|']
    for r in mainrows:
        lines.append('| '+r['R_state']+' | '+r['E_state']+' | '+json.dumps(r['R_status_counts'])+' | '+json.dumps(r['other_status_counts'])+' |')
    lines += ['', '## 判断边界', '',
              '所有尺度、来源与缺测位置均完整保留。来源限定关系可有生物价值，'
              '但阳性RNA/邻近描述不能证明细胞间直接通讯、调节T特异机制或临床增量。'
              'RNA量、检出深度、空间聚集、组织区域、分割边界与未知重复仍需后续独立设计。', '',
              '本来源取材/治疗阶段与病理终点未资格，不能称NSCLC新辅助MPR或头颈癌的同终点复现。'
              '若捕获不足，应优先解决测量深度或转经参考验证的面板适配，不据此断言细胞不存在，也不直接放宽本轮门槛。', '',
              '![k10全部组合](all9_primary_k10_description.png)', '',
              '![同图像其他T差值](all9_matched_other_T_description.png)', '']
    (out/'REPORT_zh.md').write_text('\n'.join(lines))
    print(json.dumps(dict(status='COMPLETE_SOURCE_CAPTURE_REPORT_AND_THREE_FIGURES_NO_NEW_INFERENCE')))


if __name__=='__main__':main()
