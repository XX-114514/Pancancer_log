"""Render the complete fixed source-aggregate family, never individual arrays."""
from pathlib import Path
import argparse,hashlib,json,os


def number(x):return '—' if x is None else f'{x:.4g}'


def label(key):
    parent,source,target=key.split('__')
    parent={'RegMixed':'Mixed Reg/Tex','CD4control':'CD4 control','RegMinusCD4':'Reg/Tex minus CD4'}[parent]
    return parent+' | '+source.replace('_reference',' ref').replace('_','/')+' -> '+target.replace('CD8_','CD8 ')


def main(source,outdir):
    r=json.loads(source.read_text());assert r['status']=='COMPLETE_FIXED_CONTINUOUS_MARKER_SPATIAL_CASE_FAMILY'
    assert not outdir.exists();outdir.mkdir()
    os.environ.setdefault('MPLCONFIGDIR',str(outdir/'plot_cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    primary=r['results']['primary25_grid2']
    main_keys=[k for k in primary if 'Ki67_reference' not in k]
    reference_keys=[k for k in primary if 'Ki67_reference' in k]
    assert len(main_keys)==12 and len(reference_keys)==4
    order=main_keys+reference_keys
    lines=['# 固定连续标记空间配对结果','',
      '原2026-10-14免疫图谱CM02 v1保持优先；本研究是独立乳腺来源中的连续蛋白标记空间比较，不是原RNA CM02投影。',
      '', '## 数据与检验范围','',
      f"- 临床连接且混合Reg/Tex、CD8各至少20细胞：{r['source_scope']['clinical_linked_target_eligible_cases']}病例。",
      f"- 其中纯CD4比较父类也至少20细胞：{r['source_scope']['control_parent_also_at_least20_cases']}病例。",
      '- 全成员CRC/SHA、原资格计数/病例对齐已复核；只保存来源级结果，原始payload/个体数组未落盘。',
      '- 主25微米/2x2网格，固定核信号/面积及非Ki67轴Ki67背景；2000次ROI内联合残差向量置换仅定义空间参照均值；2000次病例signflip检验，2000病例bootstrap。',
      '- 固定12主槽BH，四个一般Ki67源轴只作参照；非主设置和ER/作者标志子组无新增P。',
      '- 全局类型映射未完成；混合Reg/Tex不等于纯Treg，PD1不等于功能耗竭。','',
      '## 完整主家族及一般增殖参照','',
      '|比较|病例数|空间超额相关 Δr|病例bootstrap 95%区间|病例signflip P|BH12|BY12敏感性|状态|',
      '|---|---:|---:|---|---:|---:|---:|---|']
    robust=[]
    for key in order:
        v=primary[key]['ALL_ELIGIBLE'];ci=v.get('CI95');cis='—' if ci is None else f'[{number(ci[0])}, {number(ci[1])}]'
        lines.append('|'+ '|'.join([label(key),str(v['n_cases']),number(v.get('excess_r')),cis,
          number(v.get('patient_signflip_P')),number(v.get('BH_fixed_family')),number(v.get('BY_arbitrary_dependence_sensitivity')),v['status']])+'|')
        if key in main_keys and v.get('BH_fixed_family',1)<=.05 and ci is not None and ci[0]*ci[1]>0:
            robust.append(key)
    lines+=['','## 有限证据判断','',
      f'同时通过主BH12与病例区间不跨零的槽数：{len(robust)}。这仅是本来源条件空间关联，尚无独立来源或临床增量验证。',
      '按原固定设计保留全部方向与不可估计槽；参照/子组/敏感性不用于另选显著结果。以下三个设置只作效应及区间诊断。','',
      '## 空间敏感性完整清单','',
      '|设置|比较|病例数|Δr|95%区间|','|---|---|---:|---:|---|']
    for config,summary in r['results'].items():
        if config=='primary25_grid2':continue
        assert len(summary)==16
        for key in order:
            v=summary[key]['ALL_ELIGIBLE'];assert 'patient_signflip_P' not in v
            ci=v.get('CI95');cis='—' if ci is None else f'[{number(ci[0])}, {number(ci[1])}]'
            lines.append('|'+ '|'.join([config,label(key),str(v['n_cases']),number(v.get('excess_r')),cis])+'|')
    lines+=['','## 解释边界','',
      'Δr是观察的异类邻近边标记相关减去同ROI/父类/宏观区域置换的均值；图、数量和位置固定，所以不等于双方细胞总数相关。细胞置换的原P已因模拟误报膨胀弃用；新推断以病例整体signflip为单位，仍要求病例独立及null效应对称。空间参照均值及测量背景仍有模型假设，不能证明直接接触、LR通讯或因果。',
      'Reg-minus-CD4只在同病例同ROI交集比较；需要该差值、各自效应和敏感性共同支持，才能讨论相对CD4的差异，不据Reg一项显著断言Treg特异。',
      '本来源已有相关增殖Treg/预后文献。临床连接/预留样本覆盖受限，本轮没有临床模型，也不称新的耐药机制或独立临床验证。','',
      '## 计算证据','',f"- 真实运行：{r['started_UTC']} 至 {r['completed_UTC']} UTC。",
      f"- 数值审计：{json.dumps(r['independent_real_positive_negative_and_null_replay'],ensure_ascii=False)}。",
      f"- 源结果SHA256：{hashlib.sha256(source.read_bytes()).hexdigest()}。"]
    (outdir/'REPORT_zh.md').write_text('\n'.join(lines)+'\n')
    fig,ax=plt.subplots(figsize=(13,9.8));colors=['#146c94' if k in main_keys else '#777777' for k in order]
    for i,key in enumerate(order):
        v=primary[key]['ALL_ELIGIBLE'];effect=v.get('excess_r');ci=v.get('CI95');y=len(order)-1-i
        if effect is not None:
            ax.plot(effect,y,'o',color=colors[i],ms=5)
            if ci is not None:ax.plot(ci,[y,y],color=colors[i],lw=1.5)
        text=f"n={v['n_cases']}"
        if 'BH_fixed_family' in v:text+=f"; q={v['BH_fixed_family']:.3g}"
        elif key in main_keys:text+='; insufficient coverage'
        else:text+='; reference, no P'
        ax.text(1.02,y,text,va='center',fontsize=8,transform=ax.get_yaxis_transform())
    ax.axvline(0,color='black',lw=.8);ax.axhline(3.5,color='#bbbbbb',lw=.7)
    ax.set_yticks(range(len(order)),[label(k) for k in order[::-1]],fontsize=8)
    ax.set_xlabel('Excess edge correlation (observed minus conditional permutation mean)')
    ax.set_title('Primary spatial family and general Ki67 references\n25 um; measured-background adjustment; patient bootstrap 95% CI',fontsize=11)
    ax.grid(axis='x',alpha=.2);fig.subplots_adjust(left=.41,right=.78,bottom=.09,top=.91)
    fig.savefig(outdir/'primary_and_reference_forest.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(4,4,figsize=(16,11),sharex=False)
    configs=list(r['results']);names=['25um / grid2','50um / grid2','25um / grid4','25um / grid2 + tissue blocks']
    for ax,key in zip(axes.ravel(),order):
        for y,(config,name) in enumerate(zip(configs,names)):
            v=r['results'][config][key]['ALL_ELIGIBLE'];effect=v.get('excess_r');ci=v.get('CI95')
            if effect is not None:
                ax.plot(effect,3-y,'o',ms=4,color=colors[order.index(key)])
                if ci is not None:ax.plot(ci,[3-y,3-y],lw=1,color=colors[order.index(key)])
            ax.text(.98,3-y,f"n={v['n_cases']}",ha='right',va='center',fontsize=6,transform=ax.get_yaxis_transform())
        ax.set_title(label(key).replace(' | ','\n'),fontsize=8);ax.axvline(0,color='black',lw=.6)
        ax.set_yticks(range(4),names[::-1],fontsize=6);ax.grid(axis='x',alpha=.2)
    fig.suptitle('All fixed spatial sensitivity settings: effects and patient-bootstrap intervals; no added P tests',fontsize=12)
    fig.tight_layout(rect=(0,0,1,.95));fig.savefig(outdir/'all_spatial_sensitivities.png',dpi=160);plt.close(fig)
    (outdir/'render_receipt.json').write_text(json.dumps(dict(status='RENDERED_COMPLETE_SOURCE_AGGREGATE_FAMILY',
        source_SHA256=hashlib.sha256(source.read_bytes()).hexdigest(),main_slots=12,reference_slots=4,
        matched_primary_BH_and_CI_slots=len(robust),no_individual_or_learned_arrays=True),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--outdir',type=Path,required=True);a=p.parse_args();main(a.source,a.outdir)
