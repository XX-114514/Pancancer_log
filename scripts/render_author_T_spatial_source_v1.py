"""Render source-only coverage and complete descriptive spatial summaries."""
from pathlib import Path
import json
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render(task,version='source_figures_v1'):
    assert version in ['source_figures_v1','source_figures_v2']
    task=Path(task);out=task/version;out.mkdir(exist_ok=False)
    m=json.loads((task/'private_complete_measurement_v1.json').read_text())
    d=json.loads((task/'private_complete_description_v1.json').read_text())
    c=json.loads((task/'private_frozen_description_v1.json').read_text())
    labels=c['all_CD8_labels'];regions=c['all_qualified_region_labels'];n=m['images']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12,max(4.5,.43*len(labels)+2)))
    y=np.arange(len(labels));parent=m['images_with_concordant_fine_T'];eligible=m['images_meeting_future_20each_fourgroup_gate']
    axes[0].barh(y-.18,[parent.get(k,0) for k in labels],height=.32,color='#94a3b8',label='Author type present')
    axes[0].barh(y+.18,[eligible.get(k,0) for k in labels],height=.32,color='#2563eb',label='All four groups >=20')
    axes[0].set_yticks(y,labels);axes[0].invert_yaxis();axes[0].set_xlim(0,n+1);axes[0].set_xlabel('Native files (not verified independent patients)');axes[0].legend(loc='lower left',bbox_to_anchor=(0,1.02),frameon=False,fontsize=8)
    primary=[s for s in d['all_spatial_slots'] if s['primary']]
    for i,label in enumerate(labels):
        slot=next(v for v in primary if v['CD8_author_label']==label)
        for delta,key,color,mark in [(-.12,'global_native_graph','#94a3b8','o'),(.12,'same_author_region_MH','#2563eb','s')]:
            value=slot[key]
            if value['median'] is not None:
                axes[1].plot(value['IQR'],[i+delta]*2,color=color,lw=2)
                axes[1].plot(value['median'],i+delta,marker=mark,color=color,ms=5)
            else:axes[1].text(0,i+delta,'NA',color=color,fontsize=7,ha='center',va='center')
        axes[1].text(1.02,i,f"{slot['global_native_graph']['estimable_images']}/{slot['same_author_region_MH']['estimable_images']}",transform=axes[1].get_yaxis_transform(),fontsize=8,va='center')
    axes[1].axvline(0,color='#64748b',lw=.7);axes[1].set_yticks(y,labels);axes[1].invert_yaxis()
    axes[1].set_xlabel('Median log odds; whiskers = image IQR, not CI')
    axes[1].set_title('k=10: gray overall; blue within author region\nRight: usable overall / region files',fontsize=10)
    fig.suptitle('GSE328481 author T types: coverage and tissue description',y=.98)
    fig.text(.5,.015,'CD4_Treg vs other CD4 anchors; CD8 subtype vs other CD8 receivers. Author labels, no clinical endpoint or P values.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.045,.98,.94));fig.savefig(out/'coverage_and_primary_description.png',dpi=170);fig.savefig(out/'coverage_and_primary_description.svg');plt.close(fig)
    rlabels=['CD4_Treg']+labels
    values=np.full((len(rlabels),len(regions)),np.nan);counts=np.zeros(values.shape,int)
    for item in d['all_regional_composition_slots']:
        i=rlabels.index(item['author_type']);j=regions.index(item['author_region_label']);v=item['relative_to_same_parent_elsewhere']
        if v['median'] is not None:values[i,j]=v['median']
        counts[i,j]=v['estimable_images']
    finite=np.abs(values[np.isfinite(values)]);limit=max(1,float(finite.max())) if len(finite) else 1
    fig,ax=plt.subplots(figsize=(max(7,1.25*len(regions)),max(4.5,.45*len(rlabels)+2)))
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#eeeeee')
    im=ax.imshow(values,cmap=cmap,vmin=-limit,vmax=limit,aspect='auto')
    for i in range(len(rlabels)):
        for j in range(len(regions)):
            v=values[i,j];ax.text(j,i,f'{v:.2f}\nn={counts[i,j]}' if np.isfinite(v) else 'NA',ha='center',va='center',fontsize=8,color='white' if np.isfinite(v) and abs(v)>.6*limit else 'black')
    ax.set_xticks(range(len(regions)),regions,rotation=25,ha='right');ax.set_yticks(range(len(rlabels)),rlabels)
    ax.set_title('Author region vs elsewhere, relative to the same T parent\nMedian image log odds; n = usable files, not patients',fontsize=11)
    fig.colorbar(im,ax=ax,label='Median log odds',shrink=.75)
    fig.text(.5,.015,'Unknown physical units / treatment / patient bridge. Region names and transferred T labels are not validated function.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.05,1,1));fig.savefig(out/'complete_regional_composition.png',dpi=170);fig.savefig(out/'complete_regional_composition.svg');plt.close(fig)
    lines=['# GSE328481作者T类型与组织关系：完整来源描述','','## 测量范围',
        f"完整{n}份原生作者文件，{m['rows']}条保留观察。有限/唯一坐标{m['coordinates_usable']}；全部精确重复坐标排除{m['duplicate_xy_rows_excluded']}。这不是已审阅独立患者数。",
        '原CM02、RNA功能代理、作者类型及RNA marker状态保持不同对象。此分析未读取表达矩阵，不能验证原CTLA4/GZMB或LAG3/增殖状态。',
        f"粗细T标签不一致观察{m['fine_T_not_coarse_T']}，粗T但不属已资格CD4/CD8细类{m['coarse_T_not_typed_CD4_CD8']}；后者不是不存在或非T的判断。",
        '', '### 完整作者T类型来源覆盖','','| 作者细类 | 一致父群且可用坐标观察 | 有该类型文件 |','|---|---|---|']
    for label,value in m['concordant_fine_T_counts'].items():lines.append(f"| {label} | {value} | {m['images_with_concordant_fine_T'].get(label,0)} |")
    lines+=['','### 原区域来源覆盖','','| 作者区域 | 原观察 | 有该标签文件 |','|---|---|---|']
    field=m['fields']['bin50 structural annotation']
    for label,value in field['qualified_counts'].items():lines.append(f"| {label} | {value} | {field['images_with_label'].get(label,0)} |")
    lines += [f"未通过披露词表的区域观察{field['unqualified_rows']}，缺失区域{field['missing_rows']}；前者只在源RAM按不同原类别分层，不合并或强行命名。",
        '', '## 全部空间描述槽',
        '数值是每文件log优势比的等文件中位数及IQR，无P/FDR/患者CI。正值表示相对于其他CD4锚点，作者CD4_Treg邻域的CD8中更偏向该作者细类；不等于总CD8更多或直接接触。主k10；k5/20保留敏感性。区域MH只保留双方锚点且至少10边的同原区域stratum，不加伪计数；粗比及区域组成使用0.5。',
        '', '| 作者CD8类型 | k | 全图 n; 中位数 [IQR] | 同原区域MH n; 中位数 [IQR] | 缺资格原因 |', '|---|---|---|---|---|']
    def fmt(v):return f"{v['estimable_images']}; {v['median']:.4f} [{v['IQR'][0]:.4f}, {v['IQR'][1]:.4f}]" if v['median'] is not None else '0; NA'
    for item in d['all_spatial_slots']:
        lines.append(f"| {item['CD8_author_label']} | {item['k']} | {fmt(item['global_native_graph'])} | {fmt(item['same_author_region_MH'])} | {json.dumps(item['missing_reasons'])}; region {json.dumps(item['region_missing_reasons'])} |")
    lines+=['','## 全部区域组成槽','','| 作者区域 | 作者类型 | 可用文件; 中位数 [IQR] |','|---|---|---|']
    for item in d['all_regional_composition_slots']:lines.append(f"| {item['author_region_label']} | {item['author_type']} | {fmt(item['relative_to_same_parent_elsewhere'])} |")
    lines+=['','## 解释边界','作者CD8_Treg字面保留，不能解释为常规CD8效应T。细类须与作者粗T一致；未知标签不输出、不当生物学缺失。所有缺资格位置保留。文件来自预选含未成熟TLS组织，不估计一般LUAD组织的无偏频率；未知阶段/临床桥，不能接到MPR结局。',
        '图中IQR及留一文件中位数范围是描述，不是患者置信区间。邻域来自全部原生坐标的kNN和固定5倍中位最近距离上限，含全部距离并列，不换算µm、不称物理直接接触。区域比较依作者字段，不把未知数值TLScluster或高T密度另命名TLS。',
        '此项不能确立CM02耐药/免疫抑制、因果通讯、同一状态跨来源复现或临床增量；后续须明确状态测量及可比较患者/阶段/终点。GSE288199 Table S1按用户安排后补。']
    (out/'REPORT_zh.md').write_text('\n'.join(lines)+'\n')
    files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    (out/'private_render_terminal_v1.json').write_text(json.dumps(dict(status='COMPLETE_SOURCE_RENDER_PENDING_VISUAL_QA',files_SHA256=files,
        renderer_code_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n')
    print('SOURCE_AGGREGATE_REPORT_AND_TWO_FIGURES_RENDERED')


if __name__=='__main__':
    import sys
    render(sys.argv[1],sys.argv[2] if len(sys.argv)>2 else 'source_figures_v1')
