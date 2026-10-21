"""Source aggregate capture and pilot depth, not new effect inference."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def render(task):
    task=Path(task);out=task/'source_figures_v1';out.mkdir(exist_ok=False)
    capture=json.loads((task/'private_complete_existing_RNA_capture_audit_v1.json').read_text())
    pilot=json.loads((task/'private_specific_gene_pilot_summary_v1.json').read_text())
    groups=sorted({r['author_group'] for r in capture['all_group_panel_capture']})
    panels=sorted({r['original_panel'] for r in capture['all_group_panel_capture']})
    display={'author_CD4_Treg':'Author CD4_Treg','other_author_CD4':'Other author CD4','author_CD8_CXCL13':'Author CD8_CXCL13',
             'other_author_CD8':'Other author CD8','author_gdT':'Author gdT','other_typed_CD4_CD8':'All typed CD4/CD8 comparator'}
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(16,6),sharey=True)
    for ax,metric,title in zip(axes,['cells_any_gene','cells_ge2_genes'],['At least one panel gene','At least two panel genes']):
        mat=np.array([[next(v['source_pooled_observation_fractions'][metric] for v in capture['all_group_panel_capture'] if v['author_group']==g and v['original_panel']==p)*100 for p in panels] for g in groups])
        im=ax.imshow(mat,cmap='Blues',vmin=0,vmax=max(1,mat.max()),aspect='auto')
        for i in range(len(groups)):
            for j in range(len(panels)):ax.text(j,i,f'{mat[i,j]:.2f}',ha='center',va='center',color='white' if mat[i,j]>.6*mat.max() else 'black',fontsize=8)
        ax.set_xticks(range(len(panels)),[v.replace('_',' ') for v in panels],rotation=35,ha='right')
        ax.set_yticks(range(len(groups)),[display[g] for g in groups]);ax.set_title(title)
        fig.colorbar(im,ax=ax,shrink=.7,label='Observed fraction (%)')
    fig.suptitle('GSE328481: all original files, fixed author groups and RNA panels',fontsize=13)
    fig.text(.5,.025,'Pooled observation capture, not identity accuracy or patient effects. All original panels retained. No new hypothesis tests.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.065,1,.94));fig.savefig(out/'full_source_fixed_panel_capture.png',dpi=170);fig.savefig(out/'full_source_fixed_panel_capture.svg');plt.close(fig)
    rows=pilot['all_types_and_markers'];y=np.arange(len(rows))
    fig,axes=plt.subplots(1,2,figsize=(12,max(6,.38*len(rows)+1.5)),sharey=True)
    for ax,key,threshold,title in zip(axes,['deposited_total_quartiles','detected_gene_quartiles'],[500,250],['Deposited RNA total (not verified UMIs)','Detected genes']):
        for i,row in enumerate(rows):
            q=row['quality'][key]
            if q is None:continue
            ax.plot([q[0],q[2]],[i,i],color='#64748b',lw=2);ax.plot(q[1],i,'o',color='#2563eb',ms=5)
        ax.axvline(threshold,color='#dc2626',ls='--',lw=.8)
        ax.set_yticks(y,[f"{r['author_type']} (n={r['quality']['n']})" for r in rows]);ax.set_xlabel(title);ax.grid(axis='x',alpha=.15)
    axes[0].invert_yaxis();fig.suptitle('Same minimum-size pilot: author T rows, median and IQR',fontsize=13)
    fig.text(.5,.02,'Dashed 500/250 lines are measurement references. They do not define cross-platform identity or viability.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.055,1,.95));fig.savefig(out/'pilot_author_T_depth.png',dpi=170);fig.savefig(out/'pilot_author_T_depth.svg');plt.close(fig)
    lines=['# CM02作者组织类型与RNA/临床桥：测量资格报告','','## 已完成的范围',
        '复用全部11份原空间文件的作者类型固定RNA panel捕获，并在同一个历史最小试点新核具体marker及深度。所有新原文件/个体表达和质量数组只RAM；仅来源级汇总保存。试点不是全部11份特定基因的资格。',
        '既有GSE243013临床家族含CD8T_Tex_CXCL13，但调节端为Treg内TNFRSF9/LAG3/MKI67阳性比例，不是空间作者Treg总量。该来源是新辅助后手术材料；同癌种LUAD分层未达原资格门槛，患者逐材料身份桥仍未接受。旧相应家族没有支持把该组织候选称为同状态临床组合。',
        '', '## 全来源固定panel捕获', '', '| 作者群 | panel | 原观察数 | 至少1基因 (%) | 至少2基因 (%) | 至少2基因各值>=2 (%) |','|---|---|---|---|---|---|']
    for r in capture['all_group_panel_capture']:
        v=r['source_pooled_observation_fractions'];lines.append(f"| {r['author_group']} | {r['original_panel']} | {r['observations']} | {100*v['cells_any_gene']:.3f} | {100*v['cells_ge2_genes']:.3f} | {100*v['cells_ge2_genes_each_value_ge2']:.3f} |")
    lines+=['','## 单一试点的完整作者类型深度','','| 作者类型 | n | RNA Q1/中位/Q3 | 检出基因 Q1/中位/Q3 | 同时达500/250参考 |','|---|---|---|---|---|']
    for row in rows:
        q=row['quality'];lines.append(f"| {row['author_type']} | {q['n']} | {q['deposited_total_quartiles']} | {q['detected_gene_quartiles']} | {q['rows_RNA500_genes250']} |")
    keys=['CD3D','CD3E','TRAC','FOXP3','IL2RA','CTLA4','CXCL13','GZMB','PRF1','LAG3','TNFRSF9','MKI67','TRDC','TRGC1','TRGC2']
    lines+=['','## 指定类型的具体marker检出记录数', '原始count-like RNA不是已确认原UMI；只报告基因检出，不由缺失认定生物阴性。所有既定marker完整记录留在来源级JSON；表中只展示预定关键检出字段。','','| 作者类型 | n | marker:检出数 |','|---|---|---|']
    for row in rows:
        if row['author_type'] in ['CD4_Treg','CD8_CXCL13','gdT']:
            v={r['marker']:r for r in row['all_marker_marginals']};values=', '.join(f"{k}:{v[k]['positive_rows']}" for k in keys)
            lines.append(f"| {row['author_type']} | {row['quality']['n']} | {values} |")
    lines+=['','## 当前解释',
        '全部原来源的作者T类pan-T捕获相对其他作者粗类有方向一致的提高，支持其包含T相关RNA背景，但不是逐细胞身份验证。更细调节标签和具体marker状态需要区别。试点特定基因及低RNA深度不支持直接确认FOXP3/CTLA4调节状态、CXCL13持续表达状态或γδ身份；也不能据此说作者标签必然错误或这些细胞不存在。',
        '空间Treg–CXCL13相关CD8的同原区域描述仍保留，当前限定为作者类型推断的组织假设。源选择、参考转移、细分身份、RNA捕获和临床阶段/患者单位均影响进一步解释，尚不形成直接通讯或临床增量链。未按这些结果降低marker或病例门槛。',
        '下一应优先在具备更深RNA/匹配蛋白及明确患者阶段的材料上核实状态。该低深度来源可用于保守组织类型描述，若做组织程序/局部聚合则须明确改换估计对象、另冻完整方法，不当同细胞状态的补证。原CM02保持图谱重点；原RNA代理临床线索、当前作者空间类型及已闭临床状态家族分别报告。',
        '原附件已是已接入的HNSCC表，与肺癌临床身份桥不同。GSE288199 Table S1仍依用户安排待补。受拒生存表/401表达对象/403代码与公共推送均未重试。']
    (out/'REPORT_zh.md').write_text('\n'.join(lines)+'\n')
    terminal=dict(status='COMPLETE_CAPTURE_AND_DEPTH_REPORT_RENDER_PENDING_VISUAL_QA',completed_UTC=datetime.now(timezone.utc).isoformat(),
        files_SHA256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in out.iterdir() if f.is_file()},
        code_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'private_render_terminal_v1.json').write_text(json.dumps(terminal,indent=2)+'\n');print(terminal['status'])


if __name__=='__main__':
    import sys
    render(sys.argv[1])
