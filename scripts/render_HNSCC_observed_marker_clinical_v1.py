"""Render only already accepted source-level clinical aggregates, privately."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(exist_ok=False)
    p=args.task/'production'
    read=lambda n:json.loads((p/n).read_text())
    terminal=read('private_terminal.json')
    assert terminal['status']=='COMPLETE_CALIBRATION_ACCEPTED'
    rows=read('private_all12_marginal_clinical_results.json')
    pairs=read('private_all6_pair_clinical_conjunctions.json')
    sens=read('private_all48_descriptive_sensitivities.json')
    labels=[r['phase']+' | '+r['state'] for r in rows]
    y=np.arange(len(rows))[::-1]
    fig,ax=plt.subplots(figsize=(10.5,7.5))
    for yy,r in zip(y,rows):
        if r.get('CI95'):
            lo,hi=r['CI95'];ax.plot([lo,hi],[yy,yy],color='#6b7280',lw=1.5)
        if 'rho' in r:ax.plot(r['rho'],yy,'o',color='#31688e' if r['phase']=='pre' else '#b65c31')
        ax.text(1.05,yy,f"{r['n']:d}      {r.get('auxiliary_BH12',1):.3f}",va='center',fontsize=9)
    ax.axvline(0,color='black',ls='--',lw=.8)
    ax.axhline(5.5,color='#d1d5db',lw=.8)
    ax.set(yticks=y,yticklabels=labels,xlim=(-1,1),xlabel='Adjusted partial rank correlation with pathological regression',
        title='GSE286827 | fixed RNA states and pathological response')
    ax.text(1.05,12,'n      BH12',fontsize=9)
    fig.text(.02,.025,'95% patient bootstrap intervals; therapy, HPV and own-parent capture/abundance controls.\nExploratory reused source; separate phases share patients. No clinical prediction increment or direct communication tested.',fontsize=9)
    fig.subplots_adjust(left=.25,right=.84,bottom=.16,top=.9)
    fig.savefig(args.output/'clinical_state_forest.png',dpi=170);plt.close(fig)
    settings=['main','clinical_background_only','all_retained_RNA','parent_min10','parent_min30']
    lookup={(r['phase'],r['state'],r['setting']):r for r in sens}
    values=np.asarray([[r.get('rho',np.nan)]+[lookup[r['phase'],r['state'],s].get('rho') for s in settings[1:]] for r in rows],float)
    fig,ax=plt.subplots(figsize=(10,6.7))
    im=ax.imshow(values,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
    for i in range(len(rows)):
        for j in range(len(settings)):
            if np.isfinite(values[i,j]):ax.text(j,i,f'{values[i,j]:+.2f}',ha='center',va='center',fontsize=9)
    ax.set(yticks=np.arange(len(rows)),yticklabels=labels,xticks=np.arange(5),
        xticklabels=['Full adjustment','Therapy/HPV only','All retained RNA','Parent >=10','Parent >=30'],
        title='Fixed descriptive sensitivities | no additional P values')
    ax.tick_params(axis='x',labelrotation=25)
    fig.colorbar(im,ax=ax,shrink=.7,label='Partial rank correlation')
    fig.text(.02,.025,'Main analysis: both RNA parents >=20; retained-source RNA QC500/250/MT20.\nSensitivity estimates do not replace the fixed primary family or resolve wide intervals.',fontsize=9)
    fig.subplots_adjust(left=.25,right=.94,bottom=.22,top=.91)
    fig.savefig(args.output/'clinical_capture_sensitivity.png',dpi=170);plt.close(fig)
    report=['# GSE286827固定RNA状态—病理关联：完整探索结果','',
        '## 本轮结果','',
        '采用连续病理退缩率；正值表示状态比例较高与退缩较多同向。两阶段分开，不能当两次独立患者验证。完整临床桥、原始counts计数、观测规则、全部真实随机流R复算与两家族BH/BY接受；工作模型有限报警通过，不表示所有抽样假设均已证明。','',
        '|阶段|固定观测状态|患者数|调整后相关|95%病例区间|原始P|辅助BH12|',
        '|---|---|---:|---:|---|---:|---:|']
    for r in rows:
        ci=r.get('CI95');ci=f'[{ci[0]:+.3f}, {ci[1]:+.3f}]' if ci else '缺资格'
        report.append(f"|{r['phase']}|{r['state']}|{r['n']}|{r.get('rho',float('nan')):+.3f}|{ci}|{r.get('raw_P',1):.4f}|{r.get('auxiliary_BH12',1):.4f}|")
    report += ['','## 组合：两端是否都关联同一病理终点','',
        '此处max P是交并检验，不是互作项、协同效应、临床预测增量或配体受体检验。','',
        '|阶段|调节相关端|效应相关端|主要BH6|与外部共变方向相容（仅描述）|',
        '|---|---|---|---:|---|']
    for r in pairs:report.append(f"|{r['phase']}|{r['R_state']}|{r['E_state']}|{r['primary_BH6']:.4f}|{r['endpoint_clinical_directions_compatible_with_external_covariance']}|")
    report += ['','## 解释与下一步','',
        '本轮有足够病例计算固定模型，却没有获得多重校正后支持；所有边际区间跨零，不能表述为“无法检验”，也不能表述为已证明无生物作用。病例数约数十、区间宽，尚不足支持临床增量或高自由度交互预测。只调治疗/HPV与完整调整的符号差异也应保留，不能挑较顺眼模型。','',
        'R/E是分子阳性观测父群。R_TNFRSF9、R_LAG3、R_cycling分别在调节相关父群内读TNFRSF9、LAG3、MKI67且TOP2A/UBE2C任一；E_dysfunction为效应相关父群内PDCD1且TOX/HAVCR2任一，E_memory为TCF7/IL7R任一，E_cycling沿用同一增殖规则。这些名称不证明正式Treg功能、耗竭或抗肿瘤能力。','',
        'CM02继续作图谱重点线；独立GSE288199的来源内状态共变证据与本轮临床不支持、既有空间主检验不支持并列呈现，不拼接成已证实耐药生态位。GSE288199 Table S1缺失已记录，待用户补齐后才能在同一明确治疗阶段/终点接病例结局；当前不重试附件、不猜疗效、不改检验阈值。','',
        '下一优先是把已观察的状态共变与目前尚未验证的组织关系/临床增量分开：保留三个外部状态组合为待验证图谱线索，临床与空间不升格；补表后另冻同来源同阶段同终点问题。若仍不支持，转向更大、测量完整且终点一致的独立队列，而非在当前小来源轮换终点/阈值。','',
        '原CM、正式身份、权威Zarr和旧CM/MCP/作者计数临床家族保持。所有个体数据/随机流仅RAM，报告和图仅私有。']
    (args.output/'REPORT_zh.md').write_text('\n'.join(report)+'\n')


if __name__=='__main__':main()
