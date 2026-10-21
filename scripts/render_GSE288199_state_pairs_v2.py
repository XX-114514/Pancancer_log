"""Source-aggregate plots/report from an accepted full family; no individual I/O."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def number(x, digits=3):
    return 'NA' if x is None else format(x, '.'+str(digits)+'g')


def interval(x):
    return 'NA' if x is None else '['+', '.join(number(v) for v in x)+']'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    terminal = args.input/'private_terminal.json'
    t = json.loads(terminal.read_text())
    assert t['status'] == 'COMPLETE_FIXED18_SOURCE_STATE_ASSOCIATIONS_AND_ALL_REFERENCES_AND_SENSITIVITIES_R_AUDITED'
    for name, record in t['files'].items():
        assert hashlib.sha256((args.input/name).read_bytes()).hexdigest() == record['SHA256']
    read = lambda name: json.loads((args.input/name).read_text())
    data = read('private_complete18_state_pair_results.json')
    ref = read('private_complete18_matched_other_T_contrasts.json')
    sens = read('private_complete90_measurement_sensitivities.json')
    assert len(data) == len(ref) == 18 and len(sens) == 90
    key = lambda x: (x['phase'], x['R_state'], x['E_state'])
    assert len({key(x) for x in data}) == 18 and {key(x) for x in ref} == {key(x) for x in data}
    out = args.output
    assert not out.exists()
    out.mkdir(mode=0o700)
    labels = dict(R_TNFRSF9='Reg TNFRSF9+', R_LAG3='Reg LAG3+', R_cycling='Reg cycling RNA',
                  E_dysfunction='CD8 PDCD1/TOX/HAVCR2', E_cycling='CD8 cycling RNA', E_memory='CD8 TCF7/IL7R')
    fig, axes = plt.subplots(1,2,figsize=(14,8),sharey=True)
    for ax,phase in zip(axes,['pre','post']):
        selected = [x for x in data if x['phase']==phase]
        ax.axvline(0,color='#999999',lw=.9)
        for j,r in enumerate(selected):
            color = '#ac2f35' if r['BH18'] is not None and r['BH18']<.05 else '#2e6588'
            if r['CI'] is not None:
                ax.plot(r['CI'],[j,j],color=color,lw=1.5)
            if r['rho'] is not None:
                ax.plot(r['rho'],j,'o',color=color,ms=5)
            ax.text(1.03,j,'n='+str(r['eligible_source_units'])+' q='+number(r['BH18']),
                    transform=ax.get_yaxis_transform(),va='center',fontsize=8)
        ax.set_title(phase+' / complete prespecified family')
        ax.set_yticks(range(9),[labels[r['R_state']]+' | '+labels[r['E_state']] for r in selected],fontsize=9)
        ax.set_xlim(-1,1);ax.set_xlabel('Partial Spearman rho / 95% source-unit bootstrap CI')
        ax.grid(axis='x',alpha=.18)
    axes[0].invert_yaxis()
    fig.suptitle('GSE288199: outcome-blind RNA-pattern associations\nRed: BH18 < .05. Dots without bars: interval validity gate unmet.',fontsize=12)
    fig.text(.5,.018,'Coded source units are not reviewed independent patients; phases share cases. No treatment, clinical or direct-contact claim.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.045,.98,.92]);fig.savefig(out/'all18_state_pair_forest.png',dpi=150,bbox_inches='tight');plt.close(fig)
    settings=['alternate_parent','all_released_columns','low_capture','mid_capture','high_capture']
    lookup={(key(r),r['setting']):r for r in sens}
    values=np.full((18,6),np.nan);counts=np.zeros((18,6),int)
    for i,r in enumerate(data):
        values[i,0]=np.nan if r['rho'] is None else r['rho'];counts[i,0]=r['eligible_source_units']
        for j,s in enumerate(settings,1):
            v=lookup[(key(r),s)];values[i,j]=np.nan if v['rho'] is None else v['rho'];counts[i,j]=v['source_units']
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dddddd')
    fig,ax=plt.subplots(figsize=(12,13));im=ax.imshow(values,aspect='auto',cmap=cmap,vmin=-1,vmax=1)
    for i in range(18):
        for j in range(6):
            ax.text(j,i,('NA' if not np.isfinite(values[i,j]) else format(values[i,j],'.2f'))+' / n='+str(counts[i,j]),
                    ha='center',va='center',fontsize=8)
    ax.set_xticks(range(6),['Primary']+settings,rotation=25,ha='right',fontsize=9)
    ax.set_yticks(range(18),[r['phase']+' | '+labels[r['R_state']]+' | '+labels[r['E_state']] for r in data],fontsize=8)
    ax.set_title('All18 x primary + five measurement settings\nDifferent eligible sets: matched-set estimates are in the complete source report; no sensitivity P values.')
    fig.colorbar(im,ax=ax,fraction=.035,label='Partial Spearman rho');fig.tight_layout()
    fig.savefig(out/'all90_measurement_sensitivities.png',dpi=150,bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(14,8),sharey=True)
    for ax,phase in zip(axes,['pre','post']):
        selected=[r for r in ref if r['phase']==phase];ax.axvline(0,color='#999999',lw=.9)
        for j,r in enumerate(selected):
            if r['pointwise_CI'] is not None:ax.plot(r['pointwise_CI'],[j,j],color='#52658b')
            if r['difference'] is not None:ax.plot(r['difference'],j,'o',color='#52658b',ms=5)
        ax.set_title(phase+' / matched source units and12 controls')
        ax.set_yticks(range(9),[labels[r['R_state']]+' | '+labels[r['E_state']] for r in selected],fontsize=9)
        ax.set_xlim(-2,2);ax.set_xlabel('rho(Reg,CD8) - rho(otherT,CD8) / pointwise95% CI');ax.grid(axis='x',alpha=.18)
    axes[0].invert_yaxis();fig.suptitle('Complete matched other-T reference comparisons\nDescriptive intervals only; no P/FDR or formal specificity validation.',fontsize=12)
    fig.text(.5,.018,'Dots without bars: fewer than1900 of2000 paired bootstrap draws were estimable; intervals are withheld.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.045,1,.92]);fig.savefig(out/'all18_matched_other_T_contrasts.png',dpi=150,bbox_inches='tight');plt.close(fig)
    report=['# GSE288199：来源内T细胞RNA模式组合','',
        '本轮为新辅助HNSCC来源的结局盲探索分析。原CM02保持图谱优先；这里不是原CM02活动投影、正式细状态注释、直接通讯或临床增量。Table S1待用户补充。',
        '当前coded-case/material不是已审独立患者，pre/post共享来源病例；未知治疗臂、化学版本和其它混杂尚未消除。所有效应与区间条件于来源单位分析假设。',
        '', '## 完整主要家族','',
        '|阶段|调节相关RNA模式|CD8 RNA模式|可测来源单位|未调相关|调整相关|95%区间|BH18|BY18|状态|',
        '|---|---|---|---:|---:|---:|---|---:|---:|---|']
    for r in data:
        report.append('|'+ '|'.join([r['phase'],labels[r['R_state']],labels[r['E_state']],str(r['eligible_source_units']),
            number(r.get('unadjusted_spearman')),number(r['rho']),interval(r['CI']),number(r['BH18']),number(r['BY18']),r['status']])+'|')
    report+=['','## 如何读结果','',
        '调整相关控制双方各4项父群RNA质量/数量/组成，不能认为所有混杂已消除。CI要求2000次中至少1900次设计有效；缺失CI不补为零。所有18位置统一BH并列BY，未估计位置占原家族分母。',
        '红色或BH通过仅是当前家族的探索性支持；阶段间一致也不是独立队列复现。即使主相关通过，也须考察同材料其他T参照、捕获敏感性及未来精确临床桥，不能直接写成调节/抑制机制。',
        '','## 同材料其他T参照','',
        '|阶段|调节相关RNA模式|CD8 RNA模式|同材料数|R–E相关|其他T–E相关|差值|点区间|有效配对bootstrap|状态|',
        '|---|---|---|---:|---:|---:|---:|---|---:|---|']
    for r in ref:
        report.append('|'+ '|'.join([r['phase'],labels[r['R_state']],labels[r['E_state']],str(r['eligible_source_units']),
            number(r['R_rho']),number(r['other_T_rho']),number(r['difference']),interval(r['pointwise_CI']),str(r.get('valid_paired_bootstraps','NA')),r['status']])+'|')
    report+=['','其他T按同一RNA规则排除R/CD8/交叠池，不是经独立身份真值核实的非Treg。参照使用R/E/otherT共12项控制，与主模型的8控制不同；只能比较这两列同集同模型结果。有效配对bootstrap不足1900/2000时区间保持缺失，点估计不冒充精确区间。没有新P/FDR，点区间不自动验证特异性。',
        '','## 测量敏感性与范围','',
        '完整90项保留自身样本集及与原测量相同交集上的两套效应；图中为每个设置自身合格集，不能把不同n的差异当纯生物变化。所有结果见原production source JSON；未增加敏感性P/CI。',
        '','## 图','', '![完整18](all18_state_pair_forest.png)','',
        '![全部90敏感性](all90_measurement_sensitivities.png)','',
        '![同材料其他T参照](all18_matched_other_T_contrasts.png)','',
        '## 后续','',
        '先审全部方向、区间资格、共同RNA背景、参考细胞和测量敏感性。若有稳定阶段限定组合，待Table S1精确连接治疗/病例/阶段/病理后另立临床问题；若未支持，归档本家族，不改marker、阈值或端点救显著性。空间和临床证据仍需真正同状态/同终点资格。']
    (out/'REPORT_zh.md').write_text('\n'.join(report)+'\n')
    proof=dict(status='SOURCE_AGGREGATE_FIGURES_AND_FULL_REPORT_RENDERED_VISUAL_REVIEW_PENDING',
        accepted_terminal_SHA256=hashlib.sha256(terminal.read_bytes()).hexdigest(),
        files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()})
    (out/'private_render_receipt.json').write_text(json.dumps(proof,indent=2)+'\n')
    print(proof['status'])


if __name__=='__main__':
    main()
