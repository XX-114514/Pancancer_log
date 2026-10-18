#!/usr/bin/env python3
"""Private aggregate figures/report, conditional on actual independent numeric acceptance."""
import argparse,json,html,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

AXES=['pCM02_regulatory_cytotoxic','pCM04_memory_GZMK_effector']
LABELS=['CTLA4/GZMB-associated proxy','GZMK/IL7R/CTSW-associated proxy']
PALETTE=['#0072B2','#D55E00']
def table(frame):
 def fmt(x):return format(x,'.5g') if isinstance(x,(float,np.floating)) else str(x).replace('|','/')
 return '\n'.join(['| '+' | '.join(map(str,frame.columns))+' |','|'+'|'.join(['---']*len(frame.columns))+'|']+['| '+' | '.join(fmt(x) for x in row)+' |' for row in frame.itertuples(index=False,name=None)])
def main(a):
 audit=json.loads((a.audit/'independent_audit_public.json').read_text());assert audit['n_slots']==304 and audit['full_family_BH_BY_pass']
 assert audit['producer_terminal_sha256']==hashlib.sha256((a.run/'terminal_public.json').read_bytes()).hexdigest()
 out=a.output;assert not out.exists();out.mkdir(mode=0o700);figdir=out/'figures';figdir.mkdir()
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':6,'axes.titlesize':7,'axes.labelsize':6,'xtick.labelsize':5,'ytick.labelsize':5,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
 results=pd.read_csv(a.run/'all_frozen_model_slots_private.tsv',sep='\t');entries=[]
 if a.external:
  exaudit=json.loads((a.external_audit/'independent_audit_public.json').read_text());assert exaudit['n_slots']==8 and exaudit['full_family_BH_BY_pass']
  results=pd.concat([results,pd.read_csv(a.external/'all_frozen_model_slots_private.tsv',sep='\t')],ignore_index=True)
 def save(fig,key,caption):
  fig.canvas.draw();width,height=fig.get_size_inches()*25.4;assert abs(width-183)<1e-4 and height<=170.01
  for ext in ['pdf','svg','png']:fig.savefig(figdir/(key+'.'+ext),dpi=400)
  plt.close(fig);entries.append(dict(id=key,caption=caption))
 def forest(ax,sub,labels,title):
  ax.axvline(1,color='#888888',lw=.7);ax.set_xscale('log');ax.set_xlim(.15,6.5)
  for i,(label,row) in enumerate(zip(labels,sub)):
   if row is not None and row.get('status')=='exploratory_model_numeric_audit_pass':
    effect,low,high=float(row['effect']),float(row['CI_low']),float(row['CI_high']);colour=PALETTE[AXES.index(row['axis'])]
    ax.plot([max(.15,low),min(6.5,high)],[i,i],color=colour,lw=.8)
    ax.plot(np.clip(effect,.15,6.5),i,'D',ms=3,color=colour,markerfacecolor=colour if row['BH_q']<.05 else 'white')
   else:ax.text(1,i,'NA',va='center',ha='center',color='#999999',fontsize=5)
  ax.set_yticks(range(len(labels)));ax.set_yticklabels(labels);ax.invert_yaxis();ax.set_title(title);ax.set_xlabel('Effect per SD (95% Wald CI)')
 cancers=json.loads((a.run/'frozen_plan_before_effects.json').read_text())['TCGA_cancers']
 for endpoint in ['PFI','OS']:
  fig,axes=plt.subplots(1,4,figsize=(183/25.4,168/25.4));fig.subplots_adjust(left=.08,right=.98,bottom=.14,top=.88,wspace=.8)
  for j,(axis,model) in enumerate((x,y) for x in AXES for y in ['clinical_basic','stage_sensitivity']):
   sub=results[results.family.eq('TCGA_survival')&results.axis.eq(axis)&results.endpoint.eq(endpoint)&results.model.eq(model)].set_index('cohort')
   rows=[sub.loc[c].to_dict() if c in sub.index else None for c in cancers]
   forest(axes[j],rows,cancers if j==0 else ['']*len(cancers),('Proxy02' if axis==AXES[0] else 'Proxy04')+'\n'+('Age/sex/T background' if model=='clinical_basic' else '+ overall stage'))
  fig.suptitle('TCGA '+endpoint+': tumour-type-specific exploratory prognosis',fontsize=8)
  fig.text(.08,.035,'Complete BH/BY family: 264 slots; solid diamonds BH q<0.05. NA retains failed gates.\nPrimary solid tumours only; not an ICI treatment-response analysis. Intervals beyond limits are truncated.',fontsize=5.5)
  save(fig,'TCGA_'+endpoint+'_all33_cancers','Cox HR per SD, source-specific standardized scores, fixed clinical adjustment. Complete264-slot BH/BY; unavailable retained. CI is model-based Wald. OS limitations follow TCGA-CDR. No cancer pooling or optimized thresholds.')
 fig,axes=plt.subplots(2,2,figsize=(183/25.4,165/25.4));fig.subplots_adjust(left=.21,right=.96,bottom=.13,top=.91,hspace=.65,wspace=.65)
 for ax,ep in zip(axes.flat,['response','OS','PFS','response_increment']):
  sub=results[results.endpoint.eq(ep)&~results.family.eq('TCGA_survival')]
  rows=[];labels=[]
  for c in ['IMvigor210','Riaz2017','Hugo2016','Jung2019','Rose2021']+(['SU2C_MARK'] if a.external else []):
   for axis in AXES:
    hit=sub[sub.cohort.eq(c)&sub.axis.eq(axis)];rows.append(hit.iloc[0].to_dict() if len(hit) else None);labels.append(c+' / '+('02' if axis==AXES[0] else '04'))
  forest(ax,rows,labels,{'response':'Objective response (OR)','OS':'Overall survival (HR)','PFS':'Progression-free survival (HR)','response_increment':'Balance beyond both arms (OR)'}[ep])
 fig.text(.15,.035,'Blue: proxy02; orange: proxy04. Solid diamonds: own complete-family BH q<0.05.\nAll endpoints/gates kept. SU2C has separate frozen external families. Observational association, not treatment-specific prediction.',fontsize=5.5)
 save(fig,'ICI_all_cohorts_complete_endpoints','Patient-level fixed proxy associations, pan-T-adjusted. Objective response=CR/PR vs SD/PD; unknown remains missing. Distinct complete response, survival and increment families; SU2C external families separate. Archived/unknown sampling times do not support baseline prediction. No SD reclassification to improve effects.')
 native=pd.read_csv(a.run/'native_T_panel_expression_private.tsv',sep='\t');genes=list(dict.fromkeys(native.gene));nodes=list(dict.fromkeys(native.node))
 matrix=native.pivot_table(index='node',columns='gene',values='mean_CP10K',aggfunc='first').reindex(index=nodes,columns=genes);val=np.log1p(matrix.to_numpy());val=val/np.nanmax(val,axis=0)
 fig,ax=plt.subplots(figsize=(183/25.4,125/25.4));fig.subplots_adjust(left=.19,right=.88,bottom=.23,top=.91)
 im=ax.imshow(val,cmap='cividis',aspect='auto',vmin=0,vmax=1);ax.set_xticks(range(len(genes)));ax.set_xticklabels(genes,rotation=90);ax.set_yticks(range(len(nodes)));ax.set_yticklabels(nodes);ax.set_title('Full-native T-node mean expression: fixed panels');fig.colorbar(im,ax=ax,fraction=.025,pad=.02).set_label('Gene-relative log1p(mean CP10K)')
 save(fig,'Native_T_identity_and_panel_evidence','All saved native270-node sufficient statistics;16 T nodes shown, no cell sampling. Colour is log1p mean CP10K normalized per gene, not DEG significance or same-cell co-expression. Marker/state proposals remain unpromoted.')
 summary=[]
 for family,part in results.groupby('family',sort=False):
  good=part.status.eq('exploratory_model_numeric_audit_pass');summary.append(dict(family=family,n_planned=len(part),n_estimable=int(good.sum()),BH_significant=int((good&part.BH_q.lt(.05)).sum()),BY_significant=int((good&part.BY_q.lt(.05)).sum())))
 report=['# 优先免疫CM：TCGA与免疫治疗bulk探索','',
 '对象：2026-10-14免疫图谱CM（v1，T_whole_TME_primary_tumor）CM02/CM04。以下是事前固定多基因功能代理，不是原CM的W/S/H投影或特定细群复现。原HNSCC30家族没有重开。','',
 '## 完整范围与结果','',table(pd.DataFrame(summary)),'',
 '全部效应见all_model_slots_private.tsv；NA为按原规则不可估计，不是阴性效应。P值来自探索性LRT，区间为Wald；未做临床预测器训练或切点优化。TCGA只说明癌种内一般预后，未证明ICI特异性。','',
 '## 所有经BH的关联（含阴性方向与敏感性）','']
 significant=results[results.status.eq('exploratory_model_numeric_audit_pass')&results.BH_q.lt(.05)]
 columns=['family','cohort','axis','endpoint','model','n','events','effect','CI_low','CI_high','p','BH_q','BY_q']
 report.append(table(significant[columns]) if len(significant) else '本轮没有通过各自完整BH家族的关联。')
 report+=['','## 不可估计原因','',table(results[results.status.eq('unavailable')].groupby(['family','reason']).size().rename('n_slots').reset_index()),'',
 '## 解释边界与下一步','',
 '- TCGA PFI与OS不能互换，癌种间效应方向可能不同。基础模型与分期完整病例模型均保留，分期缺失并未填零。',
 '- bulk代理的两臂共同升高仍可能反映T细胞总量、共同炎症或细胞混合。共同出现不是两个特定亚群直接互作。',
 '- 增量检验考察固定两臂水平后的共富集/平衡；不是配体受体机制。未测纯度、治疗线数等仍可能混杂。',
 '- 独立数值审核证明计算一致，不等于临床验证；模型比例风险及复杂混杂仍需进一步诊断，不将初步HR称可用标志物。',
 '- 新生数据先核患者/时点/材料/事件/测量和跨研究重叠，不因本轮某方向较好就换签名或删家族。',
 '- 来源中病例可重复的Riaz只保留pre；Hugo的on样本排除。未知事件保持NA，Rose未对接临床的表达样本排除，材料/时点限制保留。','']
 (out/'RESEARCH_REPORT_CN.md').write_text('\n'.join(report));results.to_csv(out/'all_model_slots_private.tsv',sep='\t',index=False)
 page=['<!doctype html><meta charset="utf-8"><title>Priority CM exploration</title><style>body{font:16px system-ui;max-width:1100px;margin:30px auto}img{width:100%}p{line-height:1.6}</style><h1>Priority immune CM: exploratory functional proxies</h1><p>TCGA prognosis and ICI association are separate; not exact cellular module recovery or validated prediction.</p>']
 for e in entries:page+=['<h2>'+html.escape(e['id'])+'</h2><p>'+html.escape(e['caption'])+'</p><p><a href="figures/'+e['id']+'.pdf">PDF</a> · <a href="figures/'+e['id']+'.svg">SVG</a></p><img src="figures/'+e['id']+'.png">']
 (out/'index.html').write_text('\n'.join(page))
 manifest=dict(status='AGGREGATE_RESEARCH_FIGURES_RENDERED_VISUAL_REVIEW_PENDING',observed_utc=datetime.now(timezone.utc).isoformat(),n_figures=len(entries),figures=entries,patient_identifiers_rendered=False,summary=summary)
 (out/'render_manifest_private.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(dict(status=manifest['status'],n_figures=len(entries))))

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for name in ['run','audit','output']:p.add_argument('--'+name,type=Path,required=True)
 p.add_argument('--external',type=Path);p.add_argument('--external-audit',type=Path);main(p.parse_args())
