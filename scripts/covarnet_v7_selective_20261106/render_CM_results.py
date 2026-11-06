"""Render aggregate author CM results without changing members, labels or thresholds."""
from pathlib import Path
import sys,json,textwrap,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

root=Path(sys.argv[1]);figdir=root/'figures';figdir.mkdir(exist_ok=True)
def wrap_label(value,width):
 # Wrap at annotation separators; never split a gene such as SOX2 or IL1B.
 return textwrap.fill(value.replace('_',' '),width,break_long_words=False,break_on_hyphens=False)
p=root/'primary_v2';o=p/'CM_results'
top=pd.read_csv(o/'factor_top15_candidates.tsv',sep='\t');members=pd.read_csv(o/'CM_members_after_author_network_filter.tsv',sep='\t');edges=pd.read_csv(o/'each_edge.tsv',sep='\t');summary=pd.read_csv(o/'CM_network_topology_and_permutation.tsv',sep='\t');catalogue=pd.read_csv(p/'inputs/node_catalogue.tsv',sep='\t');parents=dict(zip(catalogue.subCluster,catalogue.majorCluster));counts=dict(zip(catalogue.subCluster,catalogue.n_cells));singleton=dict(zip(catalogue.subCluster,catalogue.one_node_parent_in_input));types=sorted(catalogue.majorCluster.unique());colors={x:plt.get_cmap('tab20')(i) for i,x in enumerate(types)}
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})

fig,axes=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
for ax,variant,title in zip(axes,['primary_v2','scope_sensitivity_v1'],['Primary v2: 164 samples / 13 sources','Scope diagnostic v1: 210 samples / 14 sources']):
 z=pd.read_csv(root/variant/'CM_results/rank_survey_and_paper_selection.tsv',sep='\t');complete=json.loads((root/variant/'CM_results/complete.json').read_text());k=complete['selected_rank'];admit=z[z.paper_admissible]
 ax.plot(z['rank'],z.cophenetic,'o-',color='#526d82',lw=1.5,ms=4);ax.scatter(admit['rank'],admit.cophenetic,s=100,facecolors='none',edgecolors='#c07335',lw=1.6,label='Admissible CCC peaks');a=z[z['rank'].eq(k)].iloc[0];ax.scatter([k],[a.cophenetic],s=140,marker='*',color='#b33d3d',label=f'Selected K={k}',zorder=3);ax.set(xlabel='NMF rank (30 starts each)',ylabel='Cophenetic correlation',title=title,xticks=range(2,21,2));ax.grid(alpha=.18);ax.legend(fontsize=8,loc='lower right')
fig.suptitle('2026-11-06 V4 CoVarNet: paper CCC selection rule',fontsize=13);fig.savefig(figdir/'rank_selection.png',dpi=180);fig.savefig(figdir/'rank_selection.svg');plt.close(fig)
svg=figdir/'rank_selection.svg';svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')

diagnostics=[]
for q in summary.itertuples():
 cm=q.cm;t=top[top.cm.eq(cm)].sort_values('rank_in_factor');m=members[members.cm.eq(cm)];e=edges[edges.cm.eq(cm)];maxw=float(t.weight.max());negligible=m.loc[m.weight.lt(maxw*1e-6),'subCluster'].tolist();small=m.loc[m.subCluster.map(counts).lt(50),'subCluster'].tolist();single=m.loc[m.subCluster.map(singleton),'subCluster'].tolist();labels=list(m.subCluster)
 diagnostics.append({'cm':cm,'n_network_nodes':len(m),'n_edges':len(e),'n_components':q.n_components,'cross_parent_network':q.multicellular_network,'max_factor_W':maxw,'min_network_W':float(m.weight.min()) if len(m) else np.nan,'n_network_members_relative_W_lt_1e_6':len(negligible),'near_zero_W_members':'; '.join(negligible),'n_network_members_lt_50_input_cells':len(small),'rare_input_members':'; '.join(small),'single_node_parent_members':'; '.join(single),'same_member_set_as':'; '.join(a for a,g in members.groupby('cm') if a!=cm and set(g.subCluster)==set(m.subCluster)) if len(m) else '', 'diagnostic_only_no_primary_filter_changed':True})
 fig,(ax,net)=plt.subplots(1,2,figsize=(19,7.6),gridspec_kw={'width_ratios':[1.12,1]},layout='constrained');bars=np.arange(len(t));ax.barh(bars,t.weight,color=[colors[parents[x]] for x in t.subCluster]);ax.set_yticks(bars,[wrap_label(x,43)+f'  [n={counts[x]:,}]' for x in t.subCluster],fontsize=8);ax.invert_yaxis();ax.set(xlabel='Raw W loading (linear scale)',title='All paper top-15 candidates; full annotation + input cell count');ax.grid(axis='x',alpha=.15)
 if labels:
  n=len(labels);angles=np.arange(n)*2*np.pi/n+np.pi/2;coords={label:(float(np.cos(a)),float(np.sin(a))) for label,a in zip(labels,angles)}
  for z in e.itertuples():
   x,y=coords[z.subCluster1],coords[z.subCluster2];net.plot([x[0],y[0]],[x[1],y[1]],color=plt.get_cmap('viridis')((z.spe-.68)/.32),lw=1.5+2*z.correlation,zorder=1)
  for z in m.itertuples():
   x,y=coords[z.subCluster];net.scatter([x],[y],s=180+700*z.weight/maxw,color=colors[parents[z.subCluster]],edgecolors='white',linewidth=1,zorder=2);net.annotate(wrap_label(z.subCluster,27)+f'\nW={z.weight:.2g}',(x,y),xytext=(0,12),textcoords='offset points',ha='center',va='bottom',fontsize=8,bbox={'facecolor':'white','alpha':.88,'edgecolor':'none','pad':1.5},zorder=3)
 else:net.text(.5,.5,'No connected CM network\nafter the author edge filter',ha='center',va='center',transform=net.transAxes,fontsize=13)
 net.set(xlim=(-1.65,1.65),ylim=(-1.4,1.65));net.axis('off');net.set_title(f'Author-filtered network: {len(m)} nodes / {len(e)} edges\n{q.n_parent_types} parents / {q.n_components} components',fontsize=11);net.legend(handles=[Patch(facecolor=colors[a],label=a) for a in sorted({parents[x] for x in labels})],loc='lower center',fontsize=8,ncol=3,frameon=False);net.text(.5,-.025,'Node size: W / factor max; edge width: Pearson r\nEdge color: specificity (purple 0.68 → yellow 1)',ha='center',fontsize=8,transform=net.transAxes)
 note=[]
 if negligible:note.append('Near-zero W network member(s): '+str(len(negligible)))
 if q.n_components>1:note.append('Disconnected components retained')
 if single:note.append('Single-node denominator parent present')
 footer='; '.join(note) or 'Statistical network; source/seed/spatial validation pending'
 fig.suptitle(f'2026-11-06 V4 {cm} | primary v2 | connectivity P={q.permutation_P:.4g}\n{footer}',fontsize=13);fig.savefig(figdir/(cm+'.png'),dpi=150);plt.close(fig)

diag=pd.DataFrame(diagnostics);diag.to_csv(root/'CM_weight_and_overlap_diagnostics.tsv',sep='\t',index=False)
comparison=[]
for variant in ['primary_v2','scope_sensitivity_v1']:
 inp=json.loads((root/variant/'inputs/complete.json').read_text());c=json.loads((root/variant/'CM_results/complete.json').read_text());s=pd.read_csv(root/variant/'CM_results/CM_network_topology_and_permutation.tsv',sep='\t');comparison.append({'branch':variant,'n_cells':inp['n_non_epithelial_cells'],'n_samples':inp['n_admitted_samples'],'n_sources':inp['n_admitted_sources'],'n_nodes':inp['n_variable_input_nodes'],'selected_rank':c['selected_rank'],'nonempty_networks':int(s.n_network_edges.gt(0).sum()),'cross_parent_networks':int(s.multicellular_network.sum()),'empty_factors':int(s.n_network_edges.eq(0).sum()),'single_parent_nonempty_networks':int((s.n_network_edges.gt(0)&s.n_parent_types.eq(1)).sum())})
pd.DataFrame(comparison).to_csv(root/'source_scope_comparison.tsv',sep='\t',index=False)

lines=['# 2026-11-06 V4 CoVarNet原文配置CM结果','', '**主分析已完成并独立核验：560,871个非上皮细胞、164样本、13来源、45节点；K2–20每K30次（570次），原文CCC规则选K14。** 14个因子中11个有网络，10个含跨父类边，1个仅髓系内部；CM01/02/05为空图。编号仅属于本次V4运行。','', '[rank选择](figures/rank_selection.png) · [全部top15候选](primary_v2/CM_results/factor_top15_candidates.tsv) · [最终成员与W](primary_v2/CM_results/CM_members_after_author_network_filter.tsv) · [全部网络边](primary_v2/CM_results/each_edge.tsv) · [输入45标签及覆盖](primary_v2/inputs/node_catalogue.tsv) · [完整拓扑与置换](primary_v2/CM_results/CM_network_topology_and_permutation.tsv) · [独立最终审计](primary_v2/CM_results/FINAL_INDEPENDENT_AUDIT.json)','', '![CCC选择](figures/rank_selection.png)','', '网络只保留18条全局边。下表完整保留原文算法输出，不把空因子、单父类或近零载荷的边硬称生物学模块。**CM10、CM12、CM14存在近零W成员；CM06/12/14最终成员完全相同。** 特别是CM12两个成员W约7.9×10⁻⁹⁴及3.0×10⁻¹⁰³，虽然碰巧落入该因子top15并有全局相关边，也不足以支持该因子代表这一生物学模块。','', '| 因子/图 | 网络节点 | 边 | 父类 | 连通分量 | 置换P | 完整成员（沿用V4） | 诊断 |','|---|---:|---:|---:|---:|---:|---|---|']
for q,d in zip(summary.itertuples(),diagnostics):
 m=members[members.cm.eq(q.cm)];flags=[]
 if q.n_network_nodes==0:flags.append('空图')
 elif not q.multicellular_network:flags.append('仅同一父类')
 if q.n_components>1:flags.append('非单一连通模块')
 if d['n_network_members_relative_W_lt_1e_6']:flags.append('成员含近零W')
 if d['single_node_parent_members']:flags.append('父类出现/缺失信号')
 if d['same_member_set_as']:flags.append('成员同 '+d['same_member_set_as'])
 lines.append(f'| [{q.cm}](figures/{q.cm}.png) | {q.n_network_nodes} | {q.n_network_edges} | {q.n_parent_types} | {q.n_components} | {q.permutation_P:.4g} | '+('<br>'.join(m.subCluster) if len(m) else '—')+' | '+('；'.join(flags) or '待稳健性验证')+' |')
lines+=['','上述“近零”仅为透明诊断：网络成员W/本因子最大W<10⁻⁶；没有据此更改作者筛选或重新选择秩。[全部权重/重叠诊断](CM_weight_and_overlap_diagnostics.tsv)。不同因子网络可重叠，不应把10个跨父类网络当作10个独立生态系统。置换在固定全局图上抽取相同节点数，只评估此节点标签零分布的连接度；全部非空网络名义P及因子BH诊断虽<0.05，并未控制来源、癌种、稀有节点、参数或选择效应。','', 'CM11包含Glial单节点父类，强连接可能携带组织/疾病和父类出现信号；HSPC、Melanocytic、Platelet同属单节点父类。输入仍含极少细胞的节点（如cycling Erythroid为2、Platelet为7），按作者算法本轮未另加最小节点计数门槛；覆盖表公开，不能以minmax后变异等价于可靠生物学证据。','', '## 来源范围敏感性','', 'GSE179994在查看CM效应前因保留TME矩阵范围未确认暂存。旧v1不覆盖，单列[范围对照](source_scope_comparison.tsv)：v1含210样本、14来源、690,421细胞，选K19；v2为164样本、13来源、560,871细胞，选K14。两套网络编号不可直接一一对应。这说明发现对来源范围明显敏感；不是来源稳健性已通过。','', '## 输入和执行边界','', '[方法](../../../../methods/CoVarNet-V7-selective-author-CM-20261106.md) · [run](../../../../runs/20261106_v7_selective_CoVarNet_author_v1.md) · [执行代码与SHA](../../../../scripts/covarnet_v7_selective_20261106/) · [研究级平台补核](SOURCE_PLATFORM_EVIDENCE.md)。V4继承亚型不等于本轮新增亚型验证；12分母父类与作者91非上皮亚群体系不同，使用NMF0.28而非作者测试版本，复用既有QC。不能称原文数据/环境完全复现。','', '来源留出/种子与节点覆盖敏感性、疾病语境以及独立空间验证尚未执行。下一步应先核对近零载荷和重复成员、组织/来源依赖，避免直接给CM命名机制。正式作者输出保留；改进版阈值/分母若推进，必须另开版本作对照。样本/患者行、矩阵、完整RNA和作业日志保持私有。','', '生产2026-11-06 12:40:49 UTC/20:40:49 UTC+8结束，12:40:57 UTC/20:40:57 UTC+8全部独立核验通过。producer初始complete仍写pending；以较晚的独立R/FINAL审计为终验依据。']
(root/'README.md').write_text('\n'.join(lines)+'\n')
manifest=[{'path':str(f.relative_to(root)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted(root.rglob('*')) if f.is_file() and f.name!='local_render_manifest.json'];(root/'local_render_manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps({'figures':len(list(figdir.glob('*.png'))),'factors':len(summary),'near_zero_weight_factors':diag.loc[diag.n_network_members_relative_W_lt_1e_6.gt(0),'cm'].tolist(),'audit':'aggregate input only; no author members or edges modified'}))
