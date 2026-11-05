from pathlib import Path
import csv,json
r=Path(__file__).resolve().parents[1]
def rows(n):return list(csv.DictReader((r/'tables'/n).open(),delimiter='\t'))
def mdtable(rs,keys):return '| '+' | '.join(keys)+' |\n|'+ '|'.join(['---']*len(keys))+'|\n'+ '\n'.join('| '+' | '.join(str(x[k]) for k in keys)+' |' for x in rs)+'\n'
text='''# V5与V7注释差异及V7通讯重算方案

取证日期2026-11-03。已完成发布层面的全量大类转移、恶性主判定转移、细分标签计数与恶性证据层计数对比。**细胞级细分转移、每个历史CM成员的细胞流向及V7下游重算仍需原始对象，当前没有真实V7 CM或通讯新结果。**

## 比较对象与规则

V5指0825冻结的`final_annotation_v6_identity_v5_malignancy_v3`；V7指0905冻结的`pancancer_annotation_v7_cnv_rerun_20260905`。两版均4676787细胞、1322原始样本、54 GSE、43癌种，V7保留旧字段，未删细胞。0906旧通讯结果沿用V5。细胞宇宙相同不表示进入严格分析的细胞相同。

| 方面 | V5 | V7 |
|---|---|---|
| 大类字段 | major_celltype_final | major_celltype_v7 |
| 细分身份 | identity_annotation_final，需与最终对象字段核对 | final_annotation_v7 |
| 显示/旧CM词表 | final_annotation_with_malignancy_v6，87节点 | final_annotation_with_malignancy_v7，不能直接替换87词表 |
| 身份来源 | evidence v5 | uncertainty-resolution conservative V2 |
| 来源标签 | 旧自动标签的独立性在后续审计中被纠正 | previous_automatic不再冒充独立作者标签 |
| 负marker | 旧硬否决规则可被单个负marker触发 | 单负marker降为provisional；至少2个一致负marker才构成竞争程序 |
| 低margin | 多落入Unknown/Ambiguous | 增加T/NK/ILC、髓系等层级边界 |
| 大类冲突与doublet | 旧QC字段保留 | 跨compartment冲突且doublet高分位可排除严格大类分析，细胞仍保留 |
| subtype unresolved | 允许未解决 | 扩展panel救援；de novo描述后缀与成熟亚型分开 |
| 计数与基因 | 已修复真实文库分母、缺测掩码、部分别名 | 继承并核实，不应当作V7独有创新 |
| 恶性 | v3逐样本CNV证据 | 按V7谱系重建candidate/reference后重跑inferCNVpy及Copykat_python |
| 无恶性检出 | non_malignant | cnv_no_malignancy_detected，明确不等于病理良性 |
| 下游 | K9/LIANA已完成 | 旧分支未自动重算 |

## 大类完整对比

大类标签改变770833细胞，占16.4821%；原标签保持3905954。这个数是同名标签比较，不是错误率。Unknown与Ambiguous合计由846526（18.1006%）降至129843（2.7763%）；减少的细胞有一部分进入边界或provisional，并非全部获得严格身份。

'''
text+=mdtable(rows('major_v5_to_v7_marginals.tsv'),['label','V5','V7','delta','retained','inflow','outflow'])
text+='''
每个旧大类到每个新大类的完整非零转移及旧类别内百分比见`tables/major_v5_to_v7_transition_long.tsv`；原始完整矩阵也保留。表中同名细胞保留不等于其细分身份或CNV不变。

V7严格大类4057641（86.76%），非严格619146（13.24%）；严格亚型2455424（52.50%），非严格2221363（47.50%）。已知大类原1073743个subtype unresolved中，91805严格救援、609412获得描述性基因后缀、372526仍无稳健marker。后两类不能自动按成熟亚型参加主检验。

## 恶性主判定完整对比

同名主判定改变1053571（22.5277%），其中包含non_malignant改名造成的语义变化，不能称为105万细胞被纠错。若仅将旧non_malignant→新cnv_no_malignancy_detected的61525细胞视为语义保持，其余变化992046（21.2121%）。

'''
text+=mdtable(rows('malignancy_v3_to_v7_marginals.tsv'),['label','V5','V7','delta','retained','inflow','outflow'])
text+='\n完整细胞流向：\n\n'+mdtable(rows('malignancy_v3_to_v7_transition_long.tsv'),['v5','v7','n_cells','percent_of_v5'])
text+='''
旧malignant 679957中390245仍为malignant，142415改为not_candidate、138410为uncertain、8484为无检出、403为不可评估。V7 malignant共515031，其中124786来自其他旧类别。这既受身份与candidate池影响，也受重跑CNV影响，现有聚合表不能分离两者各自贡献。

V7角色：candidate1156120、reference2465967、其余1054700。inferCNV成功938样本，Copykat_python成功668；跳过不可当阴性。固定Copykat_python版本1.0.0、commit c2390e8a46e21fb42aa06f05839900c0f12564c5。

## 恶性证据层完整计数对比

'''
text+=mdtable(rows('malignancy_tier_count_comparison.tsv'),['label','V5','V7','delta'])
text+='''
证据tier与cluster主判定必须分开。V7交集阳性238418中30080在cluster主判定仍为uncertain；所以严格交集并非cluster恶性集合的子集。nullable的is_malignant_primary_v7对uncertain/not_evaluable/not_candidate保持NA。

## 全部细分标签对比

以下为两版细分身份标签并集的计数对比。**相同标签数量保持不能证明逐细胞身份保持；差值也不能恢复流向。** 当前发布未提供细胞级细分V5×V7交叉表，必须读取冻结V7完整checkpoint中的历史字段后补齐，不能依据大类转移推造。

'''
text+=mdtable(rows('fine_identity_count_comparison.tsv'),['label','V5','V7','delta'])
text+='''
## V7重算与四层通讯分析

1. **V7细胞级差异审计。** 核对完整checkpoint SHA、唯一细胞索引、样本与计数逐行一致。导出样本/GSE分层的大类、细分、malignancy与strict eligibility转移；历史CM01/07/08每个成员分别计算保留/转出/流入，记录恶性身份变化。
2. **重算CM。** 首轮K=9、30次初始化及历史组成归一化用于受控比较；另做rank survey与来源留出。V7全显示标签复现旧输入语义，严格身份分支研究可解释状态；两分支的细胞/节点/样本变化单列。依据细胞映射、成员载荷与样本活动共同对齐V5模块，新模块用V7 run编号，不预设新CM07等于旧CM07。
3. **重算LIANA。** 在每个原始样本的全部合格状态对上运行rank_aggregate（每群≥20，expr_prop 0.1，seed1337），然后添加V7 CM成员关系。保留CM外丰度匹配细胞对；分母为双方均满足测量/细胞数条件的患者或样本。0.01/0.05/0.10为rank阈值，不能当p值。
4. **成员空间邻近。** 先通过V7状态识别、原始panel和恶性身份资格；在患者/区室内固定半径图与密度匹配零模型比较。细胞比例和库大小作为混杂控制，同患者多ROI不能充作独立患者。
5. **blisa热点。** 50μm主bin、25/100μm敏感性，999次置换；按LR类别预设接触/扩散模式。保留all_pval及High–High方向，对每切片LR×bin族BH校正。热点数之外报告合格bin分母、效应和病例一致性。[官方API](https://chenlaboratory.r-universe.dev/blisa/doc/manual.html)
6. **发送/接收贡献。** blisa的counts_by_group/CCI_scores追踪状态对贡献，结合逐细胞邻接图复核热点内是否确有相邻发送者与接收者。CCI分数本身不是独立p值；用matched non-CM pairs与同状态不同局部暴露作比较。
7. **接收端响应。** LIANA+衔接LR与胞内通路/TF；以独立靶基因集计算同一接收状态内响应，剔除身份定义和LR直接重叠基因。患者×接收状态汇总，校正检测深度、发送/接收丰度、组织区室、来源；患者/来源留出。面板靶基因不足记不可评估。[原论文](https://www.nature.com/articles/s41556-024-01469-w)

## 已执行状态及输入缺口

当前真实结果是注释聚合比较，表的两版总数均核对为4676787。重算配置见`config/analysis_plan.json`，尚未产生V7 CM、LIANA、blisa或响应效应。

- 本机冻结目录已有V5 W/H、本体、总结；V7细胞级checkpoint与旧counts仅有登记，本机迁移目录未确认实物。
- 旧服务器SSH本次返回`No user exists for uid 6880`，在连接前即失败，不能据此声称远端数据缺失。
- scanpy环境检测到LIANA模块，但导入加载尚待完成；decoupler缺失。
- blisa/fastLISA/SpatialExperiment仍缺，下载r-universe仍为DNS解析失败。安装与冒烟入口位于上一任务目录。
- 32个已处理空间counts对象未同时覆盖COL1A1–DDR1或LGALS1–PTPRC；需原始panel/过滤审计及合适空间输入，再进入两条CAF主轴。缺测不能作为通讯阴性。

任务仅在可写ondemand目录追加；历史冻结对象和记录库未修改。完整细胞级差异与真实重算未完成，所缺输入及环境状态应随恢复继续更新。
'''
(r/'V5_V7_COMPARISON_CN.md').write_text(text)
