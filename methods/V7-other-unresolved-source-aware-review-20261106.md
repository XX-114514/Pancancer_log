# V7其余unresolved：父身份、QC与状态分开复查

## 本次真正执行的范围

结果权威入口为[四类结果与59群完整表](../inventories/analysis/2026-11-06/v7_other_unresolved_review_v1/README.md)。固定V7全4,676,787行及95,038基因轴，继承前轮QC-v2。选择旧Erythroid、Platelet、HSPC、ILC及两个相邻边界、Ambiguous、Unknown，共298,674个当前保留细胞；不是重新运行全V7所有谱系或所有doublet工具。

从8个入口提取全基因原始整数counts及213个身份/状态marker。缺失基因先按每来源原始var核实；模块记录组件是否可测、raw≥1/≥2的正向证据。阈值是项目诊断口径，**不是文献通用分类器或概率**。不把原父标签、一次全局cluster或CellTypist多数票当最终身份。旧父类以外的四类单一候选入口只用于召回，并非确认身份。

## 实际QC操作

同一细胞须同时有两个已注册的成熟、相互竞争谱系组件，且原doublet分数在同来源/同样本/同方法中经验分位≥0.95、有效评分n≥100，才进入本轮联合候选。分位数不是doublet概率。具体[配对与例外](../inventories/analysis/2026-11-06/v7_other_unresolved_review_v1/policies/QC_pair_policy_v1.json)全部公开。

不自动排除T+NK、T+ILC、B-ASC过渡、HSPC髓系/红系/淋巴priming、红系-巨核相邻程序，也不把内皮-纤维共享或上皮-间质程序直接当doublet。HSPC对内皮的QC组件比普通候选入口严格：SPINK2/HLF、CD34/PROM1及PTPRC分别有正向表达，避免仅凭CD34/GATA2/MEIS1删内皮。血液肿瘤语境只作独立holdout依据，不作恶性确证；regex含词边界以避免Gallbladder匹配ALL。

逐对审阅来源数、各程序占全RNA份额、原评分及表达组件后，实际应用913个操作性混合RNA排除与2个血液语境暂存，写出新QC-v3 mask、私有细胞ledger和counts子集。该操作不是人工doublet的真值验证；细胞、原分数、规则及旧注释均保留。全球累计2,916排除、238暂存、4,673,633保留。每行QC与身份是否可解释分开：保留细胞仍可能需要身份复核。

## 从清理后的counts重新建模

四个模型分别9,098、2,318、14,072、45,290细胞。保留全95,038 RNA特征做全库LogNormalize(scale.factor=10,000)，不假定已获得完整去污染表达层。大多数是UMI数据，混合技术中统一称原始counts，不将所有数值称UMI。

模型候选特征须在≥80%保留来源的var中存在；不全局屏蔽globin，不先回归增殖、应激或分化程序。临时模型assay以VST选择2,000 HVG，ScaleData、50 PCA，前30 PC建图。先计算未校正res0.4，再以gse_id作Harmony协变量重新建图res0.2/0.4/0.6/0.8；gse_id只是来源，不能当已核实物理平台。seed=20261106，每模型16线程，UMAP在新图上重算。

res0.4只是统一**诊断分区**，得到18/11/16/14群；不预设生物学亚型数，也不以分辨率制造亚型。RNA FindAllMarkers在所有保留来源均可测的基因中计算，only.pos=FALSE、min.pct=0.1、min.diff.pct=0.25、logfc.threshold=0.25、return.thresh=0.01。展示每群正向adj.p<0.05的top15，完整DE留私有运行。部分群没有满足条件的DE，这不证明无身份或无差异；已在独立审计列出。来源集中与其他分辨率的群数均公开；仍需补按来源及参数的身份稳定性验收。

## marker与注释判定

DE使用RNA assay的Wilcoxon检验，是细胞级探索，未把患者或来源当独立重复进行差异推断；共同可测基因过滤也不消除来源混杂。DE用于辅助身份/程序命名，不作为来源集中小群的新亚型或患者层效应显著性证明。

首轮233 marker发现其他组织父身份，追加肝、腺泡、纤毛、髓系/骨吸收及S期背景，形成260个不同基因的E2。**补充panel是观察后的诊断扩展，不是预设统计验证**。每群对每基因同时报告可测细胞数、检出数、测量覆盖、检测比例、全库CP10K均值。独立程序逐基因从原始counts复算所有记录。图顶明确对应的身份/候选程序/状态，左侧为工作标签和n；点大小为可测细胞检出率，颜色为每基因内部相对表达，灰×为来源覆盖<80%，不解释成阴性。

- 红系先核GYPA/SLC4A1/RHAG/GYPB/EPB42与AHSP/ALAS2/KLF1。HB不能单独证明完整红系细胞；cycling只作状态，不硬分早/晚红系。
- 血小板/巨核谱系核PF4/PPBP与ITGA2B/GP9/GP1BA/GP1BB/TUBB1。NRGN共享，不能单独定血小板。RNA缺少形态/细胞核依据时保留Platelet_Megakaryocytic，不硬定有核巨核或新亚型。
- HSPC样程序核primitive组件、PTPRC、成熟谱系及组织/疾病语境。原候选与非血液细胞可共享CD34/GATA2/MEIS1。HSPC_like与髓系priming不证明正常HSC/MPP/GMP/MEP；白血病来源还需配套恶性/基因型及参考数据判断。
- ILC需先解决与T/NK及非造血父类的边界。IL7R/KLRB1/GATA3/RORC/KIT均有共享情境。正向T apparatus支持转T；TCR未检出不等于Lin阴性。多组件ILC2/3程序加PTPRC、没有当前正向竞争组件，只能进入ILC_like_candidate，不能自动接受最终ILC2/3。

每群人工决策分别保存working_parent、身份复查状态、命名marker、state_context、理由。明确大类的群给簇级工作标签，但逐细胞出现竞争程序者继续身份复查；未定群只在有直接T或非HB红系组件时作正向路由。组件未触发不作细胞身份阴性。HSPC相关priming不是混合排除理由。状态列是簇级RNA背景，**并未给每个细胞确证功能状态**。最终细分标志统一FALSE。

实际写出四类及合并工作侧表，保留旧V7标签/父类；27,074个大类工作记录无需本轮竞争复查，43,704个仍有明确复查原因。不能将前者称为全V7发表注释通过数。未重新建模的226,981个背景入口只保留路由，不冒充注释完成。CellTypist/scATOMIC继承预测按群top5列为辅助证据，本轮没有重新训练/运行分类器，不作为模型特征或硬裁决。

## 与原文的对应及界限

侧表写出后另对166个ILC-like候选跟进正向组件：PTPRC≥2，类型2程序为GATA3及PTGDR2/IL1RL1，类型3为RORC及IL23R/CCR6/IL22（各组件≥2 counts）。得到32来源、67类型2/99类型3程序候选，无同时满足两者；候选侧表已读回。这是选择后的RNA诊断，不是无偏丰度估计或最终ILC2/3注释。补充表只报告正向counts人数，没有将不可测/零表达充作Lin阴性或计算来源可测检出率。

[Lodi作者代码](https://github.com/lambrechtslab/Pan-cancer_Lodietal)支持谱系内复核QC、重建、RNA DEG和人工命名的顺序；本轮具体marker阈值、召回、配对、分位与侧表均是项目实现，不能称原文逐字复现。对应资源、逐基因来源和此前核读范围见[panel方法](V7-panel-repair-literature-workflow-20261105.md)。

[Björklund原始研究](https://pubmed.ncbi.nlm.nih.gov/26878113/)使用已富集的人CD127+ ILC/NK材料；本轮阅读范围为摘要等可获取内容，不能把其分型能力直接迁移到未分选泛癌RNA。Nature全文本轮访问失败，未宣称完整核读。[Velten原始研究](https://doi.org/10.1038/ncb3493)强调连续谱系承诺，RNA与功能证据不能互相替代；本轮沿用既有原文核读与可获取摘要/方法范围，没有声称所有正常HSPC亚型均已验证。

## 重放与审计

实际代码位于[scripts](../scripts/v7_other_unresolved_review_20261106/)。它们读取私有原生对象并拒绝覆盖已有阶段目录。先按policies将文件放入新的逻辑RUN_ROOT，A/B/C顺序运行；B候选审阅并形成明确adjudication后才能C，不能无审阅照搬旧人数到新数据。C独立审计通过后按每个refit_config调用R；真实terminal0后运行E2及独立marker审计；人工cluster_decisions审阅后写侧表，再独立审计。R与Python代码中具体run名属于本次可追溯版本，新输入须复制成新版本配置。

```bash
python scripts/v7_other_unresolved_review_20261106/run_stage.py \
  --project-root "${PROJECT_ROOT}" --stage audit_clean_inputs.py
PROJECT_ROOT="${PROJECT_ROOT}" Rscript \
  scripts/v7_other_unresolved_review_20261106/rare_lineage_refit_v1.R \
  "${RUN_ROOT}/stage_c_QC_views/Erythroid/refit_config.json"
python scripts/v7_other_unresolved_review_20261106/run_stage.py \
  --project-root "${PROJECT_ROOT}" --stage audit_postfit_v2.py --lineage Erythroid
```

独立输入审计核完整稀疏矩阵、10x缓冲及细胞/基因轴；postfit审计核所有分辨率计数、来源集中度、HVG/DE资格和260 marker统计；侧表审计核70,778行、不重叠、旧标签、工作路由、聚合数、SHA和冻结V7未改。审计不等于生物学终验。首次红系独立审计10000×int32发生溢出，改为浮点复算后通过；生产CP10K一直使用浮点缩放。失败记录保留。
