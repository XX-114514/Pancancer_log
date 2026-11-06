# 2026-11-04 TME v1 整合检查

- 用户要求：检查最新状态和进一步计划。
- 已观察（job 30654587，bigmem，6.5 h）：TME v1 整合于 2026-11-04 06:40 UTC 完成。3,383,398个细胞，71个来源，1,211个单位，978名患者，101个 Leiden 簇。
- 已观察：patch job 30658776 把3项决定写入了 zarr。covarnet_v1_include=True 的细胞为3,015,706个；GSE189903 的71,616个细胞已改为冻存。
- 问题（已观察）：
  - Harmony 使用来源 θ2 加单位 θ1 后，来源混合反而变差。20万子样本、k=90：来源 iLISI 从4.26降到3.02，同来源近邻比例从0.35升到0.48。单位混合改善：iLISI 从13.8升到17.9。
  - 按来源分面的 UMAP 中，T/NK 云内有明显的来源条带。
  - 73个不少于2,000细胞的簇中，31个（652,715个细胞）有≥80%的细胞来自同一来源。其中一部分可能有生物学原因（解释，未验证），例如胶质瘤小胶质细胞（簇51，GSE182109）。但 T/NK 簇71（GSE234933 98%）和簇68（GSE207422 99%）应视为批次残留。
  - 已有的 cLISI 是循环指标：区室本身就是在 Harmony 图上定义的。
- 另见：uncertain 细胞集中在若干上皮/黑色素细胞簇（uncertain 占比0.6–0.9）。解释：可能是未判出的肿瘤细胞，支持 CoVarNet v1 排除 uncertain 的决定。
- 下一步：job 30658795 在30万子样本上比较 Harmony 各变体，生物学标签用与嵌入无关的单细胞 marker 标签。按结果决定重跑整合的方案（必要时改用 scVI），重跑后再做细注释。

## 去批次方法比较 v1（2026-11-04 12:00 UTC / 20:00 UTC+8）

- 用户要求：尝试 BBKNN 整合，并比较 scVI 等其他去批次算法。
- 设计：从 TME v1 固定取30万细胞（seed 2026110402），所有方法用同一子样本。
  - 输入：4,000个 TME HVG 的原始计数（供 scVI）、lognorm（供 ComBat、Scanorama），以及全数据 PCA。
  - 生物学标签：只由表达计算的单细胞 marker 标签，分粗、细两级。细级包括 CD8/CD4/Treg/NK、单核/巨噬/cDC1/cDC2/中性粒、myCAF/iCAF。
- 方法：
  - unintegrated；
  - 全数据 v1 Harmony；
  - Harmony（来源 θ2）；
  - Harmony（来源 θ2 + 单位 θ1）；
  - ComBat；
  - BBKNN（来源，每批3个邻居，trim 30）；
  - Harmony→BBKNN；
  - scVI（batch=来源）；
  - scVI（batch=单位）；
  - Scanorama。
- 指标：
  - 图指标（所有方法都算）：同来源/同单位近邻比例、图上的来源 iLISI、粗/细标签纯度、各标签的图连通性。
  - scib-metrics Benchmarker：只对嵌入方法计算，用10万个有置信粗标签的细胞。
- 环境变更：
  - 已有 `envs/conda/scvi` 环境（scvi-tools 1.4.1、scib-metrics 0.5.6、torch 2.5.1/CUDA 12.6）。
  - 用 pip --no-deps 增装了 scanorama 1.7.4、annoy 1.17.3、fbpca 1.0、geosketch 1.3、intervaltree 3.1.0、sortedcontainers 2.4.0。
  - BBKNN 1.6.0 原本就在 scanpy 环境中。
- 作业：prep 30658978 → cpu 30658979、scvi 30658980（A100）、scanorama 30658981 → eval 30658982。
- 此前的 Harmony 子样本基准（30658795）因标记基因缺失报错，已被本比较取代并取消。

## 来源、队列与样本组织来源登记 v1（2026-11-04 13:00 UTC / 21:00 UTC+8）

- 用户要求：确定数据集来源、队列信息和样本的组织来源。
- 方法：
  - 2026-11-04 重新抓取71个 GEO family SOFT，放在 `${PC_ROOT}/resources/geo_soft_20261104/`，附 SHA256SUMS。
  - 用 PubMed esummary 获取97个 PMID。GSE138709 在 GEO 中没有 PMID，按 GEO citation 用标题检索得到32505533（J Hepatol 2020）。
  - 用 `scripts/source_registry/build_source_registry_v1.py` 把单位桥接到 GSM：先用 GSM 号，再用去 GSE 前缀并规范化后的标题、补充文件名、带词边界的标题、特征值。
- 结果：
  - 1,258个单位中1,103个唯一桥接到 GSM，没有重复使用的 GSM。
  - 其余155个分布在14个来源（GSE179994、GSE182434、GSE189889、GSE200997、GSE207493、GSE256326、GSE297041 等），只有来源级描述。GSE256326 已按 material_detail 人工核对了部位。
- 发现：
  - GSE182434 是 CD19+/CD20+ B 细胞与 CD19−/CD20− 非 B 细胞的分选配对（GEO overall design），违反不用富集数据的规则。待用户决定。
  - GSE189357 中 TD5、TD7、TD8 为 AIS（原位腺癌），按 WHO 2021 属前驱病变，与不用癌前病变的规则冲突。待用户决定。
  - GSE191288 的淋巴结标记是"lymph node stage"引起的误报；GSE253973 的非肿瘤胸膜和 PDO 样本没有纳入。两者已核实无问题。
- 待说明：
  - 多化学并存的来源：GSE207422（10x 5′、10x 3′、BD Rhapsody）、GSE152048（10x 3′、GEXSCOPE）、GSE222446 和 GSE256326（5′、3′）。整合时应考虑以单位或化学为批次。
  - GSE182109 的 GEO contact_institute 字段无效。
  - 16个来源的 10x 版本在方案文本中没有写明。
- 产物：
  - 已提交：`inventories/source_registry/source_registry_v1.tsv`、`unit_flags_v1.tsv`。
  - 只留在服务器：`${PC_ROOT}/source_registry_v1/unit_registry_v1.tsv`（单位级）。

## 2026-11-05 范围排除与去批次比较结论（UTC / UTC+8 见各条）

- 用户决定（2026-11-05）：排除 GSE182434（B/非B 分选配对）和 GSE189357 的3个 AIS 样本（TD5、TD7、TD8）。
  - 记录在 `inventories/tme_v1/scope_exclusions_v1.tsv`，与 GSE179994 放在同一张表；integrate_tme_v1.py、report_dataset_detail_v1.py 和 patch_tme_scope_v1.py 都读这张表。
  - run_v1 zarr 已打补丁（job 30660289）：增加 scope_excluded 列，这些细胞的 covarnet_v1_include 设为 False。
  - TME 共3,347,659个细胞；CoVarNet v1 共2,980,723个细胞（此前为3,015,706）。
  - 数据集明细表已重新生成。run_v1 zarr 仍含被排除的35,739个细胞，所以 matches_sum=false；run_v2 不再包含这些细胞。
- 去批次比较 v1 的结果（30万细胞子样本；指标以图为基础，包括 scib-metrics）：
  - Harmony（来源 θ2）：同来源近邻0.154，来源 iLISI 5.64，细标签纯度0.701；scib 总分0.661，为最高。
  - BBKNN（来源）：混合最强（同来源0.014，iLISI 7.56），但细标签纯度最低（0.633），粗标签纯度0.77。解释：71个组成很不均衡的批次被强行配邻，导致过度混合。
  - Harmony→BBKNN：同样过度混合。
  - scVI（单位）：0.212 / 4.62 / 0.703，总分0.621。scVI（来源）：0.335 / 3.67 / 0.709。两者都跑满400个 epoch，没有触发早停。
  - ComBat：0.489 / 2.86，总分0.632。
  - v1 方案（来源 θ2 + 单位 θ1）：子样本0.602，全数据0.782，最差，与此前发现一致。
  - Scanorama：首次运行因 pip 版 annoy 在非 AVX512 节点上 Illegal instruction 失败，已加 AVX512 约束重跑（30660425），之后重新评估。
  - 注意：生物学指标用的是标记基因标签，偏向未整合的表达空间（未整合的纯度最高），只适合用于方法间比较。
- 决定：全数据 TME 重整合采用 Harmony 来源 θ2，并应用排除表。run_v2 job 30660430（bigmem）完成后由 30660431 写入决定补丁。

## 2026-11-05 TME run v2 完成（05:55 UTC / 13:55 UTC+8）

- 已观察：run v2（job 30660430，Harmony 来源 θ2）于04:50 UTC 完成。3,347,659个细胞，70个来源，1,201个单位，78个簇。post job 30660431 写入决定标记：covarnet_v1_include 共2,980,723个，排除的细胞为0，GSE189903 的71,616个细胞已改为冻存。
- 指标：来源 iLISI 从3.02升到12.82；同来源近邻比例从0.483降到0.067。≥2,000细胞的单一来源簇从31个降到5个。
- Scanorama 重跑完成（job 30660424，1,989 s）。同来源0.375，iLISI 3.38，细标签纯度0.744，scib 总分0.633，低于 Harmony（0.661）。
- 发现：簇13（290,546个细胞，uncertain 占76.5%）把各来源的 uncertain 上皮/肿瘤样细胞汇在一起，自动标签为 unresolved（353,834个细胞）。
- 下一步：粗区室人工核实（marker 点图），然后按区室细注释。是否沿用10月172节点分类作为骨架仍待用户决定。

## 2026-11-05 2026年10月注释（v20261010）可否作为细注释骨架（06:30 UTC / 14:30 UTC+8）

- 已观察：`Project_v3/annotation_rebuild/20261010_expression_annotation_downstream_v1/inputs/` 中有逐细胞标签。
  - meta.tsv（4,534,433行：sampleID、majorCluster、subCluster）。
  - exact_cell_roster.npz（global_position、172个 label、1,312个样本；SHA 在 run_contract.json 中）。
- 已观察：global_position 等于 core store 的 origin_native_row。在 GSE138709 和 GSE161529 上，roster 命中率为100%，单位与 roster 样本100%一致，可以精确连接。
- 覆盖：
  - TME run v2 中1,550,200/3,347,659个细胞（46.3%）有十月标签，来自36个 core 来源。
  - CoVarNet v1 中1,352,962/2,980,723个细胞（45.4%）有十月标签。
  - 3个 immune store 来源和32个补充来源没有标签。
- 方案（建议，未执行）：以十月172节点为骨架，在各区室的 Harmony 空间中把标签转移到未标注来源，并给出置信度；同时在各区室重新聚类，找出十月数据中没有的状态（新癌种、新来源）；两者交叉核对后定名。是否采用这一骨架由用户决定。

## 2026-11-05 注释面板合并与 Lodi 式注释逻辑 v2（10:15 UTC / 18:15 UTC+8）

- 用户要求：借鉴上传的 V7 Panel 工作簿（与仓库内 `resources/annotation/v7_panel_repair_v1/V7_panel_repair_v1.xlsx` 同为65,415字节），与本项目 v1 面板合并并给出文献依据；上传完整的面板和逻辑；读取 GitHub 上最新的注释策略修改；严格学习 Lodi et al. 的代码和论文。
- 已读取：
  - 远端 d61e0d9（V7 panel 修复）、44374b8（语义审计）、e6cf388（QC 和计算口径澄清）、92aa7ca（双细胞审计）、bacb550；
  - Lodi 仓库 @35f054b：`Cancer type analysis.R`、5个 `Cell subtype annotations/*.R`、`marker_genes.txt`、Minor_2.0 词汇；
  - Lodi 论文 PMC12629804 的 Methods。
- 产物：
  - `resources/annotation/tme_annotation_panel_v2/`：179个模块，通过 V7 validator，V7 测试15/15 通过。
  - `scripts/annotation_v2/build_panel_registry_v2.py`。
  - `methods/tme-annotation-v2-lodi-workflow-20261105.md`。
  - v1 方法文档标注为已被取代。
- 文献核对：
  - 原计划用于膀胱图谱的 PMID 31462534 实为肌肉论文，未采用；改用 Crossref 核实的 Yu 2019 JASN DOI。
  - NCBI esearch 和 Europe PMC 在本次会话中被限流或不可用，改用 Crossref DOI 核对标题。
  - Lodi 引用的"Chen et al, Cell 2021"未能确认，标为 unresolved。
- 关键设计变化：
  - 按 Lodi，先在各癌种内确定大类，并识别后移出癌种特有群体；
  - 跨癌种按谱系两轮分析（QC 轮，然后干净重建）；
  - 上皮和成纤维不做跨癌种细分。
  - 因此簇13的 v1 跨癌种子聚类（job 30667416）只作诊断。
- 与 Lodi 的差异（已在方法文档中写明）：
  - 批次变量用来源而不是 Technology（依据去批次比较 v1）；
  - 样本回归用按样本去均值实现（数学上等价）；
  - stress 列表不公开，用 V7 的 Stress_response 替代；
  - DoubletFinder 未安装，需补算统一的双细胞旁证分数。

## 2026-11-06 注释 v2 启动（00:40 UTC / 08:40 UTC+8）

- 用户要求：不要照搬 Lodi 参数，要理解每个参数的理由；说明簇13中6.8万个非恶性细胞如何处理；开始推进注释。
- 方法文档新增"参数的理由与本项目取值"表（fd69b89）。不照搬的项：
  - 不回归样本 ID：会把样本间的组成差异一并去掉，而组成差异正是 CoVarNet 要分析的信号；
  - 不设6000基因上限；
  - 双细胞预期率按每个样本的细胞数计算；
  - PC 数用 Kneedle 自动肘点；
  - 批次变量为来源；
  - HVG 额外排除 IG V 和 TCR V 基因；
  - 命名用的 marker 须在各来源中可重复上调。
- 采用的 Lodi 目的：回归状态程序（S/G2M、IFN、stress、hypoxia），使状态不形成伪亚型。
- 簇13诊断（job 30667416 完成）：68,374个非恶性细胞中，约1.5–2万个在身份清楚的器官上皮/间皮子簇；约4.8万个在由应激、糖酵解、增殖或缺氧程序驱动的跨癌种混合子簇（uncertain 占63–84%）。v1 的自动标签作废。处理：细胞回到各癌种的大类分析。
- A0：每个样本的 Scrublet 双细胞分数，job 30671847，按每个样本的细胞数设预期率，只作旁证。
- A2：`major_by_cancer_v1.py`。在 ICC 上试跑通过（4分钟，15个 PC，分辨率1下26簇）。修正：Lodi 大类列改为每组中检出最高的基因（每组含 CD4/CD8A 这类互斥基因）。27个癌种的正式运行为 job 30671957，在 A0 完成后启动。
- A3：`propose_major_labels_v1.py` 生成标签建议，规则透明，供人工复核。
