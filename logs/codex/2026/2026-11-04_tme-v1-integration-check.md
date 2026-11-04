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
