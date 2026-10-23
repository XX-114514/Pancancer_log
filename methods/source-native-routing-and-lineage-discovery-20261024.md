# 来源内粗身份与谱系细分：执行合同

承接[完整方案](full-context-atlas-hierarchical-annotation-plan-20261023.md)和[粒度比较](full-context-granularity-comparison-plan-20261024.md)。本轮产生试验范围的分配建议和发现群，最终身份、跨来源细身份协调及全量映射另行验收。

## 粗身份的证据和待审池

- 使用已有来源内A图的群体背景、原生多基因支持及已冻结参考预测。两组锚点沿用已记录的限定敏感性规则；属于发现后的规则修订，不称独立验证。
- T/NK作为包含先天淋巴候选的细分分支，B/浆细胞共同进入细分分支；粗分支不等于最终细胞类型。其余髓系、基质、上皮等独立，恶性状态另列。
- 唯一原生谱系与合格参考一致：记录双方支持。参考不可用时保留原生与局部群一致、或原生单独支持的不同证据层级；后者是低置信候选。
- 原生marker不足时，只有外部参考与来源内群体背景一致才允许低置信分配。参考无匹配/多匹配、未测与相反谱系证据分别记录。
- 原生跨谱系支持、参考之间跨谱系冲突或原生与参考相反进入待审池。群体多数不能覆盖单细胞矛盾；原生与参考一致的少数谱系保留自己的候选身份。
- 免疫参考的上皮/内皮/成纤维预测不用于非免疫定类；器官参考只在相应来源使用，正常器官标签不确定肿瘤恶性身份。参考间的一致不是独立投票。
- 原作者字段保留，不以旧标签代替作者真值。分配原因、原生/参考支持掩码、局部背景及旧新对应分别保存，权威标签保持。

## 来源内细分与跨来源整合

1. A：每个来源×新候选谱系分别建图和细分，局部群保留来源命名空间。跨来源细身份协调随后用原生证据完成。
2. B：复用原冻结联合整合的六档分群；输入细胞和对照成员不变。
3. C：在新候选谱系内联合建模；待审池保留来源内A的发现群，不强行跨来源对齐。

首个预定种子、共享模型特征、PCA/Harmony/近邻/Leiden参数和六档分辨率保持。此轮先分离路由/细分顺序，谱系内重选特征为下一独立因素。至少100细胞才拟合分支；不足的分支保留粗群并明确未完成细拟合，零RNA细胞仍在档案但不入表达模型。

所有正RNA细胞必须有发现或粗群归属；A群不跨来源/谱系，C已知群不跨谱系，待审池不跨来源。序列化分群逐项读回。4项实际反例覆盖边界、待审池、低数量保留和错位轴；粗分配另有6项实际反例及完整参考字典检查。

## 比较作用域

在旧粗谱系和新候选谱系两套分层分别报告，避免改变分层口径后只看有利结果。沿用有效群数比≤1.2的全部匹配对、相同查询细胞的表达单位留出预测与共同覆盖。分辨率组合不是独立重复，原生状态不当功能诊断。

发现包含全部试验来源，参考/粗身份可能编码RNA程序；这是条件性的生物结构保留诊断，不是独立身份准确率或新流程来源留出验证。新流程的原生细群证据、稀有/待审群处理、作者/正交证据覆盖及冻结参考映射仍需完成。不得依据临床P值或旧CM02选择流程。

## 来源用途与分母

新增来源级用途/覆盖表保留已有患者身份审阅、捕获、QC、参考/作者字段和测量状态。单位含患者键不自动等于身份已核；同单位多个键、缺失及部分覆盖分开。临床阶段/终点资格另审，既有已接受分析保持。

富集、产品、选择后的非恶性子集各自限制用途；“whole_tumor”或“tumor”描述不单独证明完整组织捕获。来源键唯一不能证明研究独立，已知可能重叠来源保留同组审阅。技术审计不自动授予新组成/临床分析资格。


## Fine-group review and numerical acceptance

The new source-lineage partition at the first frozen seed and resolution1 receives the same source-native full-gene expression review, with the two specified anchors. Explicit native author fields replace the old reviewer's carried-fine fallback; missing author labels stay missing. Each group carries its provisional native routes, number of cells with an actual fine-graph fit, and number in coarse-only branches. Own-group descriptive marker rankings are not replicate-aware DEG. Identity, RNA state and malignancy remain separate.

The independent annotation audit checks all source axes, branch-fit and author denominators, old-label crosswalks and reference codes; the smallest preselected source additionally receives direct CSR marker and official reference prediction replay. New-comparison acceptance checks frozen B membership, route/source boundaries, unknown-pool A/C equivalence, both lineage contexts, every matched pair and independent centered-OLS predictions in the size-selected source. These establish computational agreement; independent identity truth and new-route source holdout are separate requirements.
