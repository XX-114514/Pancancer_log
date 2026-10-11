# 固定通讯接口的条件空间邻近检验

## 预声明阶段

创建时间：2026-10-11T11:43:34.447250+00:00 / 2026-10-11T19:43:34.447250+08:00。状态：predeclared_not_started；父证据v14，未改注释v20261010/172节点/9CM或全发现BH。

问题：在同一视野、相同细胞家族位置、RNA深度与局部细胞密度条件下，固定L+发送者和R+接收者是否比随机标签更常邻近？范围为全部9个唯一接口×5文件=45主检验，回映13记录×5=65行；不因低覆盖更换端点。

[完整固定计划](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_plan_v1.json)。原生空间v3缓存/输入已独立审计；新目录保留所有历史。新检验用作者粗T/B谱系+背景口径，缺作者/无B类别/缺基因记NA；RNA纯谱系和count≥2另作描述敏感性。

几何采用原生像素和面积推算的圆半径，距离≤2倍双方半径之和；1/3倍半径敏感性。这个定义是有限邻域，不能称真实膜接触或绝对微米标定。仅同FOV，排除自接触。条件零模型在FOV×双家族归属×RNA深度二档×全QC细胞局部密度二档内打乱联合L/R两位标签，保留两端数量及同细胞共表达。固定2000次、seed20261011；文件内FOV合计，主单侧富集p按完整45家族BH，NA仅校正内部p=1，输出仍NA。

≥20/FOV原门槛是可测性描述，本轮不以此删掉稀疏可测关系；报告null方差、可移动标签数、FOV集中度/留一FOV。患者映射未知，不能把FOV或TMA当独立患者。[原补充表结构审阅](../inventories/analysis/2026-10-11/spatial_core_supplement_review_v1.json)没有完整FOV—患者桥；原Excel/RDS/细胞行不上传。

## 执行与验收

生产代码、审计代码和正式时间后续记录；现在没有运行结果。验收包括全45/65映射、源SHA、原始几何独立计数、端点/重叠/条件null全量复算与合成边界检查。审计PASS不代表患者级验证、蛋白结合、下游反应或原9CM驱动。APP–NOTCH2原文主要支持同细胞蛋白关联，HLA-DRA是完整MHC-II代理，其余结合/糖链/细胞方向限制完整保留在计划；不能把统计邻近直接写成免疫抑制机制。

可复现入口使用环境根目录；正式代码提交后补齐调用。关键/探索节点按持续[任务日志](../logs/codex/2026/2026-10-11_milestone-log-and-cloud-sync.md)普通push并独立核验SHA。

## 代码冻结检查

生产/审计代码[SHA及合成验收](../inventories/analysis/2026-10-11/spatial_conditional_code_freeze_v1.json)保存；两份代码独立实现图与统计量，未import生产者进行审计。服务器2000次合成重放、重叠/自接触/缺失均PASS。真实执行尚未开始。

```bash
PANCANCER_PREFLIGHT_ROOT="${RUN_ROOT}" python scripts/spatial_conditional_neighborhood_v1.py
PANCANCER_PREFLIGHT_ROOT="${RUN_ROOT}" python scripts/audit_spatial_conditional_neighborhood_v1.py
```

原生细胞缓存、FOV/置换细表只留服务器；公开45接口汇总。预声明已普通push并独立核验，实际代码提交随后冻结。

生产开始：2026-10-11T11:49:22.700694+00:00 /19:49:22.700694 UTC+8；代码93ec55711f6accc113660bdeb8d126a99f540405。11:50:05取证worker/child活跃，新结果持续写入。结束/退出状态及审计待确认。

生产结束：2026-10-11T11:50:17.798455+00:00 /19:50:17.798455 UTC+8，exit0，两进程均终止。全部45主检验，29可估计、16不可测（9未获取作者、6无B类别、1缺基因），主BH<0.05为0；尚待独立审计。没有将不可测项填阴性或把未校正趋势当验证。

## 最终验收（以本节为当前状态）

取证2026-10-11T11:58:08.699340+00:00，证据v15于2026-10-11T11:54:05.524551+00:00登记；0状态冲突，生产/审计4个进程均终止。原v14已验证端点与全部历史保留；未改172节点/9CM/全发现BH或权威Zarr。

| 项目 | 结果 |
| --- | ---: |
| 固定唯一接口×文件 | 9×5=45 |
| 原记录回映 | 13×5=65 |
| 可估计条件检验 | 29 |
| 作者未获取 | 9 |
| 作者无独立B类别 | 6 |
| 缺HLA-DRA | 1 |
| 主完整45 BH<0.05 | 0 |
| 未校正p<0.05 | 3 |
| 独立FOV置换重放 | 7,074,000 |

### 结果含义

原先每FOV双方≥20的预检门槛并不等于其余关系无法检验；固定条件null使29项可以估计。16项是不可测，保留NA，不能填成阴性。主校正0并不证明没有真实通讯，只表明本轮固定表达门槛、几何和条件null没有提供可发布的显著空间富集。

SELPLG–SELL：LUAD TMA1/2、MESO TMA1/2主邻近比1.035/1.117/1.922/1.102，四文件方向均正，但MESO1最低p0.004998对应q0.224888；其留一FOV超额最小仍+3.506，其他文件方向较弱。更小邻域LUAD TMA2比0.900，方向并非对所有尺度稳健。MESO1 RNA纯谱系为47发送者/7接收者、0邻近对，作者count≥2为54/5、1对；身份/表达门槛敏感。保留为来源限定的黏附/亚型或程序线索，不能称跨癌种验证成功、耐药/抑制机制。

另两个未校正趋势为MESO2 CD55–ADGRE5（2.247倍、p0.01999/q0.389805，留一FOV超额最小−0.338）和LUAD1 APP–NOTCH2（1.753倍、p0.025987/q0.389805）；前者集中度限制、后者cis/trans结合限制未消除。LGALS9两接口有端点覆盖但本轮未达到空间富集验收；不能用覆盖替代机制证据。

### 生产和独立复算

生产UTC 11:49:22.700694—11:50:17.798455；审计UTC 11:51:08.444307—11:52:24.836873（UTC+8均加8小时），exit均0。代码预先普通push并独立核验93ec55711f6accc113660bdeb8d126a99f540405；输入/plan/代码/审计SHA见[终态](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/checkpoint.json)。

审计不import生产者：query_ball_point重新构图+CSR双线性计数，对照生产query_pairs/边求和；小合成边界及每文件第一FOV直接距离矩阵复核。全45×2000联合标签、每FOV每stratum两位直方图守恒、零自接触/零跨FOV边、全部135半径与135表达口径均核对。验收阶段另检查主表等于半径2统计、留一FOV/集中度，未将同规格随机种子复算冒充独立生物学验证。

### 产物和证据

- [全部45主结果](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/all45_interface_results.tsv)
- [65记录映射](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/all65_record_mapping.tsv)
- [全部135半径敏感性](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/all135_radius_results.tsv)、[135身份/表达敏感性](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/all135_endpoint_sensitivity.tsv)
- [独立审计](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/independent_audit.json)、[验收主表/留一复核](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/accepted_primary_consistency_checks.json)、[登记](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/registration_receipt.json)、[时间线](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/timeline.tsv)
- [邻近图](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/figures/all45_conditional_neighborhood.png)、[图源](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/figures/figure_source.tsv)、[像素审阅](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/figures/visual_review.json)

图v1“B identity absent”措辞会误导为细胞不存在，v2类别说明拥挤；各版本保留，v3准确拆行且深色文字对比修正，像素审阅PASS。公开只汇总/源SHA/代码，完整矩阵、工作簿、作者RDS、细胞键、FOV与置换细表未上传。

### 后续有界任务

身份和连续表达程序分别建模；SELL+ IL7R低细胞毒程序仍是身份内部表达现象，不能仅凭它解释为信号抑制。下一步预声明接收者programme与L+邻居暴露的关系及患者/core或独立原生/免疫治疗数据验证，同时报告家族丰度、RNA深度、密度与所有不足。来源/癌种限定效应仍保留，不要求所有癌种同效应。当前家族扩展与原CM直接边匹配为0，CMdriver需另设计。整体通讯轴验证目标仍active，未标complete。

### 原文核查索引与表字段解释

[固定9接口原始文献/结构索引](../inventories/analysis/2026-10-11/spatial_conditional_neighborhood_v1/primary_binding_sources.tsv)补录本轮及父轮已核查来源，用于限定结合/糖链/细胞方向语义，不作为本空间实验的新验证。APP–NOTCH2[原始研究](https://pubmed.ncbi.nlm.nih.gov/16118793/)的同细胞蛋白关联不能自动解释为跨细胞激活。

不可测行的端点阳性、observed_pairs、p/BH均为NA；family_union/movable_cells/paired_positive_FOVs中的0是缺身份/基因时计算空掩膜的占位，不能解读为实际无该家族细胞。必须连同status读取；本轮科学图与显著性没有把这些占位当阴性。

## 云端发布

方案、代码、生产、验收4节点均普通push/独立核验；最终验收内容1319ad60a82a7e845e4a4fc8e28a2e98e8bc5ee6于2026-10-11T12:00:20.718564+00:00（20:00:20.718564 UTC+8）push，12:00:27.061812 HTTPS核验同SHA。[收据](../inventories/analysis/2026-10-11/spatial_conditional_acceptance_git_receipt_v1.json)。收据追加提交另核验；没有force/历史重写、患者行或矩阵上传。
