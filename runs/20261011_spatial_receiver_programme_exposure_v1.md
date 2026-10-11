# 固定接收程序与配体邻居暴露

## 预声明

2026-10-11T12:07:13.114311+00:00 / 2026-10-11T20:07:13.114311+08:00，父证据v15/0冲突/无活跃本阶段作业。[完整计划](../inventories/analysis/2026-10-11/spatial_receiver_programme_plan_v1.json)。上一轮是实质进展：全部45空间接口独立复算，29可估计、主BH无显著；不能用单个名义p宣称轴验证。

本轮保留15条固定记录（13接触、2程序专用），按家族/方向/LR/精确programme去重复为14定义×5文件=70主模型，映回75行；所有不可测/不可估计保留。先核查新增7基因与原1000基因面板，HBEGF/部分exhaustion基因缺失显式保留，不能只保留可算候选。

主问题是在R+接收者中，固定程序是否随局部L+发送者比例改变。分母是ALL发送家族邻居；没有发送家族邻居记缺失，有家族邻居但无L+为0。完整nativeRNA归一化，原固定程序不正式赋状态。原作者标签、FOV、双家族归属、RNA/检测基因/面积/密度/发送家族数量调整；强共线/非变异不自动删除协变量。

新有限niche尺度为双方面积圆半径之和×6；2/12作预先固定观察效应敏感性。这是接收表达问题，不能据其结果更改上轮接触尺度2的结论或称绝对微米/膜结合。R位置/程序固定，只在FOV、家族、R状态、作者细分类别、RNA/密度二档内打乱L标记。预定10000次，以免70个双侧检验仅2000置换的p分辨率限制单项校正；完整70主BH，重复原发现背景不计独立证据。

[NCEM原文](https://pmc.ncbi.nlm.nih.gov/articles/PMC10017508/)启发身份+niche预测表达；本轮只实现固定programme统计诊断，非完整NCEM复现。[NicheNet原文](https://pubmed.ncbi.nlm.nih.gov/31819264/)用于明确需靶基因支持，当前generic程序不能当已证实下游靶。Nature跳转500不可访问后改PMC/PubMed，失败与替代来源留任务日志。

## 验收和运行

生产及独立审计待代码冻结、真实执行；先保存预声明，后正常push/核验。输入/缓存/额外原生基因/归一化、图/分母、自接触、10000置换/直方图、全70BH及所有敏感性均复算。进程退出和计算PASS不提升正式身份/状态或CMdriver；独立患者和蛋白功能仍需要证据。完整细胞/FOV/置换数据仅留服务器。

## 代码冻结与合成验收

[冻结SHA及合成检查](../inventories/analysis/2026-10-11/spatial_receiver_code_freeze_v1.json)：原生全靶RNA分母、缺失非零、同细胞/跨FOV无边、强共线、SVD/QR和全部10000标签/null复算通过；真实70模型未开始。独立审计不import生产者，原生HDF5CSR额外基因、query-ball图、另编码协变量和QR回归核对。

```bash
PANCANCER_PREFLIGHT_ROOT="${RUN_ROOT}" python scripts/spatial_receiver_programme_exposure_v1.py
PANCANCER_PREFLIGHT_ROOT="${RUN_ROOT}" python scripts/audit_spatial_receiver_programme_exposure_v1.py
```

代码a7d065e97258222ff8255d4dcc00559a0fced613已正常push并独立HTTPS核验；生产2026-10-11T12:17:57.866679+00:00（20:17:57.866679 UTC+8）派发，实际开始/终态需读取execution.json，不能只依PID或派发宣称完成。

生产12:17:57.938231—12:18:56.150535 UTC（UTC+8加8小时）exit0，两进程均消失、0冲突；全部70主模型/75映射，22可估计、无名义p<0.05/主BH显著，尚待独立审计。源plan实际SHA4051c8dbd087e2cd72e9d0b9c03cc75b18bbea1728aeeca9c33e3ee03970bbd0与预先推送文件一致；创建器stdout哈希标签误指后续ROADMAP变量，已说明且不当输入变更。

独立审计12:22:09.052317—12:23:13.049368 UTC exit0；2026-10-11T12:26:04.291574+00:00取证四进程终止、0冲突。70主模型/75映射及140/70/70诊断、1,272,290原生额外值和150,000唯一标签draw复算通过，尚待附加留一FOV/映射验收及v16登记。没有名义p<0.05或完整主BH显著，未通过科学轴验证。

## 最终验收：当前状态

登记UTC 2026-10-11T12:29:43.368494+00:00 /UTC+8 2026-10-11T20:29:43.368494+08:00；本机取证UTC 2026-10-11T12:33:37.008282+00:00 /UTC+8 2026-10-11T20:33:37.008282+08:00。证据v16，0状态冲突、active_jobs为空，生产/审计四进程已终止。v16是证据登记版本；原172子群、9CM、v20261010、权威Zarr、原全局BH和历史产物未变。

| 项目 | 数量 |
| --- | ---: |
| 固定原记录 | 15（13接触、2程序专用） |
| 精确程序/方向/LR/家族去重复定义 | 14 |
| 主模型 / 原记录回映 | 70 / 75 |
| 可估计条件模型 | 22 |
| 作者未取得 / 作者无B类别 | 14 / 6 |
| 接收者不足40 | 17 |
| HBEGF / HLA-DRB1 / HLA-DRA配体缺失 | 4 / 4 / 1 |
| 程序面板不足 / 剩余自由度不足 | 1 / 1 |
| 未校正p<0.05 / 主完整70 BH<0.05 | 0 / 0 |

### 问题和结论

这次问的是：在同一身份内的受体阳性细胞中，周围发送家族的配体阳性比例增加时，接收细胞的固定表达程序是否随之改变。它与前一轮“双方是不是更常相邻”是不同问题。

发送家族邻居为0时，暴露缺失；有发送家族邻居但没有配体阳性邻居时，暴露为0。接收表达用完整原生RNA分母归一化，调整RNA、检测基因数、面积、密度、发送家族邻居数量、原作者类别、FOV及双家族归属。主邻域6、诊断2/12在结果前固定；10000次只打乱配体标签，受体位置及程序不变，完整70主检验校正。

**当前固定通用程序没有支持配体邻居效应。** 22项中最小未校正p约0.211579，全部主BH为1。没有根据阴性结果修改程序、半径或端点，也没有把缺失填成阴性。该结果不证明生物通讯不存在；通用程序不等同受体的特异靶基因，1000基因面板、作者粗标签、真实患者/core映射缺失仍限制解释。

来源/癌种限定的细胞对生态仍保留价值，不要求所有癌种完全一致或完全脱离组成。此次家族扩展与原CM直接边匹配为0，不能称原9CM驱动机制；当前亚型/marker身份没有被替换为正式状态。

### 时间线及计算验收

| 节点 | UTC | UTC+8 | 结果 |
| --- | --- | --- | --- |
| 预声明 | 12:07:13.114311 | 20:07:13.114311 | 固定70/75范围及参数 |
| 生产 | 12:17:57.938231—12:18:56.150535 | 20:17:57.938231—20:18:56.150535 | exit0 |
| 独立复算 | 12:22:09.052317—12:23:13.049368 | 20:22:09.052317—20:23:13.049368 | exit0 |
| 补充验收 | 12:27:44.531602—12:28:42.017650 | 20:27:44.531602—20:28:42.017650 | 全70一致性PASS |
| 登记v16 | 12:29:43.368494 | 20:29:43.368494 | 0冲突 |

生产和独立审计代码预先推送a7d065e97258222ff8255d4dcc00559a0fced613。审计没有import生产者：原生HDF5 CSR额外基因1,272,290值、query-ball图、另编码协变量及完整QR系数、全部10000条件null核对；15组唯一配体draw共150000次重放，22模型各10000统计量核对。共享预先规格/种子，近浮点并列用已保存精确draw核对p；不是独立生物复现。

140半径、70命名marker剔除、70受体表型诊断都复算，诊断没有新增推断p值。补充验收复用已审计的独立图/design/QR函数核对全部70映射、分母、覆盖、诊断状态和留一FOV；不冒称第三套独立生物验证。support代码是在执行后归档，时间区别保留。

### 失败/修正和未验证项

创建器stdout曾错标计划SHA，实际计划源及预先云端文件均为4051c8dbd087e2cd72e9d0b9c03cc75b18bbea1728aeeca9c33e3ee03970bbd0；错标原因是后续ROADMAP复用了变量，文件与方法没有变化。工具展示输出截断后重新只读观察终态，没有重复运行生产。当前生产/独立复算/补充验收无失败；历史输入/作者映射失败继续保存。

完整细胞/FOV/标签/置换矩阵和真实患者表不上传。公共主表/半径表删除精确covariate_names类别列，保留列数与固定schema；[变换及源SHA](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/privacy_transforms.json)可核查。图像实际像素审阅PASS。当前未确认独立患者、蛋白结合、下游靶响应或临床效应。

### 产物

- [70主模型](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/all70_primary_programme_exposure.tsv)、[75映射](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/all75_record_mapping.tsv)
- [140半径诊断](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/all140_radius_sensitivity.tsv)、[70 marker诊断](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/all70_programme_name_sensitivity.tsv)、[70受体表型](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/all70_receptor_phenotype.tsv)
- [审计](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/independent_audit.json)、[补充验收](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/accepted_primary_consistency_checks.json)、[终态及SHA](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/checkpoint.json)、[时间线](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/timeline.tsv)
- [覆盖与效应图](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/figures/all70_receiver_programme_exposure.png)、[图像审阅](../inventories/analysis/2026-10-11/spatial_receiver_programme_exposure_v1/figures/visual_review.json)

### 后续有界任务

先审阅已有失败的身份/来源/机制分支，再核查完整基因的原生及免疫治疗队列、真实患者区组和发现集重叠。检索受体对应靶基因/TF先验，记录确切版本、方向、复合物和覆盖；覆盖不足不能自动换端点。新的推断范围、完整校正家族及独立验证单位须预先固定，不能反复更换当前通用程序来追求显著。

当前目标仍active，未找到可以宣布验证完成的通讯轴。本节点只是该固定问题的计算和验收闭环，后续继续细胞对生态及机制响应探索。
