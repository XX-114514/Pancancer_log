# 巨噬—T细胞跨来源响应探索 v1

## 状态

方案冻结，真实运行尚未开始；不是独立验证，也不替代原9CM全局阴性。

## 问题与范围

按原发肿瘤跨癌种覆盖选择巨噬↔T两个方向，5个GSE、肺/乳腺/结直肠/胶质瘤/胃；CRC分选与未分选情境分别调整，但同作者患者跨情境共用块。全4624LR、外部受体先验匹配全部TF，至少2个GSE具备原固定设计覆盖；每变体9248 LR分母。原4,676,787细胞/172节点/9CM与版本v20261010不改。

## 方法

每情境标准化L和TF；拼接全部原丰度/细胞数/非LR复杂度/受体表达调整设计，各情境保留独立完整nuisance span。CR2/Satterthwaite共同斜率，effective df≥4；亚型混合PC及全基因TF为敏感性。TF家族Bonferroni联合，再全9248 BH探索/BY敏感性；无效内部p1、对外缺失。输出各GSE方向、信息贡献、逐GSE删除结果。共同斜率不表示各癌种效应相同，来源限定关联保留价值。

## 验收与边界

全family/聚合/哈希检查；两个方向×三个变体各4个输入组共24组，直接全设计矩阵复算 pooled/source/LOSO CR2。不会声称所有系数完全独立重演。lead TF区间为边际描述，非选择后校正；既有TF近似检验校准限制、基因背景差异、未知分选/重复样本/样本替代患者继续保留。真正外部验证、受体结合及扰动证据未完成。

## 冻结文件和执行

[参数与SHA](../inventories/analysis/2026-10-12/macrophage_T_multisource_response_plan_v1.json)，生产/审计代码位于scripts同名文件。沿用scanpy环境Python，OPENBLAS/OMP各2线程；真实根目录通过PANCANCER_PREFLIGHT_ROOT、PANCANCER_ANNOTATION_ROOT提供；PANCANCER_MAC_T_POOL_PLAN、PANCANCER_CR2_CORE、PANCANCER_FROZEN_LR_RESOURCE指向冻结入口。执行`python -u ${CODE_ROOT}/macrophage_T_multisource_response_v1.py`后`python -u ${CODE_ROOT}/audit_macrophage_T_multisource_response_v1.py`；独立execution manifest，不改旧primary registry。

时间：2026-10-12T20:37:56.525513+08:00 /UTC 2026-10-12T12:37:56.525513+00:00

## 终态

真实生产与审计均exit0；完整汇总、时间、SHA/范围和限制以[结果](../inventories/analysis/2026-10-12/macrophage_T_multisource_response_v1/RESULTS.md)为准。
