# WEHI R4 / CoVarNet K15：完整任务与结果归档

本页是 **2026-11-15 WEHI R4 / K15** 模块的结果快照，不与 0906 K9 或 2026 年 10 月模块编号互换。2026-11-18 扩展目标仍在进行中。2026-11-19 用户明确指定本仓库，并接受立即公开研究结果。

- [完整任务、探索与结果总结](完整任务与探索总结.txt)：覆盖前置接收、所有分析线、阴性结果、失败修正、版本替代和未完成工作。
- [逐轮进展记录](results/CM_validation_20261118/进展索引.txt)：历史状态保留；旧 q 值须按总结中的最新版家族解释。
- [当前证据整合](results/CM_validation_20261118/integrated_evidence_review/结果整合.txt)。
- [Visium 与前置成像报告](results/WEHI_R4_spatial_20261115/review/report.html)。
- [代码模板](code_templates)：原分析脚本仅将机器专属路径替换为变量占位符，原始 SHA 见清单；需在科研环境配置路径后运行，不是无数据即可重现的独立软件包。
- [本地交付清单](local_publication_candidates.tsv)、[服务器全文件清单](server_file_inventory.tsv)、[服务器完整快照回执](server_snapshot_receipt.json)。

## 已完成的主要结果

| 分析线 | 范围 | 当前结果 |
|---|---|---|
| 通讯检出 | 全 15 CM，25,912,935 检验 | 联合校正与审计完成；大部分核心 LR 跨 CM 共享，不能当差异通讯 |
| high/low 表达 | 3,253,817 细胞级/单元级检验 | 联合校正完成；上游独立审计 14/15，CM15 尚在运行 |
| 患者级差异 LR | 三来源 45,821 非恒定检验 | 联合校正无 q<0.05；核心血管对仅肾癌满足患者数资格 |
| Visium | 4 切片、10,160 spots、15 CM | CM12 核心调整相关 0.443–0.569；CM01 较弱 |
| 黑色素瘤 legacy ST | 8 切片×N5/30/60 | 18/24 完整模型；全8个N30核心正相关，完整三方hub不一致；全N比较待完成 |
| 成像 | 乳腺与新肺癌共466合格条件检验 | 部分乳腺设置支持邻近富集；肺癌未通过联合校正 |
| 免疫浸润 | 9120非Redacted TCGA评分患者；5400检验 | CM12方向依癌种而变；参考分解主要反映壁细胞相关表达 |
| 预后 | 3690常数模型、3150分段模型；13140检验 | CM12为模型特定KIRP线索；CM13有LIHC时间效应线索，未外部验证 |
| ICI响应 | Hugo、Riaz、Rose；480检验 | 无联合q<0.05；Rose取材时序未知 |

## 当前版本与文件解释

优先使用 nonredacted TCGA 结果、466 成像检验家族、13140 生存检验家族、480 ICI 检验家族，以及全15 CM联合校正。历史目录保留供追溯，不能从旧家族挑选显著结果。

所有原生脚本、输出、数据与模型在科研服务器保存。完整现有结果快照包含2377文件、1,258,586,452字节，SHA256 `c8bc3998a97a7218566acca05e2c38fcbcc54493321948e30054c370d3aa457c`；独立逐成员校验通过。原始输入与下载来源另列清单，未重复打包；不跟随外部链接。

本公开归档收录422项本地交付与1233项服务器统计导出，并附总结和清单；全部1009个表的10,841,255行通过结构检查。它包含审查后的代码、汇总统计、图和审计投影。逐患者/细胞记录、原始输入、大型模型对象及原始完整运行日志不在GitHub；文件清单逐项区分，不能将此称为服务器字节级全镜像。机器专属路径被占位符替换，相关JSON是审计回执的发布投影；其原始SHA与发布SHA分别记录。

[server_statistics](server_statistics) 提供完整统计表；大于单文件限额的压缩表按 `.part000` 等分块，依序拼接再解压，整体SHA见 `aggregate_export_manifest.json`。这些是统计组合记录，不是独立患者行。当前快照不含未来尚未完成的结果。

## 待完成

CM15 原始 high/low 审计、余下6个N60模型及全部24模型的匹配敏感性比较、最终全目标验收。病理/分割图像核对、未核实来源的患者/基因轴和临床预测外部验证仍有证据缺口。空间邻近、通讯检出、表达差异与预测效果各自解释，不相互替代。

## 主要图表入口（最新版）

- [全15 CM核心通讯事件及共享性](results/CM_validation_20261118/communication_joint_review/core_pair_event_counts.png)
- [八张legacy ST的CM12三条边](results/CM_validation_20261118/integrated_evidence_review/CM12_three_edges.png)
- [乳腺周细胞身份敏感性，最新466家族](results/CM_validation_20261118/lung_core_review/breast_mural_latest_q.png)
- [肺癌条件邻近敏感性](results/CM_validation_20261118/lung_core_review/proximity_sensitivity.png)
- [全模块TCGA免疫浸润，nonredacted版本](results/CM_validation_20261118/immune_review_nonredacted/all_modules_purity_sensitivity.png)
- [全模块参考贡献分解](results/CM_validation_20261118/component_review/all15_reference_attribution.png)
- [预后时间效应候选及全部调整](results/CM_validation_20261118/survival_review/piecewise/time_change_candidates_all_adjustments.png)
- [三个ICI队列全部模块对照](results/CM_validation_20261118/ICI_review/three_cohort/three_cohort_descriptive_comparison.png)
