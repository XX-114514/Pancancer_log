# 巨噬—T固定接收程序验证探索 v1

## 问题与范围

前轮为真实进展：全范围共同程序控制后主/亚型0，IFNG→IRF1仅全基因敏感性q=.0497。不是完成轴验证。本轮保持五源12原发情境、全4624LR两个方向；固定四外部程序：IFNγ、IFNα、TNF/NFKB及γ减α共享基因，共36992问题/变体，四变体147968完整分母。原9CM/172身份/权威Zarr不改，用户授权继续；另一个任务的CM恢复结果保留，不将骨架诊断误称独立复现。

## 读出与控制

官方人类Hallmark字典200/97/200基因，γ与α重叠73；γ包含44 IRF1参考靶基因，53全局LR基因；五源共同符号匹配γ196/α95/TNF198/γ减α124。未匹配可能是符号别名，不写成真实未测。使用固定共同基因平均log1pCP10K并源内标准化，每个被检LR去掉自己的基因；不是GSEA/GSVA或正式状态。敏感性：全局LR排除、IRF1靶排除、亚型PC。γ减α不证明γ生物特异性。

完整原丰度/细胞数/非LR复杂度/受体控制，加发送HSF1和接收自身配体RNA；原LR可识别设计与外部先验门槛保留。程序≥20基因、≥80%库有≥10基因在至少1%接收细胞中检测，且比例有变异；至少2GSE。CR2有效df≥4；所有四读出统一完整36992 BH探索/BY。来源方向、信息贡献和LOSO另报，仍不是独立验证。

## 验收/执行

现有缓存SHA与sample/major键检查；完整family及数值24均衡输入组，加预先固定IFNG/HSPA1A/LGALS1/TGM2接口额外fixtures，直接full-X pooled/source/LOSO复算。仅这些数值fixture复演，不称所有系数独立审计。

[冻结参数和SHA](../inventories/analysis/2026-10-13/mac_T_receiver_program_plan_v1.json)、[共同基因和重叠字典](../inventories/analysis/2026-10-13/program_reference_coverage_v1.json)。沿用服务器scanpy环境，BLAS/OMP各2线程；根路径通过PANCANCER_PREFLIGHT_ROOT/ANNOTATION_ROOT、CORE/LR_RESOURCE提供，PANCANCER_RECEIVER_PROGRAM_PLAN/PROGRAM_DICTIONARY/EXPLICIT_POOL_AUDITOR指向冻结入口；分别执行mac_T_receiver_programs_v1.py/audit_mac_T_receiver_programs_v1.py。新execution manifest，不改旧primary registry。

## 当前状态

冻结未启动，后续真实执行/审计/科学结论分开。参考JSON公开下载已成功，网页相对下载链接首次404，改官方human入口；非表达数据下载。原参照来自[官方γ程序](https://www.gsea-msigdb.org/gsea/msigdb/cards/HALLMARK_INTERFERON_GAMMA_RESPONSE)，公布网址和字典SHA，不复制论文全文。

事件时间：2026-10-13T20:25:28.321502+08:00 /UTC 2026-10-13T12:25:28.321502+00:00

## 终态

生产和数值audit均exit0，详见[真实结果/范围/下一步](../inventories/analysis/2026-10-13/mac_T_receiver_programs_v1/RESULTS.md)。科学轴仍未验证。
