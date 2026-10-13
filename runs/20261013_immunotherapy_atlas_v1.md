# 免疫治疗专用单细胞对象与下游 v1

- 授权日期：2026-10-13（UTC+8）。
- 当前阶段：来源范围和方法冻结；尚未声称新生物学生产完成。
- 输入：20261012_public_preprocessed_zarr_v2完整技术合集，Project_v5及历史免疫队列的既有QC/counts与作者临床资料。
- 来源角色：独立验证/内部复现/机制参照/产品分支/隔离逐原研究登记。不得将Spatial、动物、bulk计为单细胞；CAR-T产品和免疫富集单独标记。
- 输出：原项目卷新独立run下计数权威对象、降维/图/聚类、major/fine/state/单marker编号描述、逐样本CNV证据与精确回填、CM新发现及旧9CM固定投影；持久化前需阶段实测容量与完整验收。
- 注释：复用当前P3证据与命名方法，身份/恶性/状态分别保存，未决保留；不以单marker硬给成熟亚群。
- CNV：复用malignant_inference核心，原文库原生counts和同文库参考；skips不当normal；可用时既有CopyKAT证据只单独记录。
- CM：R CoVarNet新版本，rank survey与稳定性不强制K9；患者/材料/来源内共变，未捕获谱系为缺失。旧W固定投影与新W发现分开。
- 新增公共矩阵下载：本任务不恢复，优先现有Project_v5合格数据。新授权允许下游计算及独立版本。
- 详见[持续任务日志](../logs/codex/2026/2026-10-13_immunotherapy-atlas-downstream.md)。启动代码/参数/实际范围和终态按证据追加。

## 已执行节点

- 来源盘点完成：21条来源、5,274,578细胞登记，角色分层详见[汇总](../inventories/runtime/2026-10-13/immune_atlas_v1/source_preflight_public.json)。不等于临床验收或独立验证。
- GSE305078：原生count检查PASS，32文库DD/DecontX后台运行；独立复核pending。
- 下游pilot：2026-10-13 13:47:03 UTC / 21:47:03 UTC+8开始，PID 341235，12个免疫主队列/待审来源，每来源600细胞；冻结代码aef9684。范围之外的机制参照/CAR产品/基线重叠没有混入该模型。
