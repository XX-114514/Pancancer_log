# V7 全 panel 修复与文献优先算法

## 用户请求与初始状态

用户要求推进整套基因 panel 修复、显式输出并尽量逐项附文献，算法优先采用文献逻辑。上一轮已完成 162 旧 panel/17 父类审计。本轮读取最新 Mast/Tcell_6 子集证据：Tcell_6 存在 CD4/CXCL13 等非同质程序，不能按旧 Treg 候选或 CellTypist 多数票整体转换。

已读仓库规则、README/STATUS/ROADMAP/TODO、前轮方法及最新 T6 run 证据。主目录大量既有修改保留；从最新 main 创建独立工作副本。旧 V7 标签保留，新资源与算法写新版本。

## 初始计划

1. 按论文和作者代码重建 parent/branch/state/isotype/boundary 的新 registry，逐条输出正向/对照/命名要求、适用范围及来源；旧 162 项逐条迁移，不默默丢弃。
2. 实现文献引导谱系内清理→重建→RNA DEG→marker/状态独立读出的流程，以及非赢家独占的候选证据引擎。参数属于作者实例、算法默认或项目选择时分别声明。
3. 使用已核验 V7 实测聚合表达/DEG 检验 T4/T6、最新 T6 分群与旧 40 群；未覆盖基因/对照不假装有结果。
4. 合成边界测试、实际回放、显式表格与文档、正常同步并独立核验 SHA。若完整新细胞级流程未执行，则明确边界。

## 已执行与结果

已读取Lodi谱系内作者代码、Chu功能签名及此前595条作者材料目录；新核查Ma DC正文、Rao Tph原始研究、Tirosh melanoma作者PDF、Aizarani肝图谱及Henry前列腺正文。原始材料出处及适用边界见资源来源表，不把DEG摘录当作者完整分类器。

写入resources/annotation/v7_panel_repair_v1：112模块、409基因、旧162全迁移；634条基因依据中266条有已核验材料出现证据，368条仍标记逐基因参考复核。3个模块保留reference_scope_gap。43项contamination全部退出赢家排序。

写入scripts/v7_panel_repair_v1：全候选证据引擎、Seurat谱系重建、合成fixture/监督执行/独立验证、真实聚合回放及15项反例测试。全体候选逐轴检查，不生成最终标签；alias碰撞、未测分母、命名marker与DE定义均显式处理。

执行过的主要命令为`test_evidence_engine.py`、`run_smoke.py`、`replay_aggregates.py`及工作簿构建；完整可复现入口见脚本README/run。服务器合成测试06:12:28–06:12:47 UTC / 14:12:28–14:12:47 UTC+8三个阶段exit0，独立核对1227表达行/6Treg同细胞记录通过。真实40群4480条候选、新T6 12群1344条候选均回放完成；范围重叠不相加。T6 top10 DEG按截断输入处理。

使用Spreadsheets技能生成4表Excel并逐表渲染审阅：112panel、162迁移、634逐基因依据、21来源范围；公式错误扫描0。方法、资源README、run/daily及必要状态索引同步维护。

## 失败与修正

部分文献URL404/recaptcha，使用可核实原始来源替代，无法逐项核验则保留缺口。首次R异步试验未记录子进程退出码，随后受监督重跑合成QC场景保存真实终端。组件可测字段原命名过宽，拆开“每组件有一可测基因”与“全部核心均可测”并加入反例；DE别名处理和分母非有限值拒绝亦加入测试。

独立工作副本partial clone快进第一次因受限网络无法补取对象失败，获既有授权的联网执行后成功；保留远端并发TME提交及主工作区全部无关修改。

## 未验证与下一步

新409基因全V7提取、真实各谱系全矩阵重建、来源/患者留出复现和最终人工命名未完成。R UCell未安装，该可选分支及可选Harmony分支未在本轮smoke执行。既有CellTypist/scATOMIC是旁证，本轮没有重跑。大类证据充分时可以停大类，不能用本次资源修订宣称全图谱发表级注释完成。

下一步从T/NK/ILC共同边界、T6 helper/Treg及混合RNA开始，逐谱系复核。建议commit：`Repair V7 panels with literature-guided multi-axis evidence workflow`。

## 发布核验

正在独立干净工作副本进行隐私/diff/仓库轻量检查；发布状态以追加的远端SHA核验收据为准。历史V7和并发修改均保留。

发布前检查：仓库轻量验证0错误、1项既有历史内部路径警告（非本轮文件）；新增文件路径/凭据扫描及工作簿XML审查未发现私有路径。首次暂存检查发现ALL_PANELS末尾多余空行，已修正；新TSV统一LF后重算资源SHA。

独立发布核验（2026-11-05 06:20:28 UTC / 14:20:28 UTC+8取证）：内容commit `d61e0d9e0b0eea08f5f7e3b22128b569e785ee85` 已普通push到用户指定public仓库，随后独立HTTPS `git ls-remote origin refs/heads/main`与本地HEAD完整一致。35个文件；仓库验证0错误/1历史警告、暂存diff检查通过。此核验只证明资源/代码/汇总已发布，不改变全V7最终注释未完成的状态。
