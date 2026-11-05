# 服务器分析系统梳理与GitHub交接

## 用户要求、起点与范围

用户要求详细研究本服务器已做分析，系统梳理并上传既有GitHub仓库。源main为`e6cf388ff5ca53c4e00cd4fcaa7559d72028a71e`，原`nmf_benchmark`有11个未跟踪脚本。源仓库不在可写根，建立可写隔离克隆；不修改原工作区、冻结对象、任务或用户原有改动。目标沿用已授权`XX-114514/Pancancer_log`，普通push，不强推或重写历史。

取证起点2026-11-05 00:49:42 UTC / 08:49:42 UTC+8；原事件分别按各run，不以本次mtime冒充运行完成时间。审计范围限定本仓库、独立investigations、早期Xenium试点、指定PC_ROOT及labshare收据，没有检查其他用户或重跑生产分析。

## 计划与读取

先读取AGENTS、README/STATUS/ROADMAP/TODO、DATA_POLICY、WORKFLOW、CODEX_CLOUD_DEVELOPMENT、RUNTIME_SYNC及WEHI工作流，近日日志、相关方法/run/release，后逐项比较工作区与仓库。重点读V5/V7聚合比较、V7 CM/blisa原报告/脚本/配置/完成摘要、NMF 45任务与匹配实现、候选路由试点/preflight/测试、TME/去批次/来源登记、CT75重放及空间报告、早期PDAC QC和迁移问题表。只读取所需小摘要，不复制大型文件。

## 实际执行

1. `git status --short`、`git log`、remote/worktree检查；完整克隆后以原作者仓库身份配置副本，目标remote沿用源仓库GitHub SSH别名。
2. `collect_snapshot.py`按白名单读取206项证据，补录69项配置/原W hash/迁移问题表；统计1,179运行目录和16分析线，区分历史迁移/重放/新计算/预备方案。
3. 归档43脚本：原11 NMF＋32工作区脚本。公开物理根改为环境变量，NMF新运行拒绝已有输出、smoke不再删除旧目录；不改数学拟合逻辑。逐一复核原代码字节hash43/43一致，公开副本hash另登记。
4. 读已有45个小W表，做列尺度规范与等价缩放反例，写`nmf_scale_sensitivity.tsv`。没有训练或重新拟合，不公开样本H/W矩阵。
5. 写系统报告、版本/工作线/纠正表，更新README双语入口、STATUS、ROADMAP、TODO、CHANGELOG、run/daily和相关索引。
6. 显式stage本次范围，审阅差异、全部归档代码路径改动、表头/汇总语义、体积与敏感路径；未发现新增内部物理路径、逐细胞/患者表或大对象。

## 发现与结论边界

- V7主/legacy rank survey与K9已有完成报告，25µm C05已有摘要，早期pending应保留为历史。K9是override；K4的完成结果未建立。
- 扩展72来源恶性侧表都存在，但没有重做其行级验收；历史4来源等待侧表已过时。全量路由收集根/完成收据未查见。
- 当前TME范围3,347,659，CoVarNet v1合格2,980,723；旧对象保留范围外行。run_v2只核到Harmony完成，最终summary和调度器退出未核。
- Scanorama重试已有图指标，但缺scIB行。CT75主CM未空间验证；V7弱信号依赖成员定义且目标myCAF身份不合格，不能宣布CAF机制或CM验证。
- NMF原成员余弦受等价因子尺度影响；规范后部分匹配改变，主要跨seed稳定性趋势保留。参考nsNMF不是金标准，不能将不同输入/预处理结果当纯优化器准确性排名。
- labshare指定旧对象存在，15旧迁移问题表均0问题行；本次未重hash大数据或验证全部环境迁移。

## 验证、失败与修正

- 新代码Python AST 23份、shell `bash -n` 13份，R `parse()` 10份均通过。
- marker路由4项调度mock＋3项合成收集测试全部通过；没有真的提交Slurm任务。
- `verify_snapshot.py`通过：64个公开证据hash引用、43代码hash、45个唯一NMF任务/9方法/K9，两个转移表各4,676,787，试点11,381路由守恒及未冻结身份标志，rank 2–20覆盖，空间各29份摘要，尺度反例及独立符号检验算术。
- 仓库validator实测0 errors / 1项既有绝对路径warning；release validator 0 errors。新增文件敏感路径扫描无命中，单文件均小于5MiB。
- 首次stage后`diff --check`发现TSV末列空值形成尾tab；改为明确NA、保留矩形列数并重算公开hash后通过。原SHA和科学数值不改。CSV写出LF，避免Git归一化导致交接差异。
- 首次隔离工作目录调用登录shell遇到UID无用户名提示，使用非登录shell进行本地文件操作。真实GitHub `ls-remote`仍在SSH建连前失败：当前NSS查不到numeric UID，报`No user exists for uid`。此前HTTPS/DNS检查也失败；没有修改身份、认证、代理、系统或网络设置来绕过限制。

## 传输状态与后续

科学整理材料已就绪；普通commit/push和独立远端核验以随后实际收据为准。目前没有GitHub上传成功证据。源仓库原有改动始终保留，未新开后台任务。内容建议commit：`Audit server analyses and archive reproducible evidence and workflows`。

下一步先TME run_v2终验，再刷新候选路由预检；NMF按固定尺度和rank/来源留出复核，空间补身份/面板/患者单位及多重比较资格。详细下一步已写TODO顶部。未验证真实逐细胞细分流向、全部bin级检验族、独立身份准确率、新V7 LIANA/响应端和调度器最终退出；未验证内容不补造。

本次整理日志更新：2026-11-05T00:59:44.717982+00:00 / 2026-11-05T08:59:44.717982+08:00。

## 本地提交与推送失败收据

内容提交：`fc136204796fb59368c8230941b0b31bea7b8457`，普通commit成功，140文件/6,294新增行。最终仓库validator再次0 errors / 1项既有warning，release validator 0 errors，diff及离线快照核验通过。

在2026-11-05T01:01:09.759495+00:00 / 2026-11-05T09:01:09.759495+08:00执行普通`git push origin HEAD:main`，退出128，SSH在建连前报`No user exists for uid`；没有上传成功，也无法独立核对GitHub main。不是已知远端仓库故障或密钥失效。未尝试改NSS、身份、主机校验、网络限制或持久代理。

状态：本地已提交、待GitHub同步。可写工作区另备普通终端发布脚本、增量Git bundle及传输收据；脚本核对目标、干净状态、预期提交和远端祖先关系，普通push后独立ls-remote比对SHA。若远端出现独立新提交则停止，不强推。源仓库及科学对象保留原状。
