# 迁移后项目整理与下一阶段运行计划

## 用户要求

根据此前对pancancer项目的计划，先探索目录、整理项目，再规划下一步各任务的运行计划。取证时间2026-10-31约01:30–01:50 UTC / 09:30–09:50 UTC+8 / 12:30–12:50 Melbourne（AEDT）。

## 开始前状态

仓库工作树干净，本地main与origin/main同为`2f2b289`。

## 检查过的内容

- 仓库：STATUS、TODO、ROADMAP、AGENTS、WEHI工作流、迁移规划第九节、TME注释/CM/状态goal、下一阶段研究方案、状态读出验证方法，以及10-30各线最新任务日志尾部。
- WEHI：`${PC_ROOT}`布局、`README.md`、`index/runs_by_line.tsv`（按线计：0线166、1线注释198、1线CM 70、2线78、3线632个运行）、`views/`、`envs/specs/`、`~/envs`、`quota -s`、`sinfo -s`、2026-10-30/31运行的终态字段。
- 旧服务器（SSH只读）：快照后新运行、进程列表、`/data4`占用、`migration_stage`大小、新运行的`acceptance.json`。

## 关键发现（已观察）

- 旧服务器没有分析进程在运行；快照后只新增`20261031_TME_CM_clinical_structure_review_v1`（276 KB，`COMPLETE_READONLY_DESCRIPTIVE_AUDIT`，`n_fits: 25`，`no_clinical_tests: true`）。
- 旧服务器`/data4`已用93%，`migration_stage`占305 GB。
- WEHI只有`pc-repo`、`pc-spatial`两套环境；旧服务器四套分析环境只有清单。home配额不够容纳重建环境。
- TODO（约175 KB）与STATUS中大量“当前”节已过时，但按AGENTS不得删除历史，因此只在顶部加入以新运行计划为准的入口。

## 修改的文件

- 新增`methods/next-stage-run-plan-20261031.md`：现状表、三条分析线完成/未完成对照、P0–P4运行计划、决策点和需用户操作的事项。
- `STATUS.md`、`TODO.md`、`ROADMAP.md`：顶部各加一节指向该计划；原内容未改。
- 本日志。

## 验证

运行`bash scripts/validate_repository.sh`与`git diff --check`，结果见提交前输出。

## 未验证

- T/NK父群细化链和标准A/B/C/研究留出的最终终态没有逐个读取，列为P1第1步。
- scratch清理的具体日期按“30天未访问”推算，不是系统给出的日期。
- 计划中的运行全部未执行。

## 待用户决定

持久空间申请；是否删除旧服务器`migration_stage`与R2桶/令牌；上传本机私有文件；S16主分割方案与SCP2510下载。

## 建议commit message

Add post-migration project overview and next-stage run plan

## 推送时发现远端新提交并合并

首次push被拒：远端已有旧服务器推送的`91289b2`（TME CM首轮结构与接口审查报告及计划）。两边文件不重叠，按AGENTS禁止rebase，做普通merge。据该报告修订运行计划：线①首轮阅读已完成，P1改为临床资格、免疫五轴重新测量与CM×状态关系，并注明旧服务器Codex仍在推进该线，切换计算前需对齐。

## P0执行（用户2026-10-31要求：推进P0重建环境等，删除迁移暂存区；私有文件稍后由用户上传）

详细记录见[run](../../../runs/20261031_WEHI_compute_switch_P0_v1.md)。要点：

- 删除旧服务器`migration_stage`（305 GB）前核对8项VERIFIED、318,117个文件新位置一致；硬链接树删除不影响原件；只读目录需先对暂存区目录加`u+w`。R2桶与令牌未删（未在本次授权中）。
- 增量同步1个新运行，SHA一致。
- 环境：singler_bioc、scanpy、scvi重建完成且精确版本一致（scanpy补正mpmath/sympy）；完整scrna_r因Anaconda默认源在WEHI不可达未能重建，改用conda-forge r-base同构建+CRAN归档精确版本建CoVarNet子集环境。
- 仓库检查658项在`pc-repo`有3项不过，其中2项在重建的scanpy参考环境中通过，1项为测试布局问题。
- 重放中发现迁移缺口：Project_v3 `harmonized_per_gse`（20 GB）及3个原始文件；作为M5经R2补迁（21 GB，VERIFIED，新位置111/111一致），旧服务器暂存随后删除。
- 新增bwrap旧路径沙箱和quota替代脚本；pyUCell重放与旧服务器逐字节/逐值一致；CellTypist路径等价检查通过；CellTypist跨机器全量重放和CoVarNet重放进行中。
- 失败与修正：rsync在慢链路断开（改ssh tar）；`srun`在会话作业内卡住（改sbatch）；一次范围过宽的`scancel --name`被安全检查拦下，未执行；`sbatch --export`按逗号拆分包名导致只装了NMF；`micromamba run`导致stringi configure失败；沙箱重复绑定导致df名称变化。
- P0.4完成：CellTypist v4全量重放64/64文件相同；CoVarNet共同K4在`scrna_r_covarnet2`（r-base h5074ccb_16、ICU 75.1、52个R包精确版本）中重放，表格完全相同或差≤1.1e-14，RDS内模型与网络一致。第一个子集环境因补装开发库时求解器更换r-base构建和ICU而作废。

## P0.4b：scVI、SingleR、MOFA、Harmony/UMAP重放（用户要求“先验证scVI、SingleR、MOFA、Harmony/UMAP等”）

- SingleR、MOFA一致（MOFA最大相对差6.3e-8）。
- HVG：Broadwell节点上30/4000个HVG不同；AVX-512节点上HVG/PCA完全相同。解释为NumPy在AVX-512上的不稳定排序对并列值的处理不同。Harmony嵌入相对差中位2.5e-5。
- Leiden/UMAP：边集合相同、权重差1.4e-6；Leiden 1.0/2.0相同，主分辨率0.5 ARI 0.98（0.86%换群）；UMAP不同。Skylake重放用于判断是否为CPU微架构所致。
- scVI：同型号GPU逐位可复现；V100与WEHI各GPU之间、以及WEHI GPU型号之间15近邻重合约0.77–0.82。按事先写入脚本的规则不通过；差异来自GPU架构。
- 失败与修正：重放未复制`tests/`导致Harmony脚本自检失败；两个scVI副本共用非阻塞锁，改为串行依赖；MOFA的hdf5/npz最初只按字节比较，比较工具已扩展。
- Skylake重放与Ice Lake重放逐字节相同，Leiden/UMAP与旧服务器的差异不是CPU架构所致，原因未查明。
