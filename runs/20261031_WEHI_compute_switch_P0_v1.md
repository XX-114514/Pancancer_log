# WEHI计算切换准备（P0）

依据[运行计划](../methods/next-stage-run-plan-20261031.md)P0。2026-10-31 UTC执行；时间均为UTC，另注UTC+8。

## P0.5 旧服务器迁移暂存区删除（用户2026-10-31明确授权）

- 删除前核对：WEHI上8项回执均为`VERIFIED`；新位置复核318,117个文件全部一致（`migration/verify_relocated.tsv`，n_bad与n_unmapped均为0）。
- 暂存区内容：8项的`chunks/`共220个分块及清单；`m4a_inputs_tree`、`m4b_runs_tree`、`m4c_env_tree`为指向原件的硬链接（链接数2），删除只去掉链接。
- 首次`rm -rf`遇到只读目录（冻结对象保留了原目录权限）报Permission denied；只对暂存区内的目录加`u+w`（目录不是硬链接，不影响原件；文件权限未改）后删除完成。
- 删除后：`migration_stage`不存在；抽查原件`qc_GSE294109.zarr/X/indptr`仍为只读、内容在；`/data4`由41T已用/3.5T可用（93%）变为40T/3.7T（92%）。
- 未做：R2桶与令牌删除（不在本次授权中，仍待）。

## P0.1 增量同步

旧服务器快照（2026-10-30 23:04 UTC）之后只新增`20261031_TME_CM_clinical_structure_review_v1`（2个文件，276 KB）。rsync与第一次`ssh tar`因链路断开失败，第二次`ssh tar czf`（40 KB）成功。两文件SHA256与旧服务器一致，回执`migration/incremental/20261031_sync1.tsv`；视图和索引重建，`1_TME_CM`由70变为71个运行。

## P0.3 私有wheel与仓库检查

`ATLAS_SYMPHONY_WHEEL`、`ATLAS_SYMPHONY_SHA256`写入本地配置；wheel SHA256 `6024b7d9…373a`与[依赖核查记录](../inventories/analysis/2026-10-23/full_context_annotation_v1/technical_reference_dependency_inspection_v1.json)一致。

`pc-repo`环境、4 CPU作业30608295（repo_commit `6bc893d`，dirty 2为本次新增的未提交脚本）：658项测试，1 failure、2 errors、2 skipped。

| 测试 | 原因（已观察） | 判断 |
|---|---|---|
| `test_single_source_symphony_replay_v1`（分块17行与整块结果比较） | 最大绝对差3.0×10⁻⁸，测试容差1×10⁻⁹ | 浮点累加顺序差异；需在重建的scanpy参考环境中复测 |
| `test_standard_annotation_order_primitives_v1` | `pc-repo`缺`skmisc`；旧scanpy环境有scikit-misc 0.5.2 | 测试需要参考环境，`pc-repo`不擅自加包 |
| `test_standard_annotation_reference_v2` | 按`<仓库根>/config/source_holdout_spec_v1.json`读取，该文件只在运行目录`runs/*/config/`中 | 仓库测试布局问题，与迁移无关 |

此前1项因wheel缺失的报错已消除（该测试现可加载wheel，失败原因变为上表第一行）。

参考环境复测（作业30608315，重建的scanpy环境，4线程）：`test_standard_annotation_order_primitives_v1`、`test_single_source_symphony_replay_v1`、`test_frozen_symphony_reference_v1`、`test_annotation_source_holdout_v1`共14项全部通过，包括上表前两项。因此前两项是`pc-repo`与参考运行时的差别，不是迁移数据或代码问题；第三项为测试布局问题，未修改。

## P0.2 环境重建

工具：[wehi_rebuild_env.sh](../scripts/wehi_migration/wehi_rebuild_env.sh)、[rebuild_env.sbatch](../scripts/wehi_migration/jobs/rebuild_env.sbatch)。conda explicit（不求解）→ 补装环境中缺少的PyPI精确版本（`--no-deps`）→ scrna_r复制迁移来的非CRAN R库 → 输出pip与R包差异表。环境装在`${PC_ROOT}/envs/conda/`。作业30608290（array 0–3）。

| 环境 | 结果 |
|---|---|
| singler_bioc | 47个R包版本与旧环境逐一相同（差异表为空） |
| scvi | conda explicit完成，补装29个PyPI精确版本；pip差异表为空 |
| scanpy | 补装59个PyPI精确版本；R包差异为空。首轮漏装mpmath 1.3.0、sympy 1.13.1：旧环境用pip覆盖了conda的1.4.0/1.14.0，脚本只按包名判断。已补装并改脚本为按`name==version`判断；之后PyPI精确版本与旧环境完全一致 |
| scrna_r（完整） | 失败：842个包中156个来自`repo.anaconda.com/pkgs/main`、186个来自`pkgs/r`，WEHI连不上该地址（conda-forge可连），也没有本地镜像 |
| scrna_r_covarnet2（CoVarNet子集） | conda-forge `r-base 4.3.3 h5074ccb_16`（与旧环境同一构建）、ICU 75.1 + 工具链与libxml2/zlib/liblzma开发库，一次求解建成；再按旧环境记录的版本从CRAN归档安装NMF、dplyr、reshape2、igraph、psych、jsonlite及其依赖闭包：52个包全部为旧版本 |
| scrna_r_covarnet（作废） | 第一个子集环境。后续往里补装开发库时，求解器把r-base换成`h05019bb_23`、ICU换成78.3，已编译的stringi因找不到`libicui18n.so.75`无法加载。不再使用，保留目录待清理 |

旧scrna_r中多数R包实际来自CRAN源码安装而非conda（R包清单中457个标为CRAN来源；例如conda给的igraph是1.5.1，实际加载的是CRAN的2.2.0），所以“conda-forge r-base + 按版本装CRAN包”比还原Anaconda默认源的包更接近实际运行环境。子集构建的失败与修正：首轮stringi的configure报“C++ compiler cannot create executables”与“5: Bad file descriptor”，原因是`micromamba run`不提供configure需要的文件描述符5，编译器激活脚本也未生效；改为`micromamba activate`后通过。版本比较原先按字符串，`0.5-1`与`0.5.1`被误判为不一致，已改为`package_version`比较。另外三处失败：`sbatch --export`按逗号拆分参数，`PKGS=NMF,dplyr,...`实际只传入NMF（改为冒号分隔）；当前CRAN的Matrix要求更新的R，R 4.3的包索引里查不到它的依赖，lattice没排在前面（安装脚本改为最多3轮重试）；igraph的configure用到系统libxml2头文件但env缺开发库（建环境时加入）。还有一次是作业运行中修改了正在被`Rscript`逐条读取的安装脚本，作业末尾出现语法错误；安装报告在此之前已写完，以后作业运行期间不改其使用的脚本。

工具：[wehi_install_R_exact.R](../scripts/wehi_migration/wehi_install_R_exact.R)、[build_R_subset_env.sbatch](../scripts/wehi_migration/jobs/build_R_subset_env.sbatch)。完整scrna_r（CellChat、copykat、scATOMIC等）未重建，待决定。

## M5 补充迁移：Project_v3原生h5ad

pyUCell首次重放在读基因轴时失败：旧版来源的`feature_lookup`读取`Project_v3/annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse/<GSE>_symbol_counts.h5ad`，该目录（54个h5ad及完成标记，20 GB）不在M0–M4中；另有诊断脚本从这些文件读原始行做核对，所以只给基因索引不够。随后扫描992个近期运行的config/provenance/outputs JSON中28,350处旧绝对路径（去掉JSON转义尾部后）：真正缺失的还有`Project_v3/GSE166555/GSE166555_rawdata.h5ad`和GSE189889的2个文件（共约1.1 GB）；其余为可重新下载的文献补充材料、按计划不迁移的Project_v3旧CoVarNet运行、以及slurm日志和锁文件。

M5共113个文件、21 GB，经同一R2通道（单个上传进程，桶上限9G）：旧服务器02:38 UTC建暂存（硬链接）、03:44上传完成；WEHI 03:50 `VERIFIED`，桶已清空。随后移入`${PC_ROOT}`的旧布局位置（记入`migration/relocation.tsv`），在新位置按旧服务器SHA清单复核111个数据文件全部一致（作业30608535）。旧服务器暂存随后删除，原件完好。下载日志中一条`Failed to read versioning status ... AccessDenied`来自令牌只限本桶、无权读版本状态，不影响传输。

## P0.4 数值等价重放

旧服务器配置里层层嵌套绝对路径，逐个改写会改动原件。新增[wehi_oldpath_sandbox.sh](../scripts/wehi_migration/wehi_oldpath_sandbox.sh)：用bwrap把`PATH_MAP.tsv`中的旧前缀只读绑定到WEHI位置，旧环境路径绑定到重建环境，只有重放输出目录可写。测试确认旧路径可见、写入被拒。

沙箱补充：
- `quota`替代脚本（[shims/quota](../scripts/wehi_migration/shims/quota)）。旧阶段的存储门槛调用`quota -g -v -w USER002`（旧服务器组配额）；WEHI scratch按用户计，替代脚本改为查询当前用户在WEHI的真实配额，列格式相同，数值不伪造。沙箱内用户名经winbind解析失败，因此按数字UID查询。
- 首版沙箱把输出目录在原路径上再绑定一次，结果`df`把它当成独立挂载点，报出带子路径的文件系统名，门槛因匹配不到配额行而报“Ambiguous quota”（CellTypist第二次重放）。`/vast`本身已可写绑定，去掉这次重复绑定后`df`正常，输出目录在新旧两条路径上可写，其他旧运行仍只读。
- 不能在本会话中用`srun`：会话本身运行在一个2 CPU的Slurm作业里，`srun`会在其中开step并卡住。一律用`sbatch`。

| 对象 | 原运行 | 方式 | 结果 |
|---|---|---|---|
| pyUCell | `20261030_predefined_state_UCell_smoke_v2`（2来源256细胞，126程序，2线程） | 原config/scripts/resources/vendor不改，沙箱内scanpy环境 | **一致**。首次（作业30608313）因缺M5的h5ad失败；补齐后（30608540，exit 0）16个数据文件全部一致：两个来源的UCell分数`.npy`与旧服务器逐字节相同，各单位均值/分位数和观察表逐值相同；3个JSON只有`wall_seconds`不同 |
| CellTypist（路径等价） | `20261030_TME_CellTypist_unicode_metadata_equivalence_v1`（2来源128探针细胞，6模型） | 同上，8线程 | **通过**（30608558，第三次）。验收JSON除各探针耗时外完全相同：6个模型的匹配特征数与旧服务器相同，两条推理路径的decision/probability逐位一致。注意：此检查比较的是同一次运行内的两条路径，不是与旧服务器预测值比较 |
| CellTypist（跨机器预测值） | `20261030_TME_official_CellTypist_cluster_crosscheck_v4`（续跑v3：7个来源沿用v3结果，其余重新计算；拟合细胞2,946,498个，6模型，原40分钟） | 同上，8线程120 GB | **一致**（30608571，39分钟，exit 0）：64个输出文件全部与旧服务器相同，包括各模型逐细胞最佳标签、最大概率、阳性数、来源覆盖与群投票 |
| CoVarNet共同K4 | `20261030_TME_CM_TNK_Mye_commonK4_CT_identity_v1`（67节点×628观察，43来源，原耗时约12分钟，2线程） | 同一R驱动与冻结作者代码（SHA `917c60c4…`一致），原冻结请求只改输入/输出路径；输入SHA对照冻结接受记录；`scrna_r_covarnet2`环境 | **一致**（30608656，7.5分钟，exit 0）：24个文件完全相同（rank调查、权重、S、两个种子的作者口径与原始BH口径全部边/节点表）；W、A、H、加权H共8个表最大差1.1×10⁻¹⁴；RDS中NMF模型basis/coef最大差4.4×10⁻¹⁵，网络对象的全局/逐模块/过滤/注释部分`all.equal`相同，只有内嵌的原始NMF拟合对象（含耗时与会话信息）不同。`numeric_readback_acceptance_v1.json`由原Python执行器生成，重放直接调用R驱动，故没有此文件。前三次失败（缺包、stringi无法加载）均为环境问题，见上文 |

### 结论

在WEHI上，用重建的scanpy环境和CoVarNet子集R环境，经旧路径沙箱运行原配置，三类计算（pyUCell评分、CellTypist全量预测、CoVarNet nsNMF与作者网络）的输出与旧服务器逐字节或在机器精度内一致。P0.4的切换门槛已满足。限制：只覆盖这三类代码路径；scVI（GPU）、SingleR、MOFA、Harmony/UMAP以及依赖完整scrna_r的CellChat、copykat、scATOMIC等未重放；完整scrna_r尚未重建。

比较工具[wehi_compare_outputs.py](../scripts/wehi_migration/wehi_compare_outputs.py)：TSV与TSV.GZ逐值（rtol 1e-9、atol 1e-12，报告最大差；gzip头含时间戳），JSON逐键并列出不同的键，其他文件按SHA；RDS含R会话信息，只报告差异不判定。
