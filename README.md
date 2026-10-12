# PanCancer 项目进展与工作日志

## 2026-10-12T14:15:47+08:00 全量审计编号解析需修复，合并已按门槛停止

792文库逐行检查均PASS，但GSE256326来源级文库编号对比失败。[全23880行诊断](inventories/runtime/2026-10-12/native_QC_identifier_read_diagnosis_v1.json)确认CSV自动int丢前导零导致5958表面差异；字符串读取全部匹配。原数据/模型正确保留，正在准备v4完整重验，旧失败报告封存。最终Zarr尚未写入，不能称全量已验收；新修复代码冻结后继续，本阶段新下载暂停。



## 2026-10-12T14:06:46.979472+08:00 完整控制空间关联真实运行中

2026-10-12T06:04:30.803833+00:00 UTC已启动driver144523/child144526，当前两个实际进程存活，最新finished_group GRP016 models 3166。全26994记录、原严格状态及必需controls保留，2000完整块置换/bootstrap与全家族BH计算中；独立audit尚未完成，不发布筛选结果。[启动证据](inventories/analysis/2026-10-12/CM_LR_full_control_association_v1/execution_start_v1.json)。下一检查同句柄真实终态、全量数值/重采样审计及来源/癌种一致性，不因旧等待状态重启。


## 2026-10-12T13:58:47.403118+08:00 全9源792文库方法完成，正在全量独立审计

[792检查点](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint16.md)：0失败、3,180,973保留细胞、108,635额外预测双细胞移除。独立audit已开始全行复算，全部PASS后进入稳定ID Zarr合并与再次全行验收。此时尚未生成最终合集；来源最终AC/污染可辨识性限制保留，原卷最终/备用卷临时，新下载暂停。



## 2026-10-12T13:29:22.289964+08:00 公共前处理约95%，最终证据完整性预检通过

[进度](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint15.md)：752/792、0失败、3,048,609细胞；8源已结束，仅SPECTRUM40文库继续。[证据预检](inventories/runtime/2026-10-12/durable_evidence_preflight_v1.json)36/39约42MB，只欠当前未生成的最后总结/完整审计/终态，500GB留空可满足。真实全量audit和Zarr仍待；正式来源AC待，最终原卷、备用卷仅临时，新下载暂停。



## 2026-10-12T13:20:43.440767+08:00 原CM完整LR提取与独立验收完成

332文库/8GSE/10情境，完整314432候选覆盖及6996788原行、381170目标观测、40321280值全复算，去两端活动2180KKT通过。5377案例达数据覆盖门槛，去重3118接口、598多来源有覆盖；尚未拟合关联，不是显著/验证轴。[结果](inventories/analysis/2026-10-12/CM_edge_LR_extraction_v2/RESULTS.md)。COL1A1–DDR1当前明确端点在乳腺/CRC/胃癌可测；CCL5–CCR1在CRC/GBM可测，下一全案例设计与丰度调整。CRC占约77%案例，患者/身份/分选/表达过滤限制保持。


## 2026-10-12T13:16:02.951797+08:00 公共前处理8来源结束，仅卵巢癌继续

[进度](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint14.md)：734/792、0失败、2,991,089保留细胞、102,877额外预测双细胞移除。乳腺和LBCL已完成逐文库步骤；SPECTRUM仍运行。稳定ID Zarr组装v2等待完整9源审计；最终写原卷、临时数据留备用卷，尚未完成最终矩阵。正式来源AC待，新下载暂停。



## 2026-10-12T12:39:19.265446+08:00 公共前处理超过600文库

[进度](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint13.md)：606/792、0失败、2,530,660保留细胞、90,134额外预测双细胞移除。稳定基因ID最终组装v2等待完整9源审计；尚未生成最终矩阵。原卷最终/备用卷临时的角色已落实，旧权威库保持，新下载暂停，正式来源AC待。



## 2026-10-12T12:28:42+08:00 稳定ID最终组装v2已启动

[v2启动](inventories/runtime/2026-10-12/combined_zarr_assembly_public_start_v2.json)冻结四组件及Git d8a53ee，560/792文库、0失败、2380251保留细胞。现有3源继续前处理；全9源审计PASS后进入稳定ID合并与全行独立验收。最终原卷，备用卷仅加工；旧等待器历史保留。尚未生成最终矩阵，新下载暂停，正式来源AC待。



## 2026-10-12T12:26:24+08:00 最终组装修正为稳定基因ID

真实基因预检发现名称别名会拆分同基因：72554名称将按可信稳定ID整为49985特征，400歧义/1244仅名称特征保留标记。人工别名求和和完整读回PASS；原预处理继续。旧组装等待器在写矩阵前停止并保留，下一[独立v2组装](runs/20261012_public_preprocessed_zarr_v2.md)待代码冻结启动，不能称当前已生成最终文件。最终仍原卷，新下载暂停、正式来源AC待。



## 2026-10-12T12:23:07.629760+08:00 原9CM全来源覆盖验收完成

60唯一细胞对/120模块边全部检查，11220分母行和240有向映射独立复算通过。当前10个不同情境中的34条细胞对×情境记录、15细胞对、8来源达20文库/10分析块门槛，4对在多个来源及癌种标签有覆盖；全部15对映射涵盖9CM。[结果](inventories/analysis/2026-10-12/CM_edge_coverage_v1/RESULTS.md)。这是筛查机会而非通讯效应，样本替代单位和分选限制保留。下一完整LIANA全LR提取范围已冻结，尚未派发；原数据/身份/W/H不变。


## 2026-10-12T12:18:34+08:00 最终Zarr组装已启动，等待全量前处理审计

[启动记录](inventories/runtime/2026-10-12/combined_zarr_assembly_public_start_v1.json)核实原卷新run/冻结代码/PID；既有前处理[528/792文库、0失败、2236179细胞](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint12.md)，6来源完成、3来源运行。全9源独立审计通过后执行流式合并及全行验收；目前仍待，不能把启动当作最终完成。正式来源AC和用途限制继续保留，新下载暂停。



## 2026-10-12T12:16:47.369510+08:00 统一Zarr容量可满足，组装代码准备

真实498文库压缩预检显示最终矩阵及写入余量约75.29GB，原卷有效配额593.64GB，可保留500GB空闲。[方法](methods/combined-public-preprocessed-zarr-v1.md)采用原计数、稀疏校正差值和重算X，保留基因测量掩码与来源限制。合成读回已PASS；真实9源全量审计及合并尚待，最终写原Project_v5，不更新旧权威库。[run](runs/20261012_public_preprocessed_zarr_v1.md)。



## 2026-10-12T12:10:52.892668+08:00 候选原生端点核查完成

肺癌ADAM12–ITGB1关联保留，但发送端仅1.05%RNA阳性、接收ITGA9仅1细胞且与ITGB1共检出0。多数阳性发送细胞保留T谱系marker，不自动改名或判污染；原生计数/已审全基因汇总匹配。[核查与后续](inventories/analysis/2026-10-12/ADAM12_native_endpoint_probe_v1/RESULTS.md)。本轴为来源限定表达关联，机制证据不足；没有新增验证轴，下一原9CM全谱系自动筛查与外部可比来源审阅尚未运行。


## 2026-10-12T12:07:15.946820+08:00 先验靶表达关联完成，候选尚待端点验证

全部620情境×3模型已独立验收；M1 241可估计、65诊断q阳性、10稳定区间，仅肺癌ADAM12–ITGB1→IL7R_36同时满足两项。M2 35可估计、0稳定区间，不能据此称通讯不存在。[结果](inventories/analysis/2026-10-12/target_association_v1/RESULTS.md)。当前无新增独立验证轴，CM分支登记无活跃任务。下一固定原生端点描述核查已预声明，未运行；原9CM全谱系通讯及外部验证继续按主线推进。


## 2026-10-12T11:50:31.401241+08:00 用户确认临时存储，前处理已恢复

备用卷仅作临时加工，最终统一预处理Zarr写原项目卷的新独立版本。03:45:11 UTC已按真实PID/命令核实后恢复原383文库检查点；[恢复授权记录](inventories/runtime/2026-10-12/ambient_temporary_volume_resume_v1.json)取代下方暂停历史。将分块合并稀疏counts/污染估计/明确选层的表达，保留来源、技术和缺失基因掩码，并全量独立验收。尚未合并成功；来源最终AC与污染复核限制继续保留。


## 2026-10-12T11:37:06.841860+08:00 公共预处理暂停，待输出卷确认

用户质疑备用卷项目使用权限和开销，本轮科学/记录进程已可恢复暂停；既有结果保留。实际新run **136.37GiB，约69.3%为R交换文本**；500GiB只是未分配的留空预算。原计划为下一批峰值+500GB，本轮粗估扩大成全9源+500GiB，需要纠正并按必要产物/逐批重新核算。[权限/空间说明](inventories/runtime/2026-10-12/STORAGE_PERMISSION_REVIEW_v1.md)。暂停时383/792、0失败，1,751,187细胞；未全量完成或整合。此暂停记录优先于下面的运行历史。


## 2026-10-12T11:29:47.169658+08:00 公共前处理仍后台运行

[11:28快照](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint11.md)：367/792文库完成、0失败，保留1,691,585细胞、额外移除66,451个预测双细胞。已审计3完整来源92文库579,361细胞；完整独立终态尚未到达。[终态记录组件v2](inventories/runtime/2026-10-12/preintegration_finalizer_public_start_v2.json)已真实启动、全部代码SHA及主进程核实，完成后写可读状态与准备索引。作者原有污染处理、脑谱系面板、CAR-T状态筛除及来源排重仍分别受控；未做实际整合，新下载仍暂停。


## 2026-10-12T11:26:56.773706+08:00 公共来源预处理继续后台执行

[11:20处理快照](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint10.md)：323/792文库完成、0失败，保留1,517,029细胞，额外移除60,087个预测双细胞。3完整来源92文库579,361细胞已加强独立验收；[准备索引v2](inventories/runtime/2026-10-12/preintegration_preparation_checkpoint_v2.json)保留所有整合门槛。作者原有污染处理及脑谱系面板限制已落实；原始计数保留。当前不是全量终态，未开始实际整合。


## 2026-10-12T11:15:43.446858+08:00 已下载来源前处理与整合准备

[处理快照](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint9.md)：11:11时274/792文库、0失败，保留1,340,180细胞。3个完整来源92文库579,361细胞已加强独立验收，并生成[准备清单](inventories/runtime/2026-10-12/preintegration_preparation_checkpoint_v1.json)。LBCL和乳腺作者已记录污染处理，追加估计层须作为敏感性审阅；整合[分支/历史合同](methods/public-source-integration-preparation-v1.md)明确限制。全量仍在后台，来源最终AC为0，未实际整合。


## 2026-10-12T10:54:17+08:00 公共来源前处理：3队列独立验收通过

[最新处理快照](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint6.md)：10:53时182/792文库完成、0失败，保留944,249细胞。已完整结束的CAR-T产品、GSE292046和GSE294109共92文库、579,361细胞，逐保留细胞的原始计数、ID、过滤映射、校正及标准化[独立验收通过](inventories/runtime/2026-10-12/ambient_checkpoint_audit_cached_v2.json)。这不是全量终态；污染不可辨识及过度校正复核仍保留，来源最终AC为0，未全局整合。


## 2026-10-12T10:38:45.707906+08:00 已下载公共源的双细胞与RNA污染处理

[处理快照](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint2.md)：备用卷恢复成功，10:27时44/792文库、0失败；完整验收未完成。原计数、校正估计和复核限制分开。[方法](methods/doublet-ambient-preintegration-v1.md)。新公共下载与实际全局整合仍未恢复。

## 2026-10-12T10:38:40.719500+08:00 全基因输入里程碑已独立验收

3源285,360细胞/164文库、1121原编号单位、139,558,612基因汇总值及全部336/15/1008先验映射独立复算通过；已登记v17/0活跃/0冲突。[验收证据](inventories/analysis/2026-10-12/fullgene_readiness_v17/public_checkpoint_v17.json)。321/336有配体靶先验，仅输入就绪，不代表受体特异响应/通讯轴验证。45源header/180家族覆盖独立审计通过。194旧患者状态差异已追溯至官方侧表；新[非破坏性证据覆盖表](inventories/analysis/2026-10-12/patient_block_provenance_overlay_plan_v1.json)已冻结待执行。

## 2026-10-12T10:11:30.425815+08:00 配额扩充后恢复既有分析

用户明确要求继续；实测组硬限2048G、已用约1495G。新[离线资源v4](inventories/analysis/2026-10-12/fullgene_offline_resource_recovery_v4_spec.json)复用本机已校验target，不重复公网下载；生物范围/producer/auditor不变。当前尚未实际传输/启动，先冻结代码及普通push，再核验输入和真实作业。[冻结](inventories/analysis/2026-10-12/fullgene_receiver_target_code_freeze_v4.json)；下方昨日暂停是历史节点。

## 当前计划与下载暂停：2026-10-11T23:50:48.128458+08:00

用户要求今日只分析后续任务，暂不下载；已形成[研究主线与主线内支线](ROADMAP.md)、[任务板](TODO.md)。原9CM生态与细胞对响应双入口，按完整基因/可比单位→响应→空间/临床证据推进。只读检查根GPU目录已不存在，Project_v5对应目录、Zarr及当前登记证据路径保留；未恢复下载或启动分析。[取证及限制](inventories/analysis/2026-10-11/planning_pause_checkpoint_v1.json)。后续执行以此最新计划为准，早期资源恢复和下载指令保留为历史。

## 最新QC驗收

2026-10-11 21:18 UTC+8：[四分片全量技術QC驗收](inventories/runtime/2026-10-11/QC_ACCEPTANCE.md)。993,312細胞全部對原始計數/ID/指标/文庫映射核查通過，頭頸癌7,463細胞患者映射已版本化修正。可來源內前處理；正式新增來源驗收仍待別名/基線文庫審查，產品/富集snRNA的用途限制保留。


## 最新公共下载检查

2026-10-11 18:56—19:02 UTC+8：[下载/QC/验收及下一步报告](inventories/runtime/2026-10-11/DOWNLOAD_READINESS.md)。5,341,204源细胞；4份QC分片993,312细胞；正式新增验收0。下载进程已退出，第二批8/14文件校验，未完成文件不计矩阵到库。逐来源QC/验收可继续，全量整合未就绪。下方11:22入口和动态描述保留为历史快照。


## 最新CM探索节点与时间线

2026-10-11T21:39:13.099905+08:00：[资源断流终态及45来源原生可用性](runs/20261011_fullgene_receiver_target_readiness_v1.md)。target分块反复截断已停止，准备aria2恢复；生物scope和输入验收门槛保留，未产生新统计结论。

2026-10-11T21:06:33.599638+08:00：[完整基因资源失败与恢复轨迹](runs/20261011_fullgene_receiver_target_readiness_v1.md)。本机3先验校验通过；服务器仍有部分文件，真实生产未开始。事件/取证时间及SHA已补录，后续按关键节点普通push核验。

2026-10-11 20:29 UTC+8：[70接收程序模型全量验收](runs/20261011_spatial_receiver_programme_exposure_v1.md)。22可估计、0名义p05/0主BH，v16/0冲突；当前通用程序未支持配体邻居效应，继续受体靶/TF及全基因队列审查。

2026-10-11 19:54 UTC+8：[固定45空间邻近全量验收](runs/20261011_spatial_conditional_neighborhood_v1.md)。29可估计/16不可测、完整45主BH无显著，7,074,000置换独立复算、v15/0冲突；未把计算完成当通讯轴验证。

2026-10-11 19:27:55 UTC+8：[全部65空间端点预检完成](runs/20261011_spatial_contact_native_preflight_v3.md)。正确作者181,056标签全匹配，6口径复算通过；4种LR在作者粗谱系口径有描述性覆盖，尚待邻近检验。失败和方向错误保留、显式更正。

2026-10-11 18:52:59 UTC+8：[空间输入修复完成](runs/20261011_spatial_native_input_repair_v1.md)；保留失败及旧对象，新对象原CSV全矩阵独立核对通过。下一阶段全部13×5候选端点预检；通讯/CM机制仍未验证。

2026-10-11 18:36:08 UTC+8：[空间输入检查与CSV索引故障](logs/incidents/2026-10-11_spatial-csv-index-column.md)。5份文件只读汇总完成；TMA2暂停归一化/空间推断，先修复并审计；患者映射/背景/邻近尚待验证。详见[持续任务日志](logs/codex/2026/2026-10-11_milestone-log-and-cloud-sync.md)。

2026-10-11 18:20:36 UTC+8：[关键节点证据](inventories/analysis/2026-10-11/README.md)、[完整时间线](inventories/analysis/2026-10-11/task_timeline_v1.tsv)、[正式run](runs/20261011_family_candidate_validation_v1.md)。336条参考支持记录完成检查，自动15条后续验证记录；身份内部程序检查完成，空间接触和机制尚待验证。关键/探索节点的记录与普通push规则见[RUNTIME_SYNC](docs/RUNTIME_SYNC.md)。

## 2026-10-11 当前入口：服务器运行与Codex云端开发

- **先看下载与Project_v3动态**：[本次实测报告](inventories/runtime/2026-10-11/REPORT.md)，采集于2026-10-11 11:22:53（UTC+8）。GitHub保存快照，不是实时作业面板。
- **当前状态/下一步**：[STATUS](STATUS.md)、[任务板](TODO.md)、[阶段路线](ROADMAP.md)。
- **开始云端开发**：[Codex Cloud配置与工作流](docs/CODEX_CLOUD_DEVELOPMENT.md)，[可直接提交的首个任务](prompts/cloud_next_task.md)。
- **检查代码**：Python ≥3.10环境执行 `bash scripts/cloud_check.sh`；状态采集与同步见[操作说明](docs/RUNTIME_SYNC.md)。
- 新公共源数据5,341,204细胞；608,684细胞已有QC分片；最终新增验收尚未登记。Project_v3表达标签及新CoVarNet/LIANA计算已完成，最终科学细身份仍待GSE299623审阅。
- 目标仓库已核实为**public**。2026-10-11用户明确授权向此仓库提交与普通push；仅上传已审查的代码、汇总与脱敏记录。下文此前版本入口保留为历史依据，不覆盖本次状态。

本仓库是 PanCancer 项目的轻量级、可审计记录系统，当前以 `Project_v3` 的泛癌单细胞转录组主线为主要范围。它记录当前状态、任务历史、方法、正式运行、技术决策及证据索引，是项目进展的权威入口，但不是分析数据仓库。

## 2026-09-22 投稿冻结入口

- Nature 投稿准备与注释上游锁：[`manuscript_upstream_freeze_nature_20260922`](releases/manuscript_upstream_freeze_nature_20260922/README.md)。
- V8 最新谱系审计限定：[`annotation_v8_lineage_audit_20260917`](releases/annotation_v8_lineage_audit_20260917/README.md)。
- V9 最新候选报告与开放复核表：[`annotation_v9_candidate_20260919`](releases/annotation_v9_candidate_20260919/README.md)。
- 全局投稿参考仍是 V7；V8 仅为 lineage taxonomy，V9 为 `FROZEN_CANDIDATE`，均未晋升为全局参考。
- 当前 V5-derived 通讯/空间/临床图不能无说明地称为 V7/V9-derived。

## 2026-09-12 权威入口

- 当前细胞身份/CNV 恶性参考：[V7 冻结 release](releases/annotation_v7_20260905/README.md)。
- [V8 development software candidate](releases/annotation_v8_development_20260912/README.md)
  是已索引的 `NOT_FROZEN` 全量元数据软件候选；它不替代 V7，也不授权正式 downstream
  变更或生物学完成声明。
- 当前论文通讯证据与 Prism 成图：[20260906 evidence chain](releases/communication_evidence_chain_20260906/README.md)。
- Git 维护：[维护工作流](docs/GIT_MAINTENANCE.md)。
- Prism 与 supplementary：[交付指南](docs/PRISM_SUPPLEMENTARY_GUIDE.md)。
- 版本边界：[版本政策](docs/VERSION_POLICY.md)。


## 不存储的内容

本仓库不保存原始测序数据、患者身份信息、受控访问数据、H5AD/RDS/矩阵对象、完整作业日志、缓存、压缩归档或大型分析产物。真实文件留在分析工作区，本仓库只记录稳定标识符、逻辑路径、校验信息和验证结论。详见 [数据政策](docs/DATA_POLICY.md)。

## 快速查看当前状态

1. 阅读 [STATUS.md](STATUS.md) 获取当前阶段、阻塞和下一步。
2. 阅读 [ROADMAP.md](ROADMAP.md) 获取阶段与里程碑。
3. 阅读 [TODO.md](TODO.md) 获取可执行任务。
4. 通过 [运行索引](inventories/runs.tsv) 找到正式运行记录及外部证据。
5. 通过 [术语表](docs/GLOSSARY.md) 和 [命名规范](docs/NAMING_CONVENTIONS.md) 理解标识符。

## 日常使用

```bash
# 创建当天日志
bash scripts/new_daily_log.sh

# 创建 Codex 任务日志
bash scripts/new_codex_log.sh "review phase 07 annotation"

# 创建正式运行记录
bash scripts/new_run_record.sh "phase07 annotation"

# 输出轻量级项目快照；结果目录参数可选
bash scripts/project_snapshot.sh "${RESULTS_ROOT}"

# 只检查、不修改
bash scripts/validate_repository.sh
```

新日志创建后必须补充实际内容，不能保留未替换占位符。正式运行还应同步更新 `inventories/runs.tsv`；新数据集或重要产物分别更新 `inventories/datasets.tsv` 和 `inventories/artifacts.tsv`。

## Git 工作流

```bash
git status --short
git diff --stat
git diff
bash scripts/validate_repository.sh
git add <reviewed-files>
git commit -m "type: concise description"
git push
```

本仓库当前获得了用户对常规 `commit` 和向已确认 private 远程执行普通 `push` 的长期授权；尚未配置远程时不得猜测地址。强制推送、改写历史、删除文件和破坏性 Git 操作始终禁止，详见 [AGENTS.md](AGENTS.md)。

## 文档职责

| 权威位置 | 只负责 |
| --- | --- |
| `STATUS.md` | 当前真实状态 |
| `ROADMAP.md` | 阶段、里程碑、依赖和完成标准 |
| `TODO.md` | 可执行任务 |
| `CHANGELOG.md` | 重要变化 |
| `logs/daily/` | 当天发生的事情 |
| `logs/codex/` | Codex 每次任务的检查、修改和验证 |
| `runs/` | 正式运行的输入、参数、命令、输出和结果 |
| `methods/` | 当前认可的方法 |
| `decisions/` | 重要方案选择及理由 |
| `inventories/` | 数据集、运行和产物的结构化索引 |

相同事实只在一个权威位置详细维护，其他文档通过相对链接引用。完整关系见 [目录结构](docs/DIRECTORY_STRUCTURE.md) 和 [维护工作流](docs/WORKFLOW.md)。
