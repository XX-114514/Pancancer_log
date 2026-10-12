# 当前项目状态

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


## 2026-10-12T11:15:24.098787+08:00 336候选原编号响应预检完成

全部336回映248不重复来源/情境/方向/LR接口，展开620原接收编号记录；325通过计数/变异/九列诊断设计，153文库或块不足、119来源未捕获节点、23无靶先验。独立逐细胞重放496个设计、全部620记录/336映射通过；两阶段真实exit0，v17登记活跃0/状态冲突0。[审计](inventories/analysis/2026-10-12/response_preflight_v1/independent_audit.json)。325记录对应237接口/321原程序记录，不是325独立轴或显著模型。 输入13目录保留，但原336候选实际仅T/B三家族，展开9个接收编号；4个髓系编号未在该候选框架受检，不当阴性，本轮不代表9CM全谱系通讯。

完整推断模型尚未实现：需major比例、发送端原编号混合和共同程序控制、固定响应结局/检验家族及patient block方法。已核12来源×家族count一致，必须用marker_supported_main缓存，provisional不能替代。当前没有新的系数/p值、状态或验证轴。原13编号/172节点/9CM及患者块保留；来源/癌种限定关系继续报告。

## 2026-10-12T10:54:17+08:00 公共来源前处理：3队列独立验收通过

[最新处理快照](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint6.md)：10:53时182/792文库完成、0失败，保留944,249细胞。已完整结束的CAR-T产品、GSE292046和GSE294109共92文库、579,361细胞，逐保留细胞的原始计数、ID、过滤映射、校正及标准化[独立验收通过](inventories/runtime/2026-10-12/ambient_checkpoint_audit_cached_v2.json)。这不是全量终态；污染不可辨识及过度校正复核仍保留，来源最终AC为0，未全局整合。


## 2026-10-12T10:48:39.301099+08:00 患者块证据覆盖表已验收

全1290行原字段/顺序、官方映射与原患者块经独立审计PASS；仅新增有效证据字段，194旧状态显式追溯至已审父侧表。402文库/175来源内患者块、888样本替代文库；原分组和canonical/172/9CM不改。[验收](inventories/analysis/2026-10-12/patient_block_provenance_overlay_v1/independent_audit.json)。全基因v17准备已完成，当前CM分支无登记活跃作业，没有新关联检验或验证轴。下一步按原编号子群/固定情境检查有效块、靶表达及设计可识别性，再冻结响应检验；不按单个prior匹配声称受体机制。

## 2026-10-12T10:38:45.707906+08:00 新增公共源预处理已恢复

本聊天继续已下载9源的预处理；prepare_v3已在核验过的备用卷启动，原始输入保持原位。10:27快照44文库/0失败、275653保留细胞；10:31观察66文库404664。完整独立审计/来源验收待完成，未做跨研究整合。[进度](inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint2.md)。其他研究线按其专属最新记录，不恢复公共矩阵下载。

## 2026-10-12T10:38:40.719500+08:00 全基因准备完成，下一项患者证据字段统一

生产02:18:14—02:20:18 UTC exit0；独立审计02:23:49—02:25:47 UTC exit0；02:26:59 UTC登记v17，0活跃/0冲突。285360细胞/164文库/1121单位的139558612值全部复算；336候选中249条匹配当前LR单体先验、55复合物未由平面先验组装、32当前先验无精确LR；321有配体靶先验。1008覆盖记录959可描述、45缺靶先验、4接收家族未捕获。没有新模型/p值/正式状态/轴验证。[真实交接](inventories/analysis/2026-10-12/fullgene_readiness_v17/public_checkpoint_v17.json)。

45源header/180家族覆盖独立审计通过。患者字段差异确认：194legacy状态来自已审官方父侧表，analysis_block逐194键完全匹配；402已审文库对应175来源内患者块。旧analysis_unit仍含文库键（总304个键），不用于计患者；不能把175称跨GSE完全独立人群。新覆盖表仅补有效证据标签，所有原字段/分组/身份不改，代码已[冻结](inventories/analysis/2026-10-12/patient_block_provenance_overlay_plan_v1.json)，真实执行尚未发生。

## 2026-10-12T10:20:05.804843+08:00 完整先验通过，真实全基因生产已开始

资源v4实际02:17:24—02:17:27 UTC exit0；3文件官方大小/MD5/本机SHA全一致，未重下载公网资源，旧v3失败保留。生产02:18:14 UTC实际开始，本次观察2026-10-12T02:20:05.804843+00:00状态RUNNING、worker/child命令核实，科学登记v16、活跃2；独立全基因审计尚未启动。[进程及资源证据](inventories/analysis/2026-10-12/fullgene_resume_running_checkpoint_v1.json)。

45来源header/既有覆盖已完成：4,410,998原生源细胞、1251冻结文库、13节点、19783—58336基因；至少一个细胞对分层达20样本的3来源、达10的16来源。仅计数门槛，非新统计/临床/全count验收，独立header复核待执行。[汇总](inventories/analysis/2026-10-12/reviewed45_native_headers_v1/inventory.json)。旧frozen121患者文库字段不代表后续家族临床review当前全部审阅状态，正在核对。

## 2026-10-12T10:11:30.425815+08:00 用户已恢复任务，资源v4代码准备

UTC2026-10-12T02:08:31+00:00实测组配额1495G/2048G，粗余量553G；v16/active_jobs空、原3源生产/audit输出尚不存在。旧资源v3仍失败，不改旧终态；本机3先验重新大小/MD5/SHA核实一致。[新合同](inventories/analysis/2026-10-12/fullgene_offline_resource_recovery_v4_spec.json)/[冻结](inventories/analysis/2026-10-12/fullgene_receiver_target_code_freeze_v4.json)已准备，实际传输/生产/audit尚未开始。下一步离线资源验收→3源全基因→独立审计，同时执行既定45源header/覆盖审阅。原172/9CM/336/15/1008不变；未新增机制/状态结论。昨日暂停已由用户今日“继续推进”撤销，当前不自动扩大全局整合或启动新的大规模数据下载。

## 2026-10-11T23:50:48.128458+08:00 今日转为规划，科学下载暂停

用户报告已删除整个GPU传输目录，拟2026-10-12申请配额。只读取证2026-10-11T23:45:06.536302+08:00：项目根目录gpu_transfer_20260903不存在，Project_v5对应目录和canonical Zarr存在，v16登记证据路径缺失0、登记active_jobs为空。仅确认路径存在，未重验全文件SHA或所有历史引用，未复测删除后quota，不能宣称全目录均冗余或容量已达标。代理未删除或恢复文件。

今天没有下载、科学输入传输、资源重试或新科学任务。后续以[主线/支线计划](ROADMAP.md)为准：原9CM生态和细胞对响应双入口，身份/状态、空间、免疫治疗验证为支线；仍未建立CM驱动或耐药机制。恢复须用户指令、quota/输入验收与新资源终态。下方早期运行/下载描述为历史节点。[当前只读证据](inventories/analysis/2026-10-11/planning_pause_checkpoint_v1.json)。

## 2026-10-11T22:09:31.866016+08:00 组配额失败已定位，手动清理待执行

data4分析组达到1500G硬配额，资源aria2实际21:51退出；v16/0活跃/0冲突，完整基因生产未开始。21:53旧尾行被误判活跃，现已更正。日志Git通过/tmp副本恢复普通push；约4.961GiB冗余候选完成SHA及只读preview，用户手动清理，代理未删除。 [事件](logs/incidents/2026-10-11_group-quota-and-manual-cleanup.md)。

## 2026-10-11T21:46:20.912081+08:00 aria2恢复已运行

已核实worker和子进程，旧连续前缀约25MB，日志30MiB/249MiB，传输有增长；target未通过全文件校验，科学计算未启动。生物代码/输入SHA已装入服务器，证据仍v16。 [run](runs/20261011_fullgene_receiver_target_readiness_v1.md)。

## 2026-10-11T21:39:13.099905+08:00 资源断流已处理

v2资源worker因反复短正文且完整进度停滞停止，PID已不存在/登记无活跃作业/0冲突；TF和LR完整、target仍未全文件校验。aria2 v3已准备但未派发，科学计算未开始。 [run](runs/20261011_fullgene_receiver_target_readiness_v1.md)。

## 2026-10-11T21:23:02.879031+08:00 完整基因代码已冻结

资源恢复具体worker核实活跃，TF完整校验、靶继续分块。3来源实际9/目录13原子群；代码修正缺捕获NA和TF先验标签，合成PASS、真实全基因聚合未启动，科学证据仍v16。 [run](runs/20261011_fullgene_receiver_target_readiness_v1.md)。

## 公共新增QC驗收完成：2026-10-11 21:18 UTC+8

[完整驗收報告](inventories/runtime/2026-10-11/QC_ACCEPTANCE.md)：四份993,312細胞/149文庫、全部原始計數和4×36條件PASS，已登記technical_QC_ACCEPTED_source_final_pending。腫瘤組織576,512、CAR-T產品416,800分開；GSE292046為兩患者/8文庫富集snRNA。頭頸癌7,463細胞的post2患者解析遗漏已以v2元數據覆蓋表修正並獨立驗收，原矩陣/舊記錄保留。

可以開始來源內前處理/marker審查；低doublet檢出不等於敏感性已驗證。正式新增source acceptance仍0，原研究別名與基線文庫重疊審查未閉環；全量整合仍待決定。下方running節為保留啟動歷史，其他科學線按各自最新記錄。


## 公共新增QC全量驗收已啟動：2026-10-11 21:11 UTC+8

四份993,312細胞/149文庫逐塊比對原始計數、X/counts及全部QC指標、ID/文庫/患者映射；單線程只讀，不重算过滤。目前running，尚未簽發技術或來源驗收。[run](runs/20261011_expansion_QC_acceptance_v1.md)。其他科學分支按各自記錄。

## 2026-10-11T21:06:33.599638+08:00 完整基因资源检查点

本机三先验大小/MD5通过；服务器资源v1失败、仍缺完整TF和靶；修正后合成检查PASS。真实生产/审计未开始，证据v16、无活跃登记作业。 [run](runs/20261011_fullgene_receiver_target_readiness_v1.md)。

## 2026-10-11T20:45:16.867675+08:00 完整基因准备

完整基因/靶先验准备已预声明，336/15固定记录、285360细胞；资源获取中，科学生产尚未启动。 [run](runs/20261011_fullgene_receiver_target_readiness_v1.md)。

## 2026-10-11T20:29:43.368494+08:00 接收程序验收

70接收程序模型完成原生/图/10000条件null独立复算及补充验收，证据v16/0冲突；22可估计、0名义p05/0主BH显著。通用程序未支持配体邻居效应，目标仍active，未验证通讯轴或正式状态。 [run](runs/20261011_spatial_receiver_programme_exposure_v1.md)。

## 2026-10-11T20:26:04.291574+08:00 接收程序终态检查

70主模型生产和独立复算均exit0，22可估计、没有名义或主BH显著；尚待附加验收/登记。 [run](runs/20261011_spatial_receiver_programme_exposure_v1.md)。

## 2026-10-11T20:07:13 接收程序阶段

接收程序—邻居暴露已预声明，70主模型/75映射，尚未启动；上轮0空间主BH不冒充通讯验证。 [run](runs/20261011_spatial_receiver_programme_exposure_v1.md)。


## 2026-10-11 19:54 条件空间验收

全45空间条件邻近已独立审计并登记v15/0冲突，29可估计、16不可测、主BH<0.05为0。来源/癌种限定线索保留，未验证通讯机制或CMdriver，整体目标未完成。 [run](runs/20261011_spatial_conditional_neighborhood_v1.md)。


## 2026-10-11T19:43:34 条件空间阶段

条件空间邻近已预声明，尚未启动；全45接口和65映射，单位/背景/绑定限制保留。 [run](runs/20261011_spatial_conditional_neighborhood_v1.md)。


## 空间原生端点预检：2026-10-11 19:27:55 UTC+8

[正式run](runs/20261011_spatial_contact_native_preflight_v3.md)：全部13×5/254,458细胞完成6口径独立复算，正确181,056作者注释全匹配/零孤儿；v1缓存类型失败及v2作者键方向错误完整保留，v14明确替代错向作者结果，0冲突、本里程碑作业终止。

作者粗谱系+背景口径13记录—文件组合、7原候选、4LR在3文件具备描述性覆盖；严格纯RNA谱系口径0组合达预设每FOV双方≥20门槛。不是0表达或邻近阴性；尚未进行空间p值、蛋白结合、患者级验证或CM驱动判定。原172节点/9CM/BH/身份状态冻结。下一步固定接口的物理/单位核查及条件空间零模型。下方18:52“下一空间预检”描述为保留历史。


## 公共下载/QC最新检查：2026-10-11 18:56—19:02 UTC+8

本节更新下载线，科学家族/空间线仍以各自较晚记录为准。[可读报告](inventories/runtime/2026-10-11/DOWNLOAD_READINESS.md)、[账本审计](inventories/runtime/2026-10-11/qc_ledger_readiness_v1.json)。

- 完整源库存11来源/5,341,204细胞；QC增至4份/993,312细胞，149文库的doublet和mt均完成；账本与Zarr结构一致。正式新增验收仍未登记（0），尚未达到新增450万—550万目标。
- **下载已停止**：当前无下载/前处理worker。第二批8/14文件有校验，6未完成；Zenodo39%未通过官方MD5。日志约7—7.8小时不更新，HEAD超时、aria2超时及重试耗尽有记录；v6监控的“续传中”描述不能作为运行证据。
- 可进入下一步：4份QC分片逐来源验收、5份原始计数Zarr逐文库作者QC审查/QC；TNBC/CTCL先解来源疑点。未完成counts只先审元数据。全量跨研究整合仍待验收和用户决定。
- 按用户“每一次任务及时push”要求，本节点建立独立任务日志、脱敏证据并普通commit/push，远端SHA另核验；没有创建自动定时push。


## 空间输入修复完成：2026-10-11 18:52:59 UTC+8

[正式run](runs/20261011_spatial_native_input_repair_v1.md)与[终态证据](inventories/analysis/2026-10-11/spatial_input_repair_checkpoint_v1.json)：14基因轴筛查，TMA2新对象39123×1000的原CSV全39,123,000计数独立比对通过，X/counts两层、旧文件及坐标核查通过。v1错误的no-layers断言失败已保留，v2成功；登记证据v12、0冲突、本里程碑作业全部终止。

现在解除该对象的索引列输入修复门槛；全部13候选/5文件的原生端点、背景和患者/core预检是下一任务，邻近和机制尚未验证。v12不是注释版本，不改变原v20261010/172节点/9CM/全发现BH。下方18:36的待修复描述是保留的早期状态。

## 空间输入检查补充：2026-10-11 18:36:08 UTC+8

[只读检查](inventories/analysis/2026-10-11/spatial_metadata_checkpoint_v1.json)完成全部5份CosMx文件/254,458细胞。已确认GSE299786肺腺癌TMA2的CSV索引列误入表达矩阵，39,123行原始比对一致；该文件暂停归一化和空间推断，待新版本修复审计。详见[故障记录](logs/incidents/2026-10-11_spatial-csv-index-column.md)。原13候选和证据v11不变；端点RNA、背景、患者/core映射、独立物理标定及邻近检验仍未完成。本节更新输入可用性，不晋升科学版本。

## 当前家族分析状态 2026年10月11日18时20分 UTC+8

本节仅更新CM相关扩展分析；下载数量仍以11:22:53 UTC+8库存快照为准。[本轮证据](inventories/analysis/2026-10-11/README.md)与[正式运行](runs/20261011_family_candidate_validation_v1.md)优先于下方早期CM快照。

- 全336参考支持记录完成原生表达和新增条件模型检查；112原重采样稳定，28保留新增条件关联。
- 自动15条后续记录、11种LR组合：13接触/空间表达检查、2程序特异检查。15记录身份内部3,710模型、136块级检查完成，80可估计。
- 当前本轮登记作业均终止，审计/QR复算通过，0登记冲突；尚未启动下一空间原生RNA/坐标/患者预检。其他服务器进程未由本次快照概括。
- 保留来源/癌种限定生态关系，原编号和审阅家族并列，身份与程序分开；无正式状态赋值、独立机制或CM驱动声明。证据v11不改变注释基线v20261010。

## 当前权威状态：2026-10-11 11:22:53（UTC+8）

证据入口：[实时取证快照报告](inventories/runtime/2026-10-11/REPORT.md)。本节优先于下方历史状态；时间经过不代表新进度已经完成。

| 工作线 | 实际阶段 | 剩余门槛 |
|---|---|---|
| 新公共临床单细胞收集 | 11个来源库存合计5,341,204细胞；第二批及Zenodo继续下载 | 原研究/别名、文库、基线重叠、QC和临床端点验收 |
| 逐来源前处理 | 3份QC分片608,684细胞；5份原始计数Zarr合计2,345,161源细胞；头颈癌转换进程仍在运行 | QC分片不等于正式新增验收；正式新增登记0；跨研究整合待用户决定 |
| V3表达标签 | 4,676,787细胞全范围注释写回和读回审计完成；单基因描述命名保留既有身份 | 描述marker不自动构成新细胞亚型或CNV恶性证据 |
| V3新CoVarNet/LIANA | 4,534,433 included cells；CoVarNet 1,290样本/172节点/K=9；LIANA 31,150,424行，CM内13,729,410行；图表审阅完成 | K最优性、来源/分母敏感性、生态与独立临床验证未完成 |
| CM来源内验证 | 首轮审计完成；121已审样本/45患者单位；6个可估计主模型分层，BH显著0 | 患者、癌种、组织、分选、端点映射不全；不得声称机制或预后已验证 |
| 最终细身份 | reviewed V10 assembly索引有4,676,787行、926,810细身份非空 | GSE299623 40群/54,099细胞无最终决定；NOT_READY_FOR_FINAL_IDENTITY_RELEASE |
| Codex Cloud | 仓库提供轻量代码、离线测试、状态快照和任务交接入口 | 用户在Cloud中连接此GitHub仓库并创建/Publish环境；尚未配置或执行云端任务 |

下载并非全部健康：胃癌GSE315928的tar停在12%，日志超过10分钟未更新；其他文件在快照时有0B/s读数，队列进程仍在，不能仅凭进程存在认定传输健康。见[传输事件记录](logs/incidents/2026-10-11_download-stall-observation.md)。

本次没有观察到Project_v3科学流水线进程；代理桥进程单列为基础设施。后续最终交接登记及独立审计优先于初始manifest的pending状态。

旧V7投稿冻结、V8/V9候选和20260923被降级的CoVarNet/LIANA记录仍保留。新的表达分支具有独立输入/审计，不自动晋升为全局最终身份，也不关闭旧分支的科学门槛。

## 2026-09-25 V7 broad CoVarNet/LIANA 历史分支审计

- `20260923_v7broad_covarnet_liana` 在 2026-09-23 完成计算：1,297 个样本全部终态，
  1,268 success、29 `skipped_lt2_groups`、0 failed；完成提交为 `3c2234e`。
- 后续 `20260923_annotation_master_v1` 纠正了旧 operational 输入把 strict 比较资格误作
  大类命名门槛的问题，并明确将该 CoVarNet/LIANA 分支降级为历史审计，不能作为当前正式
  V7-derived downstream release；`REL-003` 仍未完成。
- 当前分析工作树中该目录已缺失，81 个 tracked 文件显示为未提交删除；它们仍存在于当前
  Git tree（19,609,822 bytes），其中 aggregate 子树 28 files / 12,477,538 bytes。
- 当时被 `.gitignore` 排除的 native per-sample LIANA 与 metadata cache 已不在磁盘，无法从
  Git 直接恢复。外部 CoVarNet 目录仍在，35 files / 315,143,761 bytes，核心合同文件齐全。
- 删除执行者与意图没有直接证据；本次只记录事实，未恢复文件，也未提交分析仓库的删除。
  详见 [run record](runs/20260923_v7broad_covarnet_liana.md) 与
  [incident](logs/incidents/2026-09-25_v7broad-liana-worktree-deletion.md)。

## V8/V9 轻量发布审计（2026-09-22）

- 已建立 `annotation_v8_lineage_audit_20260917`，保存 2026-09-16 生成的 V8 状态和五个谱系 manifest。全局 V8 仍为 `NOT_FROZEN`，没有谱系达到 `HELDOUT_VALIDATED`。
- T/NK 为 `PILOT_FROZEN_HELDOUT_COMPROMISED`；GSE116256 参与冻结决定，不能再作为该谱系的 held-out 评估。Myeloid/DC 的 12 个后续 held-out GSE 在 split 冻结前已产生输出。
- Fibroblast、Endothelial 的 held-out exposure 为 clean，但仍保留 clustering/panel 局限；B/Plasma 为 `DEVELOPMENT_FROZEN`，45.4% development cells 回退到 competing lineage。
- 已建立 `annotation_v9_candidate_20260919`，只复制最新 V9 候选报告和 15-row 开放 L3 复核模板。含 78 个内部绝对路径 key 的原始 manifest 仅以逻辑路径、大小和 SHA-256 索引。
- 上述 release 都不改变 V7 的全局身份/CNV 权威，也不授权 V8/V9 downstream 写回。

## Nature 投稿上游冻结（2026-09-22）

- 已建立 `manuscript_upstream_freeze_nature_20260922`：V7 是唯一全局身份/CNV 投稿参考；14.2-GB raw-count Zarr、V2/V5 父对象、基因坐标、实际 runner、动态 CNV 引擎和最终对象均已用哈希锁定。
- 外部 2026-09-19 V9 状态为 `FROZEN_CANDIDATE`，不是正式冻结：15 个 L3 cluster 待人工复核，632 个样本建议 CNV 重跑但未执行，S6 主要继承 V7 CNV。V9 已隔离为候选，不能进入正式投稿分母或下游。
- V8 仍仅有 Myeloid/DC、T/NK、Fibroblast、Endothelial、B/Plasma 的 lineage-taxonomy freeze；没有全局细胞级写回。
- Nature Article 主稿、Methods 映射、初投清单、cover letter、Reporting Summary 工作表、Data/Code availability、figure-panel source 和 supplementary 骨架已经准备。
- 当前包不是可直接投稿版本。尚需完成版本晋升决定、V9/V8 科学门槛、最终身份版本上的 downstream 决策、Copykat `hg20` 核验、54-GSE source/licence manifest、历史环境补档和作者/伦理/基金信息。
- 权威入口：[`releases/manuscript_upstream_freeze_nature_20260922/README.md`](releases/manuscript_upstream_freeze_nature_20260922/README.md)。


## B/Plasma full-cohort biological freeze (2026-09-14)

- `20260914_bplasma_full_v1` completed the full frozen V7 B/Plasma universe: 981 success, 331 explicit `<50-cell` skips, 0 failed, 685,843 processed cells and 4,378 clusters.
- Development-only review excluded the sealed 12-GSE benchmark and covered 557,601 unique cells, 3,501 clusters, 783 samples and 39 GSE; 11/11 contracts passed.
- `FROZEN_V8_B_PLASMA_TAXONOMY` freezes Bcell, Plasma and Bcell-naive. Plasma-IgA falls back to Plasma under the preregistered >=90% universal rule; IFN/cycling are overlays; boundary/competing cells remain unresolved.
- P1-1 is DONE/FROZEN. P1-2 held-out reference benchmark is the sole next biological gate. Global V8 remains NOT_FROZEN and V7 remains default.


## Biological Freeze Sprint authoritative correction (2026-09-13 21:20 +0800)

- FROZEN_V8_FIBROBLAST_TAXONOMY now exists. The fixed three-GSE run completed 6/6 samples, 10,147 cells and 30 clusters; source mapping was exact for 10,147/10,147 cells.
- Freeze Fibroblast parent and Fibroblast_myCAF only: 2,748 cells map to Fibroblast_myCAF, 4,872 to parent Fibroblast, and 2,527 contamination-gate cells remain unresolved_stromal.
- iCAF, ECM-remodeling, antigen-presentation-like, IFN and cycling remain parent/overlay evidence and are not downstream nodes. P0-5 is DONE/FROZEN; P1-1 Endothelial/B-Plasma validation is next.
- Global V8 remains NOT_FROZEN; V7 remains default until the remaining scientific gates complete.

- `FROZEN_V8_MYELOID_TAXONOMY` now exists. It freezes lineage taxonomy/rules only, not the global V8 cell-level reference.
- Independent third GSE GSE274229 completed 2/2 samples and 2,117 target cells; scATOMIC was genuinely rerun on 17,109/17,109 cells and exactly mapped to all targets.
- Accepted: FCN1 monocyte, Macrophage, cDC2 identities; C1QC and inflammatory states. SPP1/TREM2 and pDC fall back to parents.
- The 54-GSE full-cohort review finished with 1,075 success, 246 explicit skips, 5,048 cluster rows and 0 failed; it is supporting recurrence evidence only (`writeback=false`, `independent_benchmark=false`).
- `FROZEN_V8_TNK_TAXONOMY` now exists. It freezes T/NK lineage taxonomy and selected state boundaries only; it does not write back a global V8 reference.
- T/NK primary run `20260913_tnk_3gse_v3` completed 6/6 samples across three GSE: 14,185 cells and 35 clusters, exact source mapping 14,185/14,185, source T/NK/other 12,534/1,258/393, and parent/boundary/resolved-unresolved/conflict 10,664/2,845/180/103.
- No T/NK fine identity passed cross-GSE recurrence. Tissue-resident-like was accepted on the state axis with combined support across 3 GSE/5 samples/6,574 cells. Genuine scATOMIC evidence mapped 4,719/4,719 cells; 4,413 were confident and 3,705 supported the lineage mapping.
- P0-3 is now `DONE/FROZEN`; P0-5 Fibroblast refinement is the next active scientific gate.
- Global V8 remains `NOT_FROZEN`; V7 remains the default until Fibroblast, light Endothelial/B-Plasma validation, held-out benchmark, V8 malignancy rebuild and final robustness pass.
- Any older `running` Myeloid snapshot below is historical and must not override this section.

## Endothelial biological freeze authoritative correction (2026-09-14)

- `FROZEN_V8_ENDOTHELIAL_TAXONOMY` now exists: 6/6 samples, 8,951 cells and 35 clusters.
- Endothelial identity passed for 7,881 cells across 3 GSE/6 samples; genuine scATOMIC independently supported 1,538/1,547 S37 cells.
- Freeze Endothelial parent, Endothelial_arterial and Endothelial_venous. Final assignments: 4,583 parent, 1,580 arterial, 1,718 venous and 1,070 unresolved. Capillary, angiogenic-tip and lymphatic fall back to parent/overlay evidence.
- P1-1 is now DONE/FROZEN after full-cohort B/Plasma validation. Global V8 remains NOT_FROZEN and V7 remains default.

# 当前权威状态（2026-09-13）


## V8 后台实时快照（2026-09-13 09:03:20 +0800）

- `20260913_084700_mydc_full_v1` 仍为 `running`；样本账本总数 1321，其中 93 success、
  41 explicit skipped、1 running、0 failed、
  981 pending、205 pending-too-few-scope。
- 已形成 terminal success/explicit-skip 的样本为 134/1321；当前数字是运行进度，
  不是注释准确率、coverage 或 biological freeze。
- launcher PID 559 与 worker PID 988 均存活；launcher PPID=1、
  PGID=559、SID=559；worker nice=10、
  RSS=2998676 KiB。当前未发现 failed sample。
- 机器可读快照：`runs/snapshots/20260913_084700_mydc_full_v1_20260913_090320.json`（SHA-256 `b50ead5b8dcfd8c405c8b1da2631382adb2fab88ffab68f017b4f9f2c65a2859`）。V7 仍是默认项目参考，V8 run 不写回。

## V8 全量后台执行（2026-09-13）

- 已启动 `20260913_084700_mydc_full_v1`：覆盖冻结 V7 中 Myeloid/DC identity 与
  Myeloid/DC-related uncertainty 的并集，共 882,415 cells、1,321 samples、54 GSE。
- 输入为与 V7 4,676,787 cell index 完全一致的 359-gene marker 对象；X 明确为基于完整
  source library denominator 的 `log1p(CP10K)`。marker-only counts layer 不用于 counts/QC。
- 真实 smoke v3 为 2/2 samples success；identity/state 已拆分，margin <0.1 的近邻候选退回
  父级，`malignant`、`CNV_negative`、`discordant`、`not_evaluable` 未进入 cell type。
- 后台 PID 559，PPID=1、PGID=SID=559，使用 `setsid + flock`；1 worker、8 numerical
  threads、nice 10、128 GiB 虚拟内存软上限。08:51 快照为 25 success、30 explicit skips、
  1 running；子进程 RSS 约 2.99 GiB。
- 首次 launcher 的旧版 `ps -o` 兼容检查在进程成功启动后误报失败；未重复提交。PID/cmdline/
  start ticks/SID 经独立核验后补写 launcher manifest，launcher 已修复并通过合同测试。
- 结果仍为 `review_only`、`writeback=false`、`independent_benchmark=false`；V7 保持默认。


## V8 实际执行更新（2026-09-12）

- 已完成 Myeloid/DC review-only pilot 的 V1–V3 迭代。V1 因 marker-only universe
  混入 epithelial/stromal/cycling cells 被否决；V2 因错误回退到 V5 uncertainty 对象且曾有
  fail-open 汇总缺陷被否决；V3 直接以冻结 V7 obs-only 对象界定候选宇宙，并采用 fail-closed
  验证。
- V3 已在 2 个 GSE、4 个样本上完成：GSE161529 两个样本与 GSE131907 的
  EFFUSION_12、NS_19；合计 22 个 review clusters、8,997 个 cluster-level cell counts。
  `candidate_monocyte_FCN1`（4/4 samples）与 `candidate_macrophage_C1QC`
  （2 GSE、2 samples）出现跨 GSE 复现；该结果不是 held-out accuracy，也不写回 V7。
- GSE131907 adapter 的两张逐细胞 label 表及合并表重跑前后 SHA-256 完全一致；V7 输入
  SHA-256 为 `e5b68d7e5f127bcdcdadc3aaa94a9bce027ec032e6b82e68e5847ea5ccf7b38e`。
- scANVI/Census 隔离运行环境已固定到 Python 3.11.14、`cellxgene-census==1.17.0`；默认镜像
  因依赖解析到 NumPy source build 失败，唯一一次官方 PyPI binary-only retry 因 TLS
  `UNEXPECTED_EOF_WHILE_READING` 失败。状态更新为
  `BLOCKED_NO_USABLE_PINNED_CLIENT_AFTER_BINARY_RETRY`；未修改既有 conda 环境、未读取
  Census expression、未训练或预测。
- GitHub SSH 已通过显式 `id_ed25519` 验证；`main` 已成功推送并核对远端提交。`REC-001`
  不再是认证阻断。
- V7 仍是 conservative backbone；以上工作均为 V8 development/review-only 证据，不构成
  biological freeze，不覆盖正式 downstream。


- 当前细胞注释与 CNV 恶性参考：V7，状态 `FROZEN_PROJECT_REFERENCE`。
- 权威 run：`20260905_v7_cnv_rerun_freeze_v1`。
- V8 full-metadata software candidate 的 rerun2 已通过完整性审计，状态
  `PASS_FULL_METADATA_SOFTWARE_CANDIDATE_V7` / `NOT_FROZEN`；它保持 4,676,787
  cells、1,322 samples、54 GSE、43 cancers，以及 128 个冻结 V7 字段合同和 V7
  cell-index hash。该候选仍以 V7 为默认项目参考，不自动晋升，也不允许正式 downstream
  mutation。
- V8 root audit 与机器可读 release gate 记录
  `PASS_DEVELOPMENT_SOFTWARE_CANDIDATE_NOT_FROZEN`、V7 default 和
  `pytest_passed: 167`。该 167 项测试是外部 root gate 所记录的结果；本记录仓库只
  复核了其外部产物路径、大小和 SHA-256，并未重跑 Project_v3 测试。
- 独立 scANVI reference 状态为
  `BLOCKED_NO_USABLE_PINNED_CLIENT_AFTER_BINARY_RETRY`：隔离安装与唯一官方 PyPI
  binary-only retry 均有日志，尚无 reference mapping、模型训练或预测。仅在获得可用的固定
  client wheel/cache 或稳定官方连接后，才继续 metadata-only source/donor/GEO、ontology
  和 license 审查；该阻断解除前不得形成独立真值或性能结论。
- 上述 V8 结论仅为字段、schema、hash、round-trip 和内部一致性的软件候选验证；它不是
  独立准确率、真值、classifier execution 或生物学完成版结论。V7 及既有 V5/V7-derived
  downstream 均未被本候选覆盖，且没有形成 biological freeze。
- scATOMIC retry2 supporting salvage 已记录为 `PASS` 的只读后处理：没有重跑算法、联网、
  安装、覆盖或患者合并，且不能作为扩展资源上界或准确率证据。其最终可复现 supportive audit v2
  为 `PASS_SUPPORTIVE_ONLY`、scope-unverified，且不反序列化 candidate pickle。
- 细胞宇宙：4,676,787 cells、1,322 samples、54 GSE、43 cancers。
- 严格大类可用：4,057,641（86.76%）；严格亚型可用：2,455,424（52.50%）。
- exact Unknown/Ambiguous：129,843（2.78%）；非严格大类总计 619,146（13.24%）。
- primary malignant：515,031；uncertain：484,919；candidate not evaluable：45,284。
- inferCNVpy：938 success、384 explicit skips、0 failed。
- Copykat_python：668 success、654 explicit skips、0 failed。
- CNS 策略 44 samples 中 inferCNVpy/Copykat_python 成功 38/33；黑色素策略
  62 samples 中成功 27/18；其余主要因单样本候选或参考不足。

权威入口：

- [V7 release](releases/annotation_v7_20260905/README.md)
- [V7 run record](runs/20260905_v7_cnv_rerun_freeze.md)
- [V7 method](methods/pancancer-5m-v7-annotation-cnv.md)
- [V8 development candidate index](releases/annotation_v8_development_20260912/README.md)
- [V8 root-audit addendum](releases/annotation_v8_development_20260912/ROOT_AUDIT_ADDENDUM_20260912.md)
- [V8 candidate run record](runs/20260912_v8_full_metadata_candidate_v7_rerun2.md)

论文证据链：

- 2026-09-06 通讯—空间—临床证据链已形成报告、汇总表和 10 组 PDF/PNG 成图。
- 该证据链复用 V5 CoVarNet/LIANA，不是 V7 重跑结果；新分析的身份/CNV 默认用 V7。
- [Communication release](releases/communication_evidence_chain_20260906/README.md)

Git 状态：

- 继续维护已有 `Pancancer_log` 历史仓库，不在 TB 级分析父目录初始化第二个竞争仓库。
- Git 仅保存轻量、已审查、可复现/可投稿内容；重型对象由逻辑路径和 SHA-256 管理。
- 当前主机没有 Git LFS，单文件上限按 5 MiB 执行。

---

以下 2026-08 状态作为历史快照保留；涉及当前注释和 CNV 时以上述 V7 段落为准。
# 当前状态

- 最近更新：2026-08-08 06:16 +08:00
- 主要范围：`Project_v3` 泛癌单细胞转录组主线及 `Project_v2` 四队列免疫治疗扩展
- 权威运行：`20260725_120105_unified_pancancer_covarnet`
- 当前阶段：主发现运行已完成且带审计警告；免疫治疗四队列进入最终验收阻断，恶性/通讯 v2 与 24-GSE 下载重试仍在运行
- 运行结论：主运行 `completed_with_audit_warnings`；扩展任务为 `blocked_final_acceptance`、`in_progress`、`in_progress_with_failure`

主运行的详细处理、参数、失败恢复和证据路径见 [正式运行记录](runs/20260725_120105_unified_pancancer_covarnet.md)。

## 三个活跃扩展任务

### 四队列免疫治疗整合

- 记录：[免疫治疗四队列运行记录](runs/20260803_immunotherapy_4cohort_analysis.md)。
- sc/snRNA 元数据、Scrublet、1,440,403-cell Harmony 图谱、43-label 分层注释、GSE316195 严格恶性合并、GSE273952 19-slice 空间整合和 40-factor cell2location reference 均已完成。
- cell2location 已完成 19/19 张切片；71,398 spots、40 factors 精确合并，missing/extra spot 均为 0，finite fraction=1、negative=0、zero-sum=0。
- Phase 6 Stage 13 于 2026-08-06 21:09 失败：22 个 GSE316195 文库中 11 个因样本内参考不足而 `skipped_insufficient_cells`，但 CopyKAT 均成功。上游 malignancy audit 同时记录 `samples_with_method_failure=0`，说明这是最终验收策略与允许跳过语义不一致，不是 CNV 程序崩溃。
- `FINAL_DELIVERABLE_AUDIT.json` 未生成且当前无该分支进程，因此整体状态为 `blocked_final_acceptance_policy_conflict`。

### 恶性识别与 CM 细胞通讯

- 记录：[恶性与细胞通讯运行记录](runs/20260731_malignancy_communication.md)。
- 当前 active v2 manifest 为 721 个生物学样本、1,849,413 个细胞；GSE155698/GSE116256 下游排除和 GSE166555 从 2 个技术组重组为 25 个生物学样本均有机器审计。
- Stage 01 已顺序完成 367/721 个样本（50.9%），第 368 个 GSM5573499/sample34 正在运行 CopyKAT；inferCNV 已成功，R 进程持续高 CPU 并按分钟写 heartbeat。
- 当前 v2 的 Stage 02–08 canonical artifacts 均不存在，会在 Stage 01 后重算；旧 attempt 的阶段历史不作为本轮完成证据。

### 新数据下载

- 记录：[24-GSE 下载与重试记录](runs/20260805_external_geo_download.md)。
- 初始 122-task 下载因 NCBI FTP connection refused 仅 4 verified、118 failed；失败历史保留。
- active Slurm retry v3 为 141 tasks，当前 90 个 terminal rows（63.8%）：76 verified、7 skipped_verified、6 unavailable_upstream、1 failed。
- 唯一失败为 GSE201347 大型 RDS，14,413,317,040/19,071,156,087 bytes 后达到 transfer attempt limit；它不是 `unavailable_upstream`，需单独恢复。
- 第 91 个 GSE274229 PCADT2 matrix 在 06:16 heartbeat 为 101,860,201/194,354,151 bytes（52.4%）；完整性验证前不能计入 verified。Slurm 26090 仍为 RUNNING，数据盘 98% used、约 1002 GiB 可用。

## 主运行已完成

以下为直接由运行产物、脚本和日志核验的事实：

- Phase 06 共 33/33 个数据集通过逐数据集 QC：2,252,118 个 QC 前细胞，2,062,230 个 QC-pass cells，移除 189,888 个细胞，覆盖 797 个 sample units；695/695 个满足条件的样本成功运行 Scrublet。
- human-only 整合计划排除非人数据集 `GSE211602`，形成 32 个数据集、2,040,048 个细胞、39 个可读 archival shards 和 11,063 个共同基因。
- 最终全量拼接对象实际使用 38 个 shards、31 个 GSE、1,963,745 个细胞和 773 个样本；`GSE278694` 在这一过渡中被移除，但现有代码和日志未记录科学理由。
- Harmony、BBKNN、scVI 三种全量整合和 CellTypist 参考标签均已完成；最终对象为 1,963,745 个细胞。
- 证据化注释使用 BBKNN 图，完成 17 个 major lineages、136 个 lineage clusters 和 106 个 final annotations；最终交付审计为 `PASS`。
- CoVarNet 在过滤后 755 个样本上以 `K=9` 完成 nsNMF，得到 106 × 9 的固定权重矩阵和 9 个细胞模块；网络含 101 个节点、374 条边。
- `GSE195832`、`GSE123813` 和 `GSE169246` 均使用冻结的 106-state mapper 与 106 × 9 权重矩阵完成外部投影，三个 `FINAL_AUDIT.json` 均通过。
- 主发现运行最近一次已核验成功写入仍为 2026-07-30 11:18；2026-08-06 的进程检查确认三个后续扩展任务均活跃，不能再概括为“当前无运行中任务”。

## 结果解释边界

- scIB 100k 子集综合分数以 scVI 最高；全图 batch mixing 指标 iLISI 以 BBKNN 最高。不存在所有评价维度上的单一最佳方法，证据注释选择 BBKNN 是下游分支选择，不等于全指标冠军。
- CellTypist 使用 `Immune_All_Low.pkl`，它是免疫参考模型；在混合泛癌对象中的标签应视为辅助证据，不能替代 marker 和人工复核。
- 三个旧外部验证均为固定模块投影而非重新拟合。BH 校正后没有显著模块，当前只能报告探索性方向，不支持确定性的疗效或生存结论。
- `GSE123813` 的主分析仅使用 7 个 sort-matched 配对；全部 15 个配对仅作为受分选混杂影响的敏感性分析。
- `GSE169246` 的 treatment arm/response 来自 legacy metadata，仍需回到原始论文或官方补充材料核验；`Prog` 表示 progression timepoint，不应写成 recurrence 或 confirmed PD。
- 免疫治疗空间分支的 cell2location 结果是 spot-level 多状态丰度；跨癌种参考只支持保守组成解释。
- CM 分支当前只完成逐样本 CNV 的 50.9%，不得把旧 attempt 的 CoVarNet/LIANA 结果写成本轮 v2 完成结果。
- 下载中文件增长只证明传输活跃；只有完整性验证通过才可计为 verified，更不等于已进入核心表达分析。

## 阻塞或审计缺口

- `manifests/stage_status.tsv` 仍把后续阶段标记为 `pending`，与主运行最终产物和审计文件不一致；本记录采用产物审计后的真实状态，源 manifest 尚未改写。
- `GSE278694` 已从 32-dataset human plan 中移除，但现有 notebook/日志只有过滤操作，没有科学理由或 ADR；在补齐依据前应视为未决范围变更。
- `GSE162498` 仍停留在 extreme-scale planning-only 分支，未进入本次最终 31-dataset 主整合。
- 最终注释仍有 56 个 low-confidence、10 个 unresolved 和 5 个 ambiguous clusters，需要人工复核。
- CM active run 的 CopyKAT 实际参数记录为 `genome=hg20`，命名异常，正式发布前需要确认是否为有意约定。
- 免疫治疗最终审计把 11 个可解释的 `skipped_insufficient_cells` 当作 required-method failures，而 malignancy audit 记为 0 method failures；必须先统一验收策略并以新 attempt 重跑 Stage 13。
- 下载分支仍面临数据盘接近容量上限和部分 SRA RunInfo 上游空响应。
- GitHub 远程已配置，但服务器认证仍阻塞普通 push；本地记录仓库可继续维护。

## 下一步行动

1. 继续监控三项扩展任务；只在各自最终 audit/summary 通过后标记完成。
2. 免疫治疗分支明确参考不足的 inferCNV skip 是否为允许的 terminal 状态，记录决策后以新 attempt 重跑 Stage 13；不得改写已通过的 19-slice 丰度结果。
3. CM 通讯分支完成剩余 354 个样本，重算合并、CoVarNet、LIANA、NicheNet，并确认实际 CopyKAT `genome=hg20` 参数是否符合预期。
4. 下载分支继续 141-task retry，为 GSE201347 设计单独可恢复 attempt，完成完整性校验和最终 summary，并在预处理前人工复核 core/sidecar 与各 GSE guardrail。
5. 对齐源 `stage_status.tsv` 与最终审计状态，并补充 `GSE278694` 排除决定及依据。
6. 人工复核 low-confidence、unresolved 和 ambiguous clusters，决定 `GSE162498` extreme-scale 分支，并在认证恢复后普通 push。

## 2026-10-11 21:46 UTC+8 新增公共队列追加预处理

用户授权双细胞、环境RNA处理与整合准备。两个完整原始文库试运行经全细胞独立验收，正式Louvain/DecontX方法已冻结。9分支792文库准备进入全量处理；尚未完成全量处理或来源排重，最终新增AC仍0。CAR-T产品、T富集snRNA及排序文库保留独立上下文。见[方法](methods/doublet-ambient-preintegration-v1.md)和[运行](runs/20261011_doublet_ambient_preparation_v1.md)。
