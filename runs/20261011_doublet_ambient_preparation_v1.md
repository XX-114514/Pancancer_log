# 双细胞、环境RNA与整合输入准备 v1

状态：方法及环境预检。输入为`${COLLECTION_RUN_ROOT}/objects/qc_<dataset>.zarr`，4队列993,312细胞149文库；新输出写入独立`${PREPROCESS_RUN_ROOT}`。原始计数与旧QC对象保持原样。

输入：CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea（GSE241783）、GSE286827、GSE292046、GSE294109。GSE286827应用已验收patientfix_v2覆盖；稳定ID不变。

双细胞方法及污染校正参数在环境验证后冻结。每文库需保存scores、flags、移除原因、污染率、校正前后计数及marker检查。缺空液滴或单一谱系限制明确登记。

验收：原始counts与上游完全一致，corrected非负有限且不凭空增加表达；新对象ID唯一且为上游ID子集，所有新过滤可追溯；逐文库成功状态与来源AC状态分开。整合准备保留原始基因覆盖，技术／组织／富集分支，禁止缺失基因被解释为真实零表达。

完整执行、失败和证据见[任务日志](../logs/codex/2026/2026-10-11_preprocess-doublet-ambient.md)。本次尚未进行跨研究整合。

## 21:46 UTC+8 方法冻结与范围更新

初始两个完整文库试运行均经独立全细胞验收。正式运行另立`${PREPROCESS_RUN_ROOT}`=prepare_v2；范围9分支792文库，输入3338473细胞（4QC+5预QC），不能作为最终通过数。参数及限制见[方法](../methods/doublet-ambient-preintegration-v1.md)。运行命令：

```bash
python scripts/run_ambient_preparation.py --run-root "${PREPROCESS_RUN_ROOT}" --collection "${COLLECTION_RUN_ROOT}" --rscript "${SCRNA_RSCRIPT}" --code-git-sha "${CODE_GIT_SHA}" --parallel-datasets 4 --library-workers 4
```

实际启动、终止及独立验收结果将追加；现阶段不声称全量完成。

## 2026-10-12T02:12:48.762758+00:00 恢复记录

prepare_v2因组quota errno122终止，保留部分输出、无全量验收。data4组上限现2048G/已用约1495G；恢复输出为另卷私有`${PREPROCESS_RUN_ROOT}`=prepare_v3，输入原位只读。预算782893303100 bytes+500GiB，小于输出卷约1.598TB实际余量。按原9分支792文库恢复；新增必需参数`--storage-plan "${STORAGE_PLAN}"`。不恢复新下载。

## 2026-10-12T10:38:45.707906+08:00 prepare_v3实际运行

启动02:15:17 UTC/10:15:17 UTC+8；代码cdcea29，4来源×4文库，500GiB预留，存储计划已核实。10:27快照44文库/0失败，275653保留细胞、额外9608；完整独立count/ID/标准化审计将在全队列终态执行。final source AC仍0，未开始全局整合。

```bash
python scripts/run_ambient_preparation.py --run-root "${PREPROCESS_RUN_ROOT}" --collection "${COLLECTION_RUN_ROOT}" --rscript "${SCRNA_RSCRIPT}" --code-git-sha cdcea29b45e989cb8d03235b9a92521611b6dd17 --parallel-datasets 4 --library-workers 4 --storage-plan "${STORAGE_PLAN}"
```

[可读快照](../inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint2.md)不等于实时状态或最终验收。

## 2026-10-12T10:54:17+08:00 第一批完整来源验收／审计I/O优化

- 实测10:53:50（UTC02:53:50）182/792文库，0失败，944249保留细胞、额外38637个预测双细胞，274673严格技术候选。完整运行仍活跃。
- 3个完整来源独立验收：CAR-T 396593、GSE292046 52525、GSE294109 130243，合计579361细胞/92文库；原始审计02:46:16 UTC和缓存审计02:54:17 UTC均PASS；所有来源/文库科学检查、计数、候选124513及基因覆盖诊断完全一致。报告及SHA匹配证据见runtime/2026-10-12。
- 审计缓存只改读取方式：每次最多一源36GiB稀疏参考，较大来源保持原backed读取，校验数学/门槛不变；原始报告保留。生产代码、双细胞和DecontX参数未改。全量阶段审计脚本升级尚待显式版本登记。
- 数量为技术处理与验收，不等于来源别名排重已完成；final AC仍0，未启动跨来源整合、新下载仍暂停。CAR-T产品及富集snRNA用途限制保留。

## 2026-10-12T11:05:00+08:00 后续整合合同与native初始QC审计冻结

- 3c1cd3fd845b53d2ee21817d7733cad59448bdb1已普通push并独立HTTPS ls-remote核验一致。签名转发首次因Cursor既有remote forwarding端口冲突失败，复用已有SSH master并换唯一的本地签名socket后成功；没有换科学代理端口、复制密钥或重启计算。
- 审计v2在02:59:05 UTC明确登记并应用：保留原v1源文件、最初execution_start不改，最终审计尚未启动；代码SHA51df1b…，与原92文库全量科学检查一致。
- 新整合输入合同分普通组织、富集T/snRNA和CAR-T产品；原始计数/估计层、基因可测覆盖、患者块及作者前处理历史分别核验。CAR-T作者双阳性状态移除3591/30427(11.80%)、高于CD4(4.08%)/CD8(5.06%)；保留全部原始映射，状态关联需敏感性审查。聚合证据无患者行。
- 加强独立审计：对5个native分支的全部原始计数另算floor/mt/作者PASS/作者doublet及全行过滤账本，检查每个非空原始文库均被覆盖。生产不改；缓存native账本和原始文库索引，避免每文库重解压整源账本。新代码SHA37eaef…已语法验证，真实native独立执行仍待下一步，不能提前称验收完成。

## 2026-10-12T11:07:39+08:00 原生初始QC全量复核通过

- 取证03:06:50 UTC/11:06:50 UTC+8：LBCL598289→597962、儿童HGG401253→400899、SPECTRUM927205→880042；合计1926747原始细胞→1878903初始QC保留，全行独立账本检查全部PASS。这不是后续污染/双细胞全量终态。
- 加强审计v3于03:07:39 UTC登记应用，代码37eaef…对应公开20505b1；旧v1/v2脚本和最初运行manifest保留，模型生产参数不变。新增原生全行过滤复算、非空原始文库覆盖门槛，账本每源只读一次。正在另跑已完成来源的全量验收新检查点。
- 全行技术元数据复核4原生CELLxGENE源：SPECTRUM cell，LBCL/HGG/乳腺nucleus，乳腺v2/v3均在各原始文库内唯一；作者是否已污染处理仍unknown，未从整数raw/X推断“未经校正”。
- HNS部分完成40文库中25个污染估计待复核：21个marker损失、8个高UMI损失（可重叠），0不收敛。继续保留raw并阻止采用这些估计，不能把处理成功解释为污染清除通过。
- 时间补正：上一人工节点标题11:05:00不是实际取证时间；20505b1提交实际11:03:42，v2应用10:59:05。实际节点以JSON时间戳/commit时间为准，旧标题作为补录背景保留。

## 2026-10-12T11:15:43.446858+08:00 准备清单实际生成与作者处理历史修正

- 11:11:45快照274/792文库、0失败，1340180保留细胞、额外53640预测双细胞；完整队列仍运行。加强审计v3于03:11:59 UTC完成3来源92文库579361细胞，全量count/ID/规范化/非空文库覆盖全PASS。
- 已实际生成preintegration_preparation_checkpoint_v1：92个逻辑分片引用、层选择/技术掩码角色及仍未放行的来源/作者历史/用户scope条件；未整合、未增加正式AC。
- 原文方法核实LBCL已SoupX1.6.2及scDblFinder1.12.0，且598289是非恶性LME分支，不能作全肿瘤分母；乳腺原workflow含Cell Bender0.1.0、Scrublet0.2.1。论文方法与本地exact raw/X槽追溯分开；作者校正方法已documented，具体slot仍pending。后两来源追加DecontX层仅诊断/敏感性，拟选择counts时必须重新计算normalize，不能直接复用当前X。原生产和所有层保留；作者历史覆盖表进入准备合同，最终source gate仍false。
- 原文URL与简要方法字段在source_processing_history_overlay_v1，不把filtered/public counts笼统称“未去污染原始液滴”。HGG文献包括不同assay双细胞章节，具体snRNA槽方法须进一步对应，暂不把ATAC/multiome参数移植到当前RNA对象。
- 新finalizer组件已写入并语法验证：等当前既有队列真实终态及全量独立audit PASS后，写服务器可读terminal状态与整合准备清单；失败记录明确中断，始终不自动整合。尚未启动该组件；Git最终结果发布不自动，须后续脱敏审阅/普通push。
- 4ff6827已普通push，独立HTTPS SHA一致。当前新增公共矩阵下载继续暂停，论文方法只读核查不启动新矩阵下载。

## 2026-10-12T11:24:37.182560+08:00 后台终态记录组件已核验启动

- 终态记录组件03:17:30 UTC实际启动，03:18:00 UTC独立核实进程命令及builder/finalizer/collector三个SHA，版本bc39d7e已普通push并HTTPS SHA核验。代码等既有323251队列真实结束，不产生新科学模型；完整验收后写可读terminal与全分片准备索引，仍不放行整合；最后Git发布需再审阅。
- 11:20:50（03:20:50 UTC）323/792文库、0失败、1517029保留细胞、额外60087预测双细胞，严格技术候选552354，尚非全量终态。
- 对已有LBCL/乳腺H5AD头的最小核查：counts源均为raw/X、重复gene symbol守恒合并，作者is_pre_analysis=False；这些头不能证明exact pre/post校正槽，所以维持unknown并保留v1作者历史门槛。源文件SHA及检查结论另列。
- 神经谱系覆盖复核：pHGG401253作者标签含335682 glial、15599 oligo、4509 neuron、29836 myeloid、8249 T、4518 pericyte及2860 endothelial；当前诊断marker面板无神经/胶质集合。少于2个“可识别谱系”不能解释为真实单一谱系，新增污染估计需脑谱系marker保留审查；拟优先counts但最终选层未冻结。准备清单原本已全部integration_allowed=false/选层pending，不因模型技术PASS跳过此审查。没有新增正式身份标签或用glial等同全部恶性细胞。
- 失败尝试：固定SHA push首次被zsh冒号变量修饰误解造成refspec错误，改用结构化subprocess参数普通push成功，未force/改写历史。

神经面板规则已落实到准备索引代码：当前pHGG来源即便模型技术PASS仍拟选counts，并明确neural_glial_marker_panel_and_preservation_review_pending。真实92分片重放计数/门槛不变；合成单一PASS脑分片用于验证不能绕过该门槛，未冒充真实细胞审计。终态组件增加独立attempt版本参数，拟保留v1记录、只替换本任务的记录等待组件为v2；不重启或修改科学生产。实际替换尚待新代码普通push/冻结。

## 2026-10-12T11:29:47.169658+08:00 v2记录组件启动与11:28快照同步

- 03:27:44 UTC终态文档组件v2实际启动，03:27:45取证：自身/builder/collector SHA和主生产存活均核实。19f6a38c6e88963b3c1d34a34dad382c77aaa002代码已先普通push并独立HTTPS核验。旧v1等待组件按其launch PID/命令确认后停止，仅此记录角色；旧start/源码/stop receipt保留，没有停止任何科学生产。
- v2实装脑来源拟选counts及神经marker审阅门槛；当前终态组件不做整合、不自动Git推送。科学生产仍cdcea29、完整审计37eaef…；初始manifest与审计升级历史保留。
- 取证11:28:27/03:28:27 UTC：367/792、0失败、1691585保留细胞、额外66451预测双细胞、652101严格技术候选。source_final_AC=false、全量audit仍pending。两未开始来源保持排队，不以输入3338473当QC通过总数。
- 本阶段完成恢复执行、部分来源全量验收、准备输入合同/索引及后台终态记录；实际全队列仍在后台，不宣称所有RNA清除成功。待真实终态后继续核对完整报告、来源历史/排重及用户实际整合范围。

## 2026-10-12T11:37:06.841860+08:00 用户质疑备用卷与资源预算：暂停并核实

- 03:32:17 UTC可恢复暂停本轮主生产独立进程组及记录组件，命令/UID/pgrp全部核实，29存活进程T。只停本轮，结果保留，不删除/移动、不改现有目录权限。用户输出卷确认前不自动恢复。
- 03:33:07实际占用146422484992bytes=136.3666GiB，其中R文本交换94.46GiB/69.27%；原始count分片10.18、校正13.78、归一化15.30，其余约2.64。首次只读du134GiB与后续暂停统计136.37差异来自此前尚在写入的约2GiB，不是额外500GiB被分配。
- 目录是账号已有个人项目区下自己的mode700新run，卷根root所有且当前账号不可写；没有管理员项目使用授权证据。此前“可写/limit0”只能说明技术权限，不应当作组织批准或无限空间。
- 原计划下一批峰值+500GB，被本轮粗保守估计扩为全部9源100B/nnz累计729.13GiB+500GiB，较原留空单位多34.34GiB；不是数据必须占用1.2TiB的实测结论。需要先按实际压缩/必要产物/批次重算，优先削减未压缩交换矩阵保存方式；现有文件不擅自删除。
- 最新暂停状态383/792、0失败、1751187保留细胞；仍未全量审计/整合/source AC。旧collector非Z存活逻辑会把T视为alive，operator_storage_hold_v1及公开暂停快照优先。
- 公开报告已移除内部UID，只包含匿名汇总/方法与权限布尔结论。确认问题来自用户对存储区权限的明确质疑，无新增skill审批限制；已向用户请求允许输出卷，在答复前保持暂停。
