# V7全数据双细胞处理溯源

用户要求检查整个数据此前是否实际剔除双细胞，尤其区分K-09本轮245个复查候选与历史工具的筛选结果。本任务先只读追溯，不将本轮95分位门槛倒用为历史工具阈值，不自动删除或改写V7。

开始取证2026-11-05 06:04:56 UTC / 14:04:56 UTC+8。主工作区15项tracked并发改动及其他未提交内容保留；读取项目规则、当前状态、V7方法/run和历史QC记录。服务器使用既有跳板及身份，读取祖先禁止删除规则；原数据和模型不变。

检查范围：冻结V7全4676787细胞、同轴全基因core、早期统一队列QC、新增队列准备及合并代码/过滤账本；核对当前Mast、T6及K-09与历史判定的交集。十月新增队列的独立DoubletDetection/DecontX分支另作版本范围核查，不将其结果自动迁移到九月V7。

初始计划：找到实际输入/输出及终端→核验工具和执行阈值→对比过滤前后ID集合与filter_reason→全量统计V7中的旧预测、保留、未评估和缺失→对245个新候选追溯旧判定→输出可复核分阶段表。代码存在或文档写“QC完成”不足以证明当前对象真的排除了旧候选。

已读初步证据：早期统一run记录33队列、695个合格样本Scrublet成功，2252118→2062230细胞；总移除189888包含多类QC，不能全部说是doublet。更晚新增9分支预处理采用独立DoubletDetection并记录保留/移除，但与V7是否同一输入尚需精确区分。下一步读取服务器原始账本与当前全量metadata，实际结论待核验。

命令采用结构化SSH调用和只读Zarr/HDF5/TSV读取；内部位置与完整stdout留私有工作目录。新公共文件仅方法、汇总与脱敏记录。科学处理未重跑，无新定时监控。

验证、失败、已完成结论、未验证范围及同步收据将在本日志连续追加。建议commit：`Audit historical doublet filtering into frozen V7`。

## 全量实际结果

06:12:30.799098 UTC / 14:12:30.799098 UTC+8完成54来源全量成员核验：4676787个V7细胞精确对应原QC过滤后对象，原分数含缺失逐项一致；54对象均有原predicted/removed字段，当前工具阳性与已剔除标记均0。当前source_obs_name与原obs index对应，逐来源当前/原QC后数量全部相同。冻结V7与core index也精确一致，非抽样。

06:13:37.529793 UTC / 14:13:37.529793 UTC+8完成原移除汇总：原33源5609个工具阳性排除，其中当前V731源5518（5177仅doublet，341与其他QC原因重叠），新增23源13821。当前54来源合计19339，实际8% cap触发0，不可将两批全部QC排除754222都说是doublet。当前旧1963745、新2713042细胞，合计4676787。

K-09的245个本轮联合候选原来全部已运行Scrublet且未被判阳性：78旧直接方法、167新增方法；不是初始未检测、cap保留或已删除后重入。K-09全1027中7个未运行，T6总156631中786，Mast45036中218。没有执行新候选的剔除或重聚类。

## 独立复核与纠正

06:16:51.715577 UTC / 14:16:51.715577 UTC+8第二套实现PASS：不复用首次indexer，核对完整ID排序多重集及来源内唯一键；全4676787分数/flags重查，旧2148018行原始cell_qc账本逐行重算，5518 doublet原因/5177 sole原因与pass数相符。新增13821通过另一个sample_qc_thresholds表交叉核对。没有重跑Scrublet或访问表达矩阵X。

06:19:50.670927 UTC / 14:19:50.670927 UTC+8完整性PASS：checkpoint SHA仍为e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e；4647788有限分数、28999 NaN、0 inf；NaN与not_run逐行精确一致。原31源102个小单元、新23源13个小单元，共115单元未运行。旧文档“只打分不剔除”不符实际代码，真实返回NaN及False占位，不能视为有效阴性。

读取的核心原文件：原33源production validation及31源cell_qc/doublet_summary/filter_reason，新增dataset_qc_summary/doublet_summary/sample_qc_thresholds，54 filtered H5AD obs，8月merge脚本，9月core构建脚本及master/restore provenance。执行器body SHA和原账本SHA入公开证据，内部路径/完整stdout不公开。

失败/限制：初始一个相对搜索路径相对cwd不成立，已用绝对本地路径纠正；采集器实际exit均0，仅有Pandas类型FutureWarning。没有以设计文档替代执行账本。未重跑模型、未验证真实残留doublet率、未恢复新增队列全部历史排除条码、未验收候选混合RNA的生物学身份。新TME独立版本的代码/文档只作范围对照，不拿其宣称的继承QC替代V7全行实证。

发布使用独立干净worktree；fetch发现并发新TME提交，已普通fast-forward到432b8a6并保留全部变更。主工作区15项tracked并发修改保留。新增[结果目录](../../../inventories/analysis/2026-11-05/v7_doublet_provenance_audit_v1/README.md)、方法/run/daily及4份实际执行body代码；T6证据页追加历史处理澄清。隐私/diff/仓库检查和普通push收据在完成后追加。

发布前实际检查：28个显式文件约409KB（含既有STATUS全文件大小），源码body SHA逐项记录、Python AST通过、隐私扫描无本机/服务器私有路径或凭据。完整仓库轻量验证0 errors/1既有warning，无断裂相对链接；git staged diff whitespace检查exit0。所有患者/样本/细胞行、矩阵和完整stdout留私有；源级汇总和代码单独暂存。普通push和独立远端SHA须在成功后另记，当前准备记录不冒称已发布。

内容commit92aa7ca的首次普通push因并发main前进而non-fast-forward拒绝。fetch到790a51c后保留两边历史普通merge；部分检出首次在默认网络下缺promisor对象，网络授权后补取成功。STATUS与run索引均为双方尾部追加冲突，已逐段同时保留并断言两块文本完整，无删除/历史改写。合并后仓库再验0 errors/1既有warning。本任务相对最新远端的28文件diff whitespace为0；merge暂存差异另报告远端已存在的annotation_gene_panel_v1.tsv空末列tab，未擅改并发panel。科学证据与四份执行body不变，后续普通push和独立SHA另记。
