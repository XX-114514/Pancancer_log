# 当前研究工作总览：证据阅读、写作与仓库集成

## 用户要求与任务边界

系统阅读Pancancer_log，在历史0906叙事与六分页基础上新增中文`docs/homepage/07_current_work.md`，按研究工作线解释最近进展、阴性/受阻结果、版本关系、依赖与P0/P1/P2计划；更新README导航、检查事实/链接并按仓库规则同步。此次仅文档与只读取证，不启动大型计算、不重拟合CM、不修改数据、结果或分析脚本。

## 起点与安全隔离

已读取AGENTS、README中英文、首页01–06、STATUS/TODO/ROADMAP/CHANGELOG、CODEX_CLOUD_DEVELOPMENT、RUNTIME_SYNC、WEHI_SERVER_WORKFLOW及ADR-0003。目标GitHub经连接器确认public且可写。

本机主checkout为`main`，HEAD `63f78a2fb0f55935e7fe3cb21ed185ac347d3cd5`，存在大量跟踪修改和未跟踪内容（包括README、AGENTS、科学状态和脚本）。初始缓存origin为`e1900f8`，不能据此判断最新远端。建立不含未提交内容的独立克隆，将origin显式设置为用户指定目标，fetch后从远端创建`docs/current-work-20261105`；没有移动原checkout的分支、修改其文件或清理其内容。

阅读基线为远端main `7617443f90cb42bf4e772d21774731cc9fc020c6`。原checkout与此基线41个本地独有/191个远端独有提交；`--cherry-pick`识别到40个无等价patch的本地提交。检查原九来源scVI冻结方法及日志差异，确认部分本地交付未以同一SHA公开；公开已有工作线看板和分页可说明其范围。本任务没有把这些本地提交合并、发布或把私有材料写入页面。

克隆最初origin指向本地源，fetch GitHub输出的ref替换反映两种来源的不同祖先，并不是本任务执行force push。未改写任何历史。

## 阅读计划、范围与关键命令

计划按研究问题归并来源/图谱、注释、组成CM、功能状态、空间、通讯/临床、复现七线；先读最近审计及原表，再追读重要commit实际diff与方法，最后写作、数值/链接/隐私检查、普通提交与远端SHA核验。

- `git status --short --branch`、`git remote -v`、`git rev-parse`、`git rev-list --left-right --count`、`git log --left-only --cherry-pick`：核本地与远端边界。
- `git log --since=2026-10-06`：基线可达771条记录（含合并），日期覆盖至11月5日；该窗口没有10月6–10日的新提交。重点看最近50条的日期、主题和变更文件（304个唯一文件），近7天按10月30日—11月5日归并科研事件。
- 深入的diff/实现包括`fc13620`审计与代码归档、`3b5d997`基准、`0138a29`Harmony变量、`e978e41`统一排除规则、`4b8088f`CSR读取、`b0ac58d`CNV主导、`0429295`超时恢复；结合`d22c8c0`空间报告、`bee3325`数据层、`929d633`持久对象更新、`91289b2`TME科学审阅追读结果。近30天回溯`1319ad6`/`00c9e1d`阴性空间/响应结果、免疫图谱和CM02证据、来源资格及受控注释设计。
- 读取11月5日系统报告、snapshot/source_evidence/version_crosswalk/state_corrections、TME进度/旧summary、V5/V7聚合表、两seed V7对应/rank/完成报告、NMF原比较与尺度表、blisa汇总/panel资格、annotation_pilots/preflight/current_presence。
- 检查实际实现：NMF compare_methods与check_nmf_scaling、marker_policy与singler_hpca、V7 aggregate_cm_blisa、batch evaluate及新改动的整合/恶性逻辑；核查算法语义与报告限制，不执行科学拟合。
- 阅读integration_v2、analysis-ready、恶性方法/相关run，CT75报告/预设方案/WEHI重放，10月31日CM结构/CoVarNet对照报告、workstream board、状态读出/增量协议、旧免疫CM02首页及公开证据口径；最近Codex/daily日志包括11月5日三组审计/注释记录。
- 用既有跳板SSH只读检查旧服务器时间、当前用户Slurm队列及进程类型；没有改主机校验、密钥、路由、任务或持久配置。

## 时间、实时观察与无法验证的范围

本任务在2026-11-05（Australia/Melbourne）执行，阅读基线中的WEHI取证起点为00:49:42 UTC / 08:49:42 UTC+8 / 11:49:42 Melbourne；各科学事件采用原run时间，不把本次读取日期当计算完成日期。

旧服务器连接探针返回01:21:13 UTC / 09:21:13 UTC+8；队列观察01:21:37 UTC / 09:21:37 UTC+8 / 12:21:37 Melbourne：9个作业均PENDING，原因4个DependencyNeverSatisfied、1个Dependency、4个JobHeldUser，没有RUNNING条目。进程类型表还有python3/codex/node，未读完整参数，不能据类型表证明没有非调度科学进程。未将内部主机、用户名、物理路径或凭据复制进公开记录。

本次未建立WEHI实时会话；其运行状态只用归档收据。TME run_v2仅核到Harmony完成，调度状态未验证、最终summary缺失，故不用“当前运行中”。本次没有重新核服务器大矩阵、全部恶性侧表、空间bin级检验族、完整环境迁移或患者身份。

## 关键发现与写作决定

- 新TME、V7重算、CT75、十月免疫CM02和0906 V5是不同版本；每个CM数字都带范围/日期/seed语境。
- integration_v2技术完成与TME run_v2终验待补并存；来源排除后对象历史行数和合格范围不同，患者键与独立患者分开。
- 全量注释未获验收，单来源候选和文献整理不能写成最终身份；HPCA顺序/适用域及High/Low依赖关系保留。
- V7 rank/K9已有完成收据，K9是override。45个NMF任务完成，但尺度审计改变部分匹配；不作准确率排名。
- CT75主验证未通过；V7弱信号依赖成员定义且myCAF不可识别。候选边、不可测、阴性和机制证据分开。
- 预定义、数据驱动和TME状态分支分别说明；RNA增量不冒称临床或通讯验证。
- 保留旧免疫CM02的限定正向证据与0906临床不一致结果，新CM不继承其效应。
- labshare已有指定对象收据，更新旧scratch叙述；不宣布所有环境/资源完成持久化。
- `fc13620`/`7617443`已从GitHub fetch取得，因此原日志当时推送失败不等于现在仍未发布。

## 修改文件及原因

1. `docs/homepage/07_current_work.md`：七线研究叙事、快照/看板、近期科研事件、科学不确定性、依赖/验收标准、版本与证据导航；复用已归档CT75图并注明事后性质。
2. `README.md`：顶部增加中文总览入口；历史0906正文与六分页保留。
3. `README.en.md`：增加明确标注Chinese的入口；未生成不完整英文配套页。
4. 本日志及同日daily：记录真实阅读/验证/限制；不重复修改STATUS/TODO中的科学状态。

## 验证、失败尝试与待完成检查

初次sandbox Git fetch报DNS不可解析；按权限流程的正常联网fetch成功。初次sandbox SSH被网络权限限制；按权限流程沿既有跳板的只读连接成功，没有绕过身份/主机校验。以上是本次访问条件，不归因服务器科学任务失败。

首轮链接检查发现草稿引用一个不存在的恶性运行日志名，已改为实际已核查的恶性worker源码；两条本任务日志链接在日志创建前暂不存在。逐事实复核将“counts整数”收紧为读回JSON实际支持的“抽查counts为整数”。新增文档不会改变科学状态或已有对象。

后续验证结果与发布收据在同一日志追加。建议提交：`docs: add evidence-based current PanCancer research overview`。本次之外的下一步为TME终验、候选路由/身份准入、NMF固定尺度复核；未在本任务执行这些科学任务。

写作取证节点：2026-11-05T01:29:40.730082+00:00 / 2026-11-05T09:29:40.730082+08:00。

## 最终内容检查

- `bash scripts/validate_repository.sh`：0 errors / 1项既有绝对/内部路径warning；warning来自原有文件，新增文档无对应命中；全仓库本地相对Markdown链接检查通过。
- `bash scripts/validate_releases.sh`：0 errors。
- `python3 scripts/server_audit/verify_snapshot.py`：PASS，复核64项公开证据hash、43项代码hash、45个NMF任务及4,676,787细胞转移计数、空间符号检验算术；不宣称生物学验证。
- 本任务临时只读文档检查：五个交付文件249处相对链接、37处锚点、页面40个不同commit引用全部通过；旧中英文README删除新增导航段后与基线逐字一致，01–06及英文配套页没有改动。
- 对照JSON/TSV复核integration_v2规模/抽样计数验证、TME当前进度和缺失终态、旧TME分层守恒、NMF任务数及尺度不变性/改变配对数；关键结果PASS。
- `git diff --check`、本次文件凭据/内部路径/IP模式检查、单文件体积检查通过。已目视审阅引用的CT75归档图、图注和事后控制边界；没有新增图件或数据表。
- 发布前再次fetch，远端仍为阅读基线，没有并发新提交。仅显式stage五个文档文件。

检查范围为文档正确性、已有公开收据一致性与轻量仓库规范；未执行全量科学测试套件、服务器重算或新的独立身份/机制验证。科学未决问题均列在页面，不作为本次文档交付失败。

## GitHub发布收据

内容提交：`270d4e1a586f6a017cbb0715b3236dc8c2474134`。普通`git push origin HEAD:main`实际exit 0；随后独立`git ls-remote origin refs/heads/main`返回相同完整SHA，远端发布已核实。内容提交后隔离工作区干净，原本机分叉checkout未修改。

远端核验完成记录时间：2026-11-05T01:32:40.152525+00:00 / 2026-11-05T09:32:40.152525+08:00。本收据另作常规追加提交；不改历史、不改页面的阅读基线。本次最终交付仅五个文档文件。
