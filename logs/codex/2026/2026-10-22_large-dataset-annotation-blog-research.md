# Large dataset 注释博客研究与设计

## 用户要求与范围

用户要求访问 Pancancer_log GitHub 或服务器，为 large dataset 注释博客整理、研究和设计。本次只读项目历史与方法、检索外部方法来源，并在仓库外建立研究方案和正文初稿；未改变注释、表达对象、科学模型、当前运行或既有任务。

用户在工作过程中确认受众为「单细胞研究者：方法与踩坑复盘」，与初稿取向一致，已写入研究方案。

## 起始状态和计划

取证日期2026-10-22 Australia/Melbourne。首次记录的精确时间为2026-10-22T04:25:38Z UTC /2026-10-22T12:25:38+08:00 UTC+8，这是过程中取证时间，不是任务开始时间。

本地起始HEAD `8e38c9f28a6515cfa408345884635979b821fd77`。起始git status有其他任务对STATUS/TODO、daily、既有Codex日志、CM02方法/run的修改和若干未跟踪脚本/收据；全部保留。STATUS明确此前public push拒绝不重试。本任务没有发起push或绕过该限制，新增本地日志待同步。

计划依次核对版本、阅读注释历史原文、研究方法角色、整理证据矩阵和配图/补实验方案、写正文、轻量检查。

## 已检查的来源和操作

- 阅读AGENTS、README、STATUS、ROADMAP、TODO、RUNTIME_SYNC和相关daily/task记录。
- GitHub连接器实际读取STATUS.md与20261010_uniform_expression_annotation_v1.md；后者blob SHA `d3daf9bdccd80ae6f70612389c2ab4ca0706bf85`。
- 阅读旧V7方法、10月10日历史系统审阅及旧分支下游禁用原文、9225目标hub审阅、模型artifact/训练重叠审阅、GSE299623 v4重建与当时policy、V10 core attachment全行审计、免疫图谱交付合同。
- 使用pwd、rg --files、git status/remote/log/rev-parse、sed/cat作最小范围本地检查；未连接服务器，未读矩阵。
- 网络检索及打开CellTypist官方说明/代码、scANVI、scArches、scIB、Seurat sketch、生物重复假阳性论文、GSE299623原论文。
- 使用write-page写作skill整理仓库外本地文稿，未创建或修改云端Page。

## 关键发现与解释

已观察：历史存在marker-only分母、测量缺失零化、旧图截取、来源自动证据重复计票等纠错；本地记录新于当前读取的GitHub STATUS。172节点不等于172成熟亚型，旧V7 strict subtype与V10 reviewed fine定义不同。hub诊断为特定目标，不能推断所有算法运行失败。

编辑建议：以来源/测量→分层证据→局部图→身份与状态→独立验证→可追溯交付为正文主线。建议三张核心图和三个最小实验，全部明确未新增执行。本轮项目数字是记录引用，未声称新独立复算或全量复现。

## 产物与修改

- `${WORKSPACE}/blog_large_dataset_annotation_20261022/研究与博客设计.md`：证据表、方法角色、章节、配图、补实验和来源。
- `${WORKSPACE}/blog_large_dataset_annotation_20261022/博客初稿.md`：研究者读者取向的中文实战初稿。
- 本任务新增日志；未修改既有状态/daily，避免干扰并发科学任务。

## 失败和未验证

部分论文open触发访问限制或PDF大小限制，已在研究方案注明摘要级用途；GEO网页触发验证。原癌肉瘤论文、scArches、Seurat与CellTypist入口成功。未将访问失败解释为来源不存在。

没有新服务器证据、真实配图、性能benchmark或独立身份真值测量。历史快照中的pending不直接转述为今天状态。外部独立验证和新补实验仅设计。

## 验证与下一步

两份文稿均已检查本地引用，失效本地链接为0。初稿4273字符（含链接等），关键规模和版本口径已对照来源。git diff --check通过；仓库validate_repository.sh实际exit0，0 errors /1 warning。该warning为全库已有内部路径模式提示，本次新增日志未加入内部服务器路径。核对本任务新增内容，无患者行、矩阵或凭据。当前本地交付；未commit/push，不宣称远端已同步。

建议下一步：用户审阅文章角度，按固定历史版本制作三张真实证据图，再编辑公开来源附录。建议Git commit message：`docs: research large dataset annotation blog and preserve version boundaries`。
