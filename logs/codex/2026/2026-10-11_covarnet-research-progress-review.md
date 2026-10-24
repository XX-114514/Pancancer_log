# Codex task PanCancer CoVarNet 文献和项目进展审阅

- 日期：2026-10-11，Asia/Taipei。
- 用户要求：研究当前 PanCancer CoVarNet 项目，检索文献，判断进展、阶段、问题及推进方向。
- 初始状态：记录仓库存在 inventories/runtime、云端检查脚本和 tests 等既有未跟踪改动；保持这些文件不变。STATUS/ROADMAP 多为九月状态，十月新分析与独立终端审计需分别核对。

## 计划与检查范围

阅读 README、STATUS、ROADMAP、TODO、AGENTS、相关 methods/docs、九月最近相关日志及正式 run。重点核对十月本地校验副本中的全行身份审计、表达注释完成登记、新 CoVarNet/LIANA 审计、来源内复现、患者修复、免疫治疗、空间及接收响应结果。检索 CoVarNet、EcoTyper、泛癌生态图谱、LIANA、组成分析和患者级推断原始文献。

主要命令：rg --files、rg -n、受限 sed/cat、Python JSON 顶层摘要、git status --short；未读取凭据，未操作科学运行或修改权威表达对象。

## 已观察事实

- 十月 V10 接受身份组装有 4,676,787 行、926,810 非空 fine；GSE299623 的 40 群/54,099 细胞尚无最终身份决定。
- final_goal_handoff_registration 已登记全量表达和描述性标签写入、UMAP、新 CoVarNet/LIANA 及图形审阅完成。部分更早 ROOT_HANDOFF/run_register 副本仍写 running 或 pending，不能覆盖终端证据。
- 当前 CoVarNet 为 172 节点、1,290 样本、K=9；最佳 K 未成立。下游排除后的身份输入为 4,534,433 细胞/1,312 样本。LIANA 有 1,248 success、64 explicit skip，共 31,150,424 条记录。
- 原始 CoVarNet Pearson p 值和 specificity 背景问题已在十月修正网络中处理；修正后全局 148 节点/299 边。不能再次作为尚未修复的问题归因。
- 172 节点中 92 个单 GSE，35 accepted fine，137 broad + 描述性 marker；marker 名称不证明新亚型。
- 来源内复现、患者修正后大类通讯、作者标签下 ICI 配对以及若干接收端响应模型没有提供通过各自主要校正门槛的跨来源机制轴。肾癌空间三条趋化候选仍为小效应、单来源线索。
- 最新处理上下文模型中，CRC 接收程序有 2 条 BH 阳性记录，发送端成分没有通过；不是两个独立复现。GSE243013 肺鳞癌大类初筛记录 4 条候选/3 个分子组合，完整临床与联合检验家族仍有门槛；T 状态 3,884,160 个配置独立审计没有候选，Bassez 固定配对 1 项可估计且不显著。没有把项目概括为所有分支均无阳性。
- 直接复算当前节点表：35 accepted fine、137 broad-marker、92 single-GSE 节点；后者约占候选细胞 7.3%。四个 major 各只有一个节点，按现有 within-major 定义的存在样本频率恒为 1，这是新模型优先审阅的分母设计问题。
- 当日运行快照报告仓库可见性为 public；该事实与旧文档所述 private 默认不同。本任务不进行提交或推送。

## 解释与建议

项目处于大规模发现和技术交付完成后、生态可迁移性及机制/临床验证未闭环的阶段。优先按患者、来源、组织、分选确定可比较队列；以来源留出和完整训练链冻结检验模块；围绕覆盖充分的固定端点推进接收响应与空间/干预验证。新增数据优先作真正独立验证，不以扩大细胞总数作为主要完成指标。

## 修改与验证

本任务新增此日志、仓库外研究报告 `${WORKSPACE}/reports/pancancer_covarnet_research_review_20261011.md` 和本地证据哈希清单 `${WORKSPACE}/reports/pancancer_covarnet_review_evidence_20261011.json`。后两项为本地工作材料，未作为云端仓库附件发布。不修改既有 STATUS、ROADMAP、TODO 或历史结论，不更改科学对象，不运行重型分析。

报告为 12,215 字符；29 个链接中的本地目标均存在。节点数、单来源数、identity_status 和 marker_confidence 已直接从268节点目录筛选1–172编号并复算。证据清单保存所引用本地文件的大小和实际 SHA256，不复制患者数据。仓库验证完成，0 errors、1 个既有路径模式 warning；新文件未引入该类路径。文献比较覆盖 CoVarNet 2025、EcoTyper 2021、泛癌生态图谱 2024、Spatial EcoTyper 2026、LIANA 2022、LIANA+ 2024、Tensor-cell2cell、scCODA 及官方 MultiNicheNet。

## 未验证与限制

依据本地可检查审计副本与 2026-10-11 运行快照；未重新读取远端全部 raw 矩阵、执行所有科学模型或确认快照之后的后台进展。新方案和工期为建议，不是已经完成的计算或正式统计预注册。

## 失败尝试

初始宽范围文件检索及直接显示大型 JSON 造成工具输出截断，随后改为目标目录和 Python 顶层摘要。部分文献页面返回内部错误；后续几次网络请求持续未返回，主动终止，未把未完成检索当作新增证据。第一次报告验证的相对工作目录错误，改为明确绝对根目录后通过。没有声称文献检索穷尽或重跑全部科学模型。

## 后续与建议 commit message

先统一最新进展索引，再冻结同癌种多来源验证合同。建议提交信息：records: review current CoVarNet evidence and literature。当前仅准备审阅文件，不执行提交或推送。
