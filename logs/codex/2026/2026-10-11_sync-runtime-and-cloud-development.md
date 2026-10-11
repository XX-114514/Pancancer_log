# Codex：运行状态同步与Cloud开发仓库

## 用户要求和授权

整理下载状况与Project_v3动态，git commit/push到明确指定的Pancancer_log；后续采用Codex云端任务。用户已授权背景下载及逐来源前处理，跨研究整合待决定。

## 开始状态与计划

公共origin/main为fabd2d9；服务器records有未推送8c6e440历史事故提交。先读取记录规则与状态，最小只读核查服务器证据，保留未推送历史，生成脱敏新快照，设计Cloud开发入口，检查后普通提交推送。

## 阅读与命令

读取AGENTS、README、STATUS、ROADMAP、TODO、数据/版本/Git政策、相关September run/incident/task记录和两份校验脚本。读取服务器CURRENT_STATE、三个新run的register/independent/final handoff/summary及新增公共队列证明与来源审计；读取server大队列/GEO/clinical skill原则。采用openai-docs skill并核查当前官方Cloud环境页面。

主要命令：`git status --short`；Git bundle保留8c6e440，`git fetch`和`git merge --ff-only FETCH_HEAD`；只读snapshot exporter；snapshot validator和report renderer；`bash scripts/cloud_check.sh`；`git diff --check`；全范围diff/privacy review；普通commit/push和远端SHA核验。真实主机/根路径只在外部参数，不写此日志。

## 关键发现

[报告](../../../inventories/runtime/2026-10-11/REPORT.md)为本次数量和状态权威：源5,341,204、待来源验收QC608,684、正式新增0；计算产物完成与最终身份发布分开。GSE29962340群/54,099细胞仍无最终决定。V3当前只有代理桥而无科学进程。胃癌tar日志停滞尚未诊断。

## 修改及原因

当前入口和任务板采用2026-10-11新节，保留全部旧记录；新增Cloud指南、同步合同、ADR-0003、Cloud首任务、四个run记录及源级索引。新增只读exporter、公开语义/隐私validator、可读report renderer、五个合成合同测试与离线验收命令。保存初始v1和改进v2快照，v2修正基础设施角色及实际LIANA审计路径。没有复制患者行/大型对象，也没有修改科学结果。

## 失败尝试与处置

GitHub连接器仅有pull权限，未用其写入。SSH22以及直接HTTPS dry-run超时，SSH443使用服务器已有身份认证成功；不换代理端口、不复制私钥。最初系统Python3.6不支持capture_output，改用服务器scanpy Python3.10。一次LIANA查询误加runs层目录，随后按实际annotation_rebuild目录核实；缺失探查不是科学结果失败。初始producer pending状态保留，最终handoff/独立audit是较晚依据。检查与实际推送结果见下方追加验收节。

## 未验证与后续

没有在Cloud发布环境/运行任务，没有重算全量科学分析、正式新增acceptance或修复下载停滞。尚未冻结每个历史worker精确commit与环境。下一个Cloud任务MAP-001只开发来源/临床映射validator；实际患者表在服务器执行。新整合、最终细身份、CNV及机制/临床发布各自通过门槛。

建议commit：`feat: publish runtime snapshot and Codex Cloud workflow`。

## 提交前验收

- 实际解释器Python3.13.9；本机python3为3.9.6，因此通过`PYTHON=python bash scripts/cloud_check.sh`运行。Cloud合同要求≥3.10。
- 五项离线合同测试全部PASS；v2 snapshot validator 0 errors；仓库检查0 errors/1 warning；历史release校验0 errors；git diff --check通过。
- warning来自8c6e440历史清单中的`annotation_rebuild@GitSHA:文件`，触发user@host的通用模式；人工核查为逻辑Git证据，不是内网地址或真实路径，因此保留原记录。
- 本次没有新增删除；完整origin/main以来的范围包含服务器原有8c6e440及本次提交。保留另一会话新建的CoVarNet文献审阅日志，不纳入本次stage。
- 最新服务器CURRENT_STATE和三个run目录重新核对；计数快照截止03:22:53 UTC，下载/转码仍可继续变化。
