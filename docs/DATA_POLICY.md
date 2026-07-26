# 数据与 GitHub 安全政策

GitHub 目标为 private 仓库。Private 不等于可以无审查上传；本政策仍按最小披露原则执行。

## 允许进入仓库

- Markdown 状态、方法、决策和日志；
- 小型 TSV/CSV/YAML 索引或配置；
- 数据集 accession、run ID、代码 commit、逻辑路径；
- 参数、命令和经过脱敏的错误摘要；
- 文件大小、行数、对象形状、hash 等轻量级验证信息；
- 不包含个人身份或受控字段的汇总结果。

## 禁止进入仓库

- 密码、API key、token、SSH 私钥、Cookie、数据库凭据；
- 患者身份信息、可重新识别个人的字段、受控访问数据；
- 原始人类测序数据；
- H5AD、HDF5、RDS、pickle、矩阵、BAM/CRAM/FASTQ 等大型或敏感对象；
- 完整作业日志、完整环境导出、缓存、模型 checkpoint；
- 未经审查的样本级元数据和内部共享链接；
- 无必要的真实用户名、绝对服务器路径、内部主机名。

## 路径与标识

文档使用 `${PROJECT_ROOT}`、`${DATA_ROOT}`、`${RESULTS_ROOT}`、`${RUN_ROOT}`。只有在本地故障定位确实需要且确认不会 push 时，才可临时记录真实路径；提交前必须脱敏。

数据集 ID、样本 ID 和运行 ID本身不自动等同于安全。若其编码包含个人信息，必须替换为非识别性映射并把映射保留在仓库外。

## 大小和格式

- Git 跟踪文件建议小于 5 MiB；更大文件必须说明理由并人工审查。
- `inventories/` 中的小表格应保留跟踪，不得被通用 ignore 规则误伤。
- 大型结果只记录 URI/逻辑路径、hash 和生成 run。

## Push 前检查

1. `bash scripts/validate_repository.sh`
2. `git status --short`
3. `git diff --cached --stat`
4. `git diff --cached`
5. 确认 remote 是用户指定的 private 仓库
6. 确认没有原始数据、凭据、敏感绝对路径和大型文件

自动检查不能替代人工审查。
