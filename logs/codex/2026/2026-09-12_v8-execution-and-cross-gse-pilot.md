# 2026-09-12 V8 实际执行与跨 GSE pilot

## 请求

继续实际推进 V8，并在服务器 SSH 配置后实时提交和推送可验证进展。

## 执行

- 审计并否决污染严重的 V1 与使用错误 V5 fallback 的 V2。
- 在冻结 V7 universe 上运行 V3，并扩展到 GSE161529、GSE131907 共 4 个样本。
- 生成 cluster evidence matrix、跨样本/跨 GSE recurrence 表、机器 gate 和中文报告。
- 重跑 GSE131907 adapter 并核对输出 SHA-256 完全一致。
- 在隔离环境固定 Census client 版本；保留安装与 TLS 失败日志，不污染既有环境。
- 使用显式 GitHub SSH identity 完成普通 push，并核对远端 `main`。

## 验证

- V3 post-hoc validator：PASS。
- cross-GSE validator：PASS；`writeback=false`、`independent_benchmark=false`。
- Python adapter/validator：`py_compile` PASS。
- 项目记录仓库验证在提交前执行，结果见本次提交输出。

## 科学边界

当前结论只支持候选标签的跨 GSE 复现，不支持 accuracy、ground truth 或 V8 freeze。
V7 继续作为 conservative backbone。
