# 正式运行复核提示

请审计指定 run：

- Run ID:
- 外部 run root:
- 预期目标:

至少核对：

1. 输入 ID、逻辑路径、输入验证和矩阵语义；
2. 代码仓库、branch、commit 与 dirty 状态；
3. Conda、Python/R 版本和调度器信息；
4. 参数、随机种子与完整命令；
5. 开始/结束时间、退出状态和失败尝试；
6. 输出 ID、shape/行数/hash 等独立验证；
7. 结果、警告、错误、解释与未验证内容；
8. `inventories/runs.tsv` 和 `artifacts.tsv` 是否同步。

不要复制大型日志或对象。证据不足时将状态写为“未验证”，并创建 Codex task log。
