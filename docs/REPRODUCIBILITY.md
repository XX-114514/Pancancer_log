# 可复现性规范

## 最小复现单元

每次正式运行由 run ID 标识，并至少固定以下信息：

| 类别 | 必需内容 |
| --- | --- |
| 时间 | 开始、结束、时区 |
| 代码 | 仓库、branch、commit；工作树是否 dirty |
| 环境 | Conda 环境、Python/R 版本、关键依赖或环境文件标识 |
| 计算 | host 或 scheduler job ID；必要资源参数 |
| 输入 | dataset/object ID、逻辑路径、大小或 hash、输入验证 |
| 参数 | 配置文件标识、覆盖参数、随机种子 |
| 命令 | 可重新执行的完整命令，秘密值必须脱敏 |
| 输出 | artifact ID、逻辑路径、大小或 hash |
| 验证 | shape、行数、状态、关键一致性检查及结果 |
| 解释 | 结论、警告、失败和下一步 |

## 证据链

```text
dataset ID
→ input artifact ID / checksum
→ run ID
→ code commit + environment
→ command + parameters
→ output artifact ID
→ validation evidence
→ STATUS / method / decision 引用
```

`inventories/datasets.tsv`、`runs.tsv`、`artifacts.tsv` 提供跨文档连接；run record 保存执行细节。不要仅凭文件名或修改时间声明复现成功。

## 工作树状态

代码仓库不是 clean 时必须：

1. 记录 `dirty` 状态；
2. 记录相关 diff 或 patch 的安全位置；
3. 说明该状态如何影响复现；
4. 后续用正式 commit 重跑或明确接受偏差。

## 环境记录

优先记录已有锁文件或环境文件的路径与 hash，不在本仓库复制庞大 `conda list`。若运行依赖少量关键版本，可在 run record 中列出。Python/R 版本必须来自运行环境实际命令输出。

## 输入输出验证

验证应与数据类型相符，例如：

- AnnData：shape、obs/var 唯一性、layer 语义、稀疏性、关键列；
- TSV：列名、列数一致性、主键唯一性、行数；
- 分片：计划 shard 数、可读 shard 数、总 cell 数；
- 模型对象：shape、HVG/PC 数、embedding/cluster 键、缺失值；
- CoVarNet：eligible samples、退出状态、核心表格和图是否存在。

“命令退出码为 0”是必要证据之一，但不能替代输出验证。
