# 目录结构与实际工作区关系

## 记录仓库

```text
project-records/
├── README.md / STATUS.md / ROADMAP.md / TODO.md / CHANGELOG.md
├── docs/          # 项目定义、策略、命名与复现
├── methods/       # 当前认可的方法
├── decisions/     # ADR
├── logs/          # daily、Codex 和 incident 历史
├── runs/          # 正式运行的轻量级记录
├── inventories/   # 数据集、运行、产物索引
├── templates/     # 记录模板
├── scripts/       # 生成与验证工具
└── prompts/       # 可重复使用的任务提示
```

这是独立 Git 仓库。它不管理上级分析工作区中的代码、数据和结果，也不修改上级根目录的 `README.md` 或 `AGENTS.md`。

## 实际分析工作区

| 逻辑变量/位置 | 当前含义 | 是否进入本仓库 |
| --- | --- | --- |
| `${PROJECT_ROOT}` | PanCancer 分析工作区根目录 | 否，仅记录逻辑变量 |
| `${PROJECT_ROOT}/Project_v1` | 初始分析、重要历史数据和脚本 | 否 |
| `${PROJECT_ROOT}/Project_v2` | 后续分析、免疫治疗/空间内容、CoVarNet 代码 | 否 |
| `${PROJECT_ROOT}/Project_v3` | 当前主分析数据集、脚本、运行、报告 | 否 |
| `${RUN_ROOT}` | 某次正式分析运行的真实输出目录 | 否，仅索引 |
| `project-records/runs/` | 正式运行的轻量级可复现记录 | 是 |

## 证据引用规则

1. 文档中优先使用 `${PROJECT_ROOT}`、`${DATA_ROOT}`、`${RESULTS_ROOT}`、`${RUN_ROOT}`。
2. 数据集用稳定 accession（如 `GSE...`）标识，不因目录美观而改名。
3. 外部大型文件只记录逻辑路径、大小、必要校验和验证结论。
4. 样本级或患者级清单必须先进行隐私审查；默认不复制。
5. 外部文件发生变化时，通过新日志记录观察结果，不在本仓库伪造其 Git 历史。

## 已观察的外部权威入口

- 当前背景：`${PROJECT_ROOT}/Project_v3/docs/project_current_background_v2_downstream_v3_20260725.md`
- 当前运行索引：`${PROJECT_ROOT}/Project_v3/run_catalog/run_index.tsv`
- 当前主运行：`${PROJECT_ROOT}/Project_v3/runs/20260725_120105_unified_pancancer_covarnet/RUN_RECORD.md`

这些文件仍属于分析工作区；本仓库通过 `inventories/` 和 run record 建立索引。
