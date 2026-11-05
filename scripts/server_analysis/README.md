# 2026-11服务器分析代码交接

配套[系统报告](../../docs/reports/20261105_server_analysis_audit/README.md)。这里归档原可写工作区中的代码，目录保留`任务/scripts/`层级，以维持脚本的相对根目录逻辑。不是新发布的冻结身份流水线。43个归档文件（含另一个[NMF目录](../nmf_benchmark/README.md)）的原/公开hash及修改说明在[imported_code.tsv](../../inventories/server_analysis/2026-11-05/imported_code.tsv)。

| 目录 | 角色 | 已有验证与边界 |
|---|---|---|
| `historical_0906_blisa/scripts/` | 早期0906追溯、面板审查及blisa准备入口 | 保留早期尝试；后续真实V7运行用另一个目录，不以早期缺依赖代表当前环境 |
| `v5_v7_comparison/scripts/` | 读取发布聚合表并生成比较 | 数据源为发布表，不产生逐细胞流向 |
| `v7_cm_blisa/scripts/` | V7频率输入、参考抽样、panel分类、blisa及CM聚合 | 历史数值实现；需外部冻结输入及原运行布局，先在新运行目录配置 |
| `celltypist72_pilot/scripts/` | 两免疫模型技术试点 | 全细胞免疫模型方案已由分流方案取代，保留历史 |
| `marker_routed_annotation/scripts/` | marker路由、SingleR/CellTypist候选、批处理及收集 | 4个mock编排测试＋3个合成收集测试；不自动赋最终细分身份 |

## 配置与复现

物理路径改为环境变量：`PC_ROOT`为既有数据根，`CODE_ROOT`为本仓库，`SPATIAL_ROOT`为空间数据根，`XENIUM_PILOT_ROOT`为早期Xenium试点。变量含义与已配置本地根一致，不必将真实路径提交GitHub。Python 3.10以上；科学依赖沿用原运行环境，不在本次整理时安装或更新。

V7运行还依赖原`inputs/v7_codes`编码文件、固定CellChatDB、config与患者/重复映射，以及原始canonical计数对象。它们留在服务器。CoVarNet主包装器已在[仓库原位置](../run_TME_CoVarNet_scheme_v1.R)，原作者冻结代码位于外部运行`RUN/scripts/frozen_author_CoVarNet_v1.R`，hash见[补充证据](../../inventories/server_analysis/2026-11-05/supplementary_evidence.tsv)。本次不另复制作者核心或发布患者映射。

需要重建工作副本时，先创建新的外部`RUN`，复制本目录对应脚本到`RUN/scripts/`，用既有本地输入/配置绑定逻辑根，并检查原配置中的rank_override、seed和输出位置。不要在冻结运行目录重执行会写输出的脚本。这里的Slurm入口只作交接，本次没有提交它们。

marker批处理还需在外部任务根放置最新`sources.txt`、已审查的`input_inventory.tsv`和可复用试点输出。`submit_all.py`必须从正常HPC终端运行；4个历史producer job ID只对应当时缺侧表的来源，重新提交前核对它们，优先按当前已存在侧表重新预检。旧试点重用依赖worker/policy/model/sidecar hash，代码变动导致不复用时应保留旧结果，另算新分支，不篡改旧收据。

## 本次验证命令

```bash
export CODE_ROOT="$PWD"
"${PC_ROOT}/envs/conda/scanpy/bin/python" -I \
  scripts/server_analysis/marker_routed_annotation/scripts/test_submission_flow.py
"${PC_ROOT}/envs/conda/scanpy/bin/python" -I \
  scripts/server_analysis/marker_routed_annotation/scripts/test_collection.py
```

所有归档Python以AST解析，shell以`bash -n`，R以`parse()`核对；这些不是重新运行真实注释、blisa或NMF。归档代码只作路径可配置化；NMF入口另加入新目录保护。原算法限制与本次尺度审计见系统报告，不能因语法检查通过称科学结果已通过验证。
