# 2026-11-02 论文式首页v2与六个分页

时间：2026-11-02 23:13 – 2026-11-03（Australia/Melbourne，UTC+11）。

## 用户要求

参照`Pancancer_log`首页模式重写README，保留现有版本；以论文形式、侧重展示而非下结论：从TME生态位/互作hub、CoVarNet的CM与算法引入，提出两个问题（CM与发生发展/免疫治疗/预后；CM内通讯与免疫浸润），展示数据组织、整合处理、UMAP与注释，展示0906版本的注释与处理以及CM07/CM08/CM01，并展示空间与TCGA验证；另建六个分页（数据集、LIANA与CM、注释逻辑、nsNMF局限、空间优先验证、CM发掘改进），利用库内日志与研究服务器图表绘制发表级图并详细说明。

## 开始前状态

`main`相对`origin/main`分叉（本地41、远端49个提交）；工作树已有14个与本任务无关的未提交修改（AGENTS/README/STATUS/ROADMAP/TODO等）。本任务保留这些修改，未执行任何git写操作。

## 检查与取证

- 仓库：`releases/communication_evidence_chain_20260906`、`releases/annotation_v7_20260905`、V8/V9发布、`docs/reports/20261031_TME_CM_review`、相关methods与免疫图谱RESULTS。
- 研究服务器（经已配置跳板，只读）：Project_v3 `runs/20260818_pancancer_5M_merge_prep`的03/04/05/08/10阶段报告与表（K=9成员权重、癌种活动、200k UMAP审阅表、元数据审计、0829 CM–浸润相关与新旧对齐）；Project_v5 `communication_evidence_chain_20260906/tables`、`spatial_expansion_20260903`的Visium跨切片汇总与Moran's I、`bulk_cohorts_20260904`的TCGA与ICI检验表。拉取的文件保存在会话临时目录，未写入仓库；服务器未做任何修改。

## 产出

- `README_20261026_CM02.md`：原首页原样存档（仅首行语言链接改为存档说明，正文逐字节一致，已用cmp核对）。
- `README.md`：论文式首页v2（摘要、引言、问题、数据、处理与注释、K=9图谱、CM01/07/08通讯、空间、bulk/ICI、开放问题、复现与版本边界、参考文献）。
- `docs/homepage/01–06`：六个分页。
- `scripts/render_homepage_v2_figures.py`：13张图（7主图+6分页图），PNG 300 dpi + PDF；`assets/publication/20261102/`及`figure_manifest.json`（SHA）。
- `inventories/publication/2026-11-02/homepage_v2/`：图中数值对应的样本/癌种/队列级汇总表；细胞级UMAP坐标只用于成图，不导出。

## 验证

- 全部图逐张目检并修正重叠/字形问题；所有相对链接与锚点脚本检查通过（0缺失）。
- 关键数值与源表交叉核对：TCGA分层Cox（CM07 1.11、CM08 1.06、CM01 0.97）、癌种内CM–浸润Spearman、ICI应答与生存表、重点轴支持比例、CosMx T3/T3b。
- `bash scripts/validate_repository.sh`结果见下。

## 边界与未验证

- 图4–7均为V5身份上的CoVarNet K=9/LIANA分支，未在V7重算；正文已注明。
- 分页4的min–max模拟为说明性合成数据，不是项目数据结果。
- 发布方式（2026-11-03，用户确认提交并推送）：本地`main`与`origin/main`分叉且有无关未提交修改，因此在基于最新`origin/main`的临时worktree中只提交本任务文件并快进推送；本地`main`、其41个未推送提交与工作树修改均未改动。`README_20261026_CM02.md`与推送前远端`README.md`除第3行外逐字节一致。

## 建议commit message

`docs: paper-style homepage v2 with six sub-pages and 0906 figures`
