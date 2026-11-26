# 2026-11-27 首页论文式重写（v3）

事件/取证时间：2026-11-26T13:00Z–2026-11-27（UTC）；墨尔本 2026-11-27 00:00 起。执行者：Claude（云端会话，无研究服务器访问）。

## 用户要求

按现有首页模式重写README（保留现有版本），以论文形式、侧重展示而非下结论：生态位/互作枢纽 → CM定义与算法 → 两个问题（CM与发生发展/免疫治疗/预后；CM内通讯与免疫浸润）→ 数据组织 → 整合处理、UMAP与0906注释 → CM07/CM08/CM01 → 空间与TCGA验证；六个分页（数据集、LIANA与CM、注释逻辑、nsNMF局限与cNMF/状态模块、空间优先验证、CM发掘改进）；充分使用库内日志和结果，绘制发表级图表并详细解释。

## 开始前状态

main @ 0f20c3a。README顶部积累了2026-11-02之后的逐日状态通告，其下为2026-11-02论文式首页；`docs/homepage/0[1-6]_*.md`为对应分页。

## 做了什么

- 原README与README.en.md原样复制为`README_20261102_0906.md`、`README_20261102_0906.en.md`；旧分页未改动。
- 新README（中文）按论文结构重写，主线为0906 K9；新增第9节将V4 K14、WEHI R4 K15的读出并列展示，第10节为开放问题与分页索引。
- 新分页 `docs/homepage/v3/01–06`，每页链接上一版分页。
- 新图14幅（主图1–8、补图S1–S6），`scripts/homepage_v3/`渲染，PNG+PDF与SHA见`assets/publication/20261126/figure_manifest.json`；作图所用聚合值另存`inventories/publication/2026-11-26/homepage_v3/`。
- 图3c–e的UMAP散点复用2026-11-02渲染栅格（细胞坐标仅在服务器）。

## 关键发现 / 纠正

- `runs/20261106_CM07_DDR1_exploration_v1.md`记录：CosMx T3b“CD8距离”坐标单位为像素而非μm（推断约0.18 μm/px，待核），差异主要来自密度；旧首页未更正。新图6c改为展示观测比值与保持密度零模型期望（RCC3 4.17 vs 3.48；RCC4 1.04 vs 1.01；RCC5 5.54 vs 5.13），正文说明。
- 图4b“未命名权重”口径为1−named fine state（含恶性分层节点），因此CM07为78%（旧图按review_flag为66%）；图注写明口径。

## 验证

- 所有新Markdown相对链接存在性检查通过；`scripts/validate_repository.sh`：0 error、1项既有警告（2026-10-13 JSON路径）。
- 每幅图渲染后逐一目视检查并修正重叠；数值与源TSV/README交叉核对（scIB、Cox HR、ICI δ/q、Visium z、V4 K14种子/LOSO、R4偏相关）。

## 未验证 / 限制

- 无法访问研究服务器，未检索服务器上的稿件图源；所有图均来自仓库内汇总表。
- 图8f读出汇总的四档分类为对已记录结果的人工归纳，不是新统计。
- README.en.md未同步重写（仍为2026-11-02英文版）。

## 建议的commit message

Rewrite homepage as paper-style v3 with new figures and subpages
