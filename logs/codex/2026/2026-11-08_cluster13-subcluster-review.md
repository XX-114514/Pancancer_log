# 簇13子聚类 job 30667416 人工复核 — 2026-11-08

取证时间：2026-11-08T07:44:39Z / 15:44:39 UTC+8。

## 用户要求
检查 job 30667416；完成则读取 subclusters.tsv 和图、按 v1 方案人工复核、起草约6.8万个非恶性细胞的标签、列出疑似残留肿瘤子簇、复制聚合表、写日志、提交推送。

## 开始前状态
job COMPLETED（exit 0，41:47，MaxRSS 17 GB）。本地 main 落后 origin；存在未跟踪的 scripts/nmf_benchmark，未改动，提交经临时 worktree 完成。origin/main 的 v2 方案（2026-11-06）已将本运行定为诊断，v1 自动标签作废，细胞回到各癌种分析。

## 执行
读取 methods/tme-de-novo-annotation-v1-20261105.md、v2 方案对应章节、subcluster_v1.py 判定代码；查看 umap_subcluster、umap_auto_label、umap_tme_flag、dotplot_panel；用 de_top50 和 gene_dot 核对争议子簇；由 scripts/annotation_v1/cluster13_manual_labels_draft_v1.py 生成草案。

## 已观察事实
- 74个子簇，290,546个细胞，非恶性68,374个，与 v2 记录一致。脚本断言总数通过。
- 子簇33、18 中 CSF3R、FCGR3B、CXCR2 检出 ≤2%。子簇10 TP63 检出1%，SAA1、CCL28、LTF、MMP7 上调。
- subcluster_v1.py 未实现方案中双细胞共表达>30%的条件；residual_tumour_suspect 标记64/74个子簇。

## 解释
- 子簇33、18 不是中性粒双细胞；子簇10 更符合乳腺管腔祖细胞。与 v2 诊断文字不同，详见 inventories/annotation_v1/cluster13/README.md。

## 修改
新增 inventories/annotation_v1/cluster13/（聚合表、manual_labels_draft.tsv、residual_tumour_suspects_for_user.tsv、README）、脚本、本日志和 daily 日志。v2 方案和恶性判定未改。

## 验证
总数断言通过；上传文件为聚合表，无细胞级行。

## 未验证与待决定
双细胞未做细胞级共表达检验；未做子簇与 v2 标签交叉表；残留肿瘤子簇待用户决定（多数已被 2026-11-06 决定排除）；未解决：是否修订 v2 文字。

## 建议 commit message
Cluster 13 subcluster diagnostic: manual label draft and residual-tumour review

## 文字更正（2026-11-08，用户要求）
更新 methods/tme-annotation-v2-lodi-workflow-20261105.md 的簇13诊断段：子簇10改为管腔祖细胞样，子簇33、18 不再列为中性粒双细胞，计数按复核分类修正，并加更正说明保留原文措辞。未改细胞标签、恶性判定或 CoVarNet v1 范围。
