# V7 panel 语义审计证据

取证/复算：2026-11-05 04:15–04:28 UTC（12:15–12:28 UTC+8），随后完成文档与发布核验。历史科学产物时间由来源路径及历史日志定义，本轮时间不是历史 V7 运行时间。

本轮静态审查 162 个候选 panel 及 17 组父类正向/排除 panel，并复算 40 个历史 unresolved 群的 373 个展开候选、2,052 条基因证据。详见[方法报告](../../../../methods/V7-panel-semantics-audit-20261105.md)与[逐 panel 完整表](all_panel_reviews.md)。V7 冻结标签保持原样；审计意见不是已发布新注释。

## 文件说明

- `all_162_panel_reviews.tsv`：每个旧 panel 的全部 marker、canonical、建议角色、具体审查意见与历史第一候选出现范围。范围不是误标计数；同群跨阶段可能出现多行，不可求和。
- `all_17_major_panel_reviews.tsv`：全部父类正向/排除组合与局限。
- `candidate_reconstruction.tsv`：全部 373 个候选的旧分数、通过基因、排名；`minimal_gene_gate` 仅重现旧“≥2 基因及≥1 canonical”，不代表生物学通过。
- `cluster_rule_audit.tsv`：40 群审查；`shadowed_gene_supported_candidates.tsv` 为 9 群、323,669 个细胞的规则缺陷暴露范围，不代表全部可纠正。
- `within_lineage_panel_overlaps.tsv`：共享 marker 供人工检查；共享本身不是错误，也不是污染证明。
- `sources/`：历史代码与 cluster×panel×gene 聚合证据。**历史代码中的注释和科学解释属于历史表述，不代表本轮认可。**这些文件用于取证，禁止作为新注释流程直接重跑。
- `inputs/`：原先已核验的聚合 gate/crosswalk 与 359 个 marker 列名。无细胞条形码、患者行或矩阵。
- `source_manifest.json`：原始来源相对路径与 SHA256。两个 prior-stage 脚本复用上一轮保存的源码，其余本轮只读取回。原始未公开的表对应 SHA 仍为历史原文件，不要与显式 NA 的公开表混淆。

## 重现

在仓库根目录、安装 numpy/pandas 的 Python 环境中，选择尚不存在的新输出目录：

```bash
python scripts/v7_panel_semantics_audit_20261105/reconstruct_old_candidates.py \
  --input-dir inventories/analysis/2026-11-05/v7_panel_semantics_audit_v1 \
  --output-dir /tmp/v7_panel_audit_reproduction_new
python scripts/v7_panel_semantics_audit_20261105/render_panel_review.py \
  --input-dir inventories/analysis/2026-11-05/v7_panel_semantics_audit_v1 \
  --result-dir /tmp/v7_panel_audit_reproduction_new
```

重现覆盖保存的基因效应统计、通过门槛、候选分数及第一名；不是从 467 万细胞重新计算表达，也未重新计算所有候选的样本复现、模型预测或新的细胞注释。评论文本来自显式人工审查字典 `panel_review_notes.py`。

## 尚未解决

T/NK/ILC 父身份与分支重核、候选专属跨来源验证、功能/组织背景命名验证及新版本细胞级验收仍需后续完成。别名已有历史修复，不能将所有 unresolved 简化为基因缺失或分辨率不够。
