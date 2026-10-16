# DIALOGUE完整聚合证据副本

[独立本机验收](independent_public_acceptance_v2.json)：公开文件SHA与完整输入验收链、14项BH、2385个带方向基因重叠独立复算通过。完整53个模型分数与官方值重建相符（此处未另读私有逐细胞原矩阵）。

- [14分层拟合](state_after_DIALOGUE_models_20261016_v2/public_fit_status.tsv)。
- [全部56端点槽状态](state_after_DIALOGUE_frozen_projection_20261016_v1/public_slot_scope.tsv)和[53槽重建](state_after_DIALOGUE_frozen_projection_20261016_v1/public_reconstruction_checks.tsv)。
- [全部28槽/14可估计内部留出](state_after_DIALOGUE_internal_holdout_20261016_v1/public_internal_holdout_results.tsv)，2BH记录；缺失保留。
- [全部MCP基因权重](state_after_DIALOGUE_frozen_projection_20261016_v1/public_fixed_gene_weights.tsv)。
- [全参考基因对应](state_after_DIALOGUE_programme_bridge_20261016_v1/all_MCP_reference_gene_overlap.tsv)。
- [科学解释和边界](../state_after_DIALOGUE_internal_review_v1/RESULTS.md)。

全部为聚合模型/基因/效应，不含细胞/患者分数或原始矩阵。拟合初始inventory的heldout_validated=false是当时阶段状态，后续内部留出另有独立终态，不改旧清单。内部留出不等于独立外部验证，MCP不等于LR轴或正式状态。

## 从仓库独立复算

本机工作布局的v2审计保留；仓库路径可用v3重算公开SHA、14项BH、2385基因重叠及范围：

```bash
python3 scripts/accept_state_after_DIALOGUE_public_v3.py \
  inventories/analysis/2026-10-16/state_after_DIALOGUE_completed_public_v1 \
  inventories/analysis/2026-10-14/full_state_programme_dictionary_v4/frozen_programme_dictionary.json \
  inventories/analysis/2026-10-16/state_after_DIALOGUE_frozen_validation_v1/input_acceptance_inventory.json \
  /tmp/dialogue_public_recheck.json
```

需要NumPy/pandas，不需要原始矩阵或患者数据；输出拒绝覆盖。
