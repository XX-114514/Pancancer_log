# 巨噬—T共同程序/自身表达敏感性 v1

方案冻结，真实运行未开始。上轮双向合并有8主BH/2亚型BH关联，原始研究核对发现部分PPI并非胞外通讯，热休克共同应激是需检验假设。原全局阴性及上轮探索结果保持。

本轮对同12原发情境、全部4624LR、全部兼容TF、每变体全9248分母运行，绝非只重跑8阳性。原丰度/细胞数/复杂度/受体、亚型敏感性、源内标准化、重复患者共块和source/LOSO不改；再加入发送者HSF1连续程序（原ULM排除全部LR基因，五源62—66实测targets）及接收者自身同一配体RNA（复合物取组分最小meanlog）。HSF1是一个程序，不能代表所有状态；不会赋予正式细胞状态标签。控制潜在真实中介也可能导致衰减，所以结果不能单独判为技术伪影。

[冻结参数/输入与代码SHA](../inventories/analysis/2026-10-12/macrophage_T_state_sensitivity_plan_v1.json)。沿用scanpy环境Python、BLAS/OMP2线程，环境PANCANCER_PREFLIGHT_ROOT/ANNOTATION_ROOT/CR2_CORE/FROZEN_LR_RESOURCE；PANCANCER_MAC_T_STATE_PLAN指向新plan；分别执行`python -u ${CODE_ROOT}/macrophage_T_state_sensitivity_v1.py`和`python -u ${CODE_ROOT}/audit_macrophage_T_state_sensitivity_v1.py`。新输出及execution manifest，不改旧registry、权威Zarr/172标签/9CM。

验收：现有缓存SHA、源unit/sample/major键对齐、完整family聚合和24均衡pooled/source/LOSO数值复算；仅此fixture数值复演，不声称全系数独立审计。真实状态尚待，事件UTC/UTC+8：2026-10-12T20:59:37.779726+08:00 /UTC 2026-10-12T12:59:37.779726+00:00

## 已完成

生产/审计均exit0；[真实结果/限制/后续](../inventories/analysis/2026-10-12/macrophage_T_state_sensitivity_v1/RESULTS.md)。尚未独立生物验证。
