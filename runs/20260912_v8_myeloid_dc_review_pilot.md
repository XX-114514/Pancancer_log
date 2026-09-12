# V8 Myeloid/DC review-only pilot（2026-09-12）

## 状态

`PASS_REVIEW_ONLY_CROSS_GSE_RECURRENCE`；非 held-out accuracy、非 V8 biological
freeze、未写回 V7、未修改正式 downstream。

## 输入与范围

- 冻结 V7 obs-only：4,676,787 cells，SHA-256
  `e5b68d7e5f127bcdcdadc3aaa94a9bce027ec032e6b82e68e5847ea5ccf7b38e`。
- 表达输入：GSE161529 与 GSE131907 counts 对象。
- 样本：GSE161529 两个样本、GSE131907 `EFFUSION_12` 与 `NS_19`。
- 候选宇宙：V7 Myeloid/Macrophage/Monocyte/Dendritic/Granulocyte 与审计 cell IDs 的并集；
  V5 fallback 禁用。

## 迭代与否决

- V1：marker-only broad filter 混入 epithelial、stromal 和 cycling clusters，否决。
- V2：发现错误使用 V5 uncertainty fallback，且一次 sample exception 被错误汇总为 success；否决。
- V3：直接使用 V7 obs-only，fail-closed 汇总和 post-hoc schema validator 通过。

## 结果

- 2 GSE、4 samples、22 review clusters、8,997 个 cluster-level cell counts。
- `candidate_monocyte_FCN1`：2 GSE、4 samples、13 clusters、7,377 cells。
- `candidate_macrophage_C1QC`：2 GSE、2 samples、3 clusters、1,012 cells。
- cDC2、pDC、SPP1 等仅在一个 GSE 出现，不能声称跨 GSE 复现。
- GSE131907 adapter 重跑的两张 label 表及合并表 SHA-256 与首次运行一致。

## Census/scANVI 状态

隔离环境固定 Python 3.11.14 与 `cellxgene-census==1.17.0`。默认镜像安装因依赖源构建失败；
唯一官方 PyPI binary-only retry 因 TLS `UNEXPECTED_EOF_WHILE_READING` 失败。状态为
`BLOCKED_NO_USABLE_PINNED_CLIENT_AFTER_BINARY_RETRY`。未修改项目 conda 环境、未读取 Census
expression、未训练 scANVI。

## 证据位置

- `${PROJECT_ROOT}/Project_v3/v8_upgrade/reports/12_V8_EXECUTION_PROGRESS_20260912_CN.md`
- `${PROJECT_ROOT}/Project_v3/v8_upgrade/pilots/unresolved_myeloid_dc_20260912/`
- `${PROJECT_ROOT}/Project_v3/v8_upgrade/pilots/unresolved_myeloid_dc_cross_gse_20260912/`
- `${PROJECT_ROOT}/Project_v3/v8_upgrade/reference/runtime_audit_20260912/`

## 下一门控

在查看第三 GSE 结果前预注册 cluster 接受/退回规则；随后验证尚未复现的 DC/macrophage
labels。完成独立 held-out benchmark 前不得将 recurrence 写作 accuracy 或写回 V7。
