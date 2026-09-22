# V9 candidate quarantine

V9 是截至 2026-09-22 分析工作区中最新的全量注释尝试，但不是投稿权威版本。

## 已观察状态

- 外部 manifest 状态：`FROZEN_CANDIDATE`，不是 `FROZEN`。
- 对象范围：4,676,787 cells、1,322 samples、54 GSE、43 cancers。
- L3 只完成 B/Plasma 和 Myeloid 试点；15 个 cluster 等待人工书面复核。
- Myeloid L3 的所有分辨率 bootstrap ARI 均未达到预注册门槛；多个 L2/L3
  双种子稳定性门槛未通过。
- 恶性轴没有重跑 CNV：3,353,123 个细胞继承 V7 CNV；632 个样本被建议重跑，
  但尚未执行。
- V9 manifest 含内部绝对路径，因此只记录脱敏逻辑路径和哈希，不直接复制文件。

## 证据锁

| Artifact | Bytes | SHA-256 |
| --- | ---: | --- |
| `${PROJECT_ROOT}/Project_v3/v9_annotation/V9_FINAL_REPORT_20260919_CN.md` | 11,835 | `f19c9e94482c67aa40cbe35d36f8fa31ed97aa424301021c13118cad31cc2ae0` |
| `${PROJECT_ROOT}/Project_v3/v9_annotation/preregistration/V9_FREEZE_MANIFEST_v2.json` | 203,565 | `af83bec2c4dc11cdc96e07182a72f3da37bc2febe3f16825f93bd21b50406ac6` |
| `${PROJECT_ROOT}/Project_v3/v9_annotation/preregistration/manual_review_decisions_L3_pilot_TEMPLATE_v2.tsv` | 2,962 | `0fc20096f4360ee5a0a10cbe663bc13f545409126f9960e1adc6114adb7c4475` |

## Publication rule

V9 may be discussed as a candidate/method-development analysis only. It must not
replace V7 in the main denominators, title, summary paragraph, figures or Data
Availability statement until the manual-review queue, CNV rerun decision,
stability gates, downstream impact audit and explicit promotion are closed in a
new append-only release.
