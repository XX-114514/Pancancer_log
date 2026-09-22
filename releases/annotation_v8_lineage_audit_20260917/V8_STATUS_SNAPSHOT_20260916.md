# V8 状态总览（自动生成）

> 由 `freezes/build_freeze_manifests.py` 于 2026-09-16 02:26 生成，请勿手工编辑。
> 详细说明见 `reports/13_V8_UPGRADE_AUDIT_AND_REFACTOR_20260916_CN.md`。

- 全局 V8：**NOT_FROZEN**；V7 仍是默认参考与正式下游输入。
- 任何谱系均未完成 P1-2 held-out benchmark，因此没有 `HELDOUT_VALIDATED`。

| 谱系 | 证据等级 | 决策 GSE 数 | 决策 GSE | held-out 暴露 | 已知问题数 |
|---|---|---:|---|---|---:|
| Myeloid_DC | `PILOT_FROZEN` | 3 | GSE131907, GSE161529, GSE274229 | HELDOUT_OUTPUTS_BEFORE_SPLIT | 5 |
| T_NK | `PILOT_FROZEN_HELDOUT_COMPROMISED` | 4 | GSE116256, GSE131907, GSE188737, GSE274229 | HELDOUT_USED_IN_DECISION, HELDOUT_OUTPUTS_BEFORE_SPLIT | 5 |
| Fibroblast | `PILOT_FROZEN` | 3 | GSE131907, GSE184880, GSE188737 | CLEAN | 5 |
| Endothelial | `PILOT_FROZEN` | 3 | GSE222315, GSE256136, GSE274229 | CLEAN | 6 |
| B_Plasma | `DEVELOPMENT_FROZEN` | 39 | GSE115978, GSE127465, GSE131907, GSE132465, GSE138709, GSE148071, GSE152048, GSE155698, GSE159115, GSE161529, GSE162025, GSE163558, GSE166555, GSE169609, GSE176031, GSE178341, GSE178481, GSE179784, GSE183904, GSE184362, GSE184880, GSE188711, GSE188737, GSE189357, GSE189889, GSE197461, GSE200997, GSE222315, GSE244101, GSE256136, GSE261774, GSE270680, GSE273127, GSE274229, GSE274934, GSE288479, GSE290925, GSE299340, GSE302903 | HELDOUT_OUTPUTS_SEALED_BY_POLICY | 5 |

## 已知问题

### Myeloid_DC
- **V8-KI-MY-001**（major）：The 54-GSE full-cohort review (finished 2026-09-13 11:48) produced outputs for all 12 held-out GSEs before the held-out split was frozen (2026-09-14 22:06). Not cited in the decision, but those GSEs are unblinded for Myeloid.
- **V8-KI-MY-002**（moderate）：FCN1_monocyte was accepted although its third-GSE top1 clusters failed the frozen contamination gate (validation caveat recorded in the decision table).
- **V8-KI-001**（moderate）：Review clusters are built on the 359-gene marker object that contains every scoring panel (partial circularity). Only Fibroblast/Endothelial have been sensitivity-tested; T/NK, Myeloid, B/Plasma not yet.
- **V8-KI-002**（moderate）：Identity is decided per cluster, so every cell in a cluster inherits the call; within-cluster heterogeneity was not recorded by the historical runners (lineage_review now writes cluster_core_dispersion).
- **V8-KI-003**（minor）：Historical runner scripts were not hashed into run provenance; lineage_review now records engine and spec SHA-256.

### T_NK
- **V8-KI-TNK-001**（critical）：Held-out GSE116256 was the preregistered T/NK validation cohort (2 samples, 3,377 cells) and contributes to the frozen T_or_NK boundary recurrence ('3 GSE'). T/NK cannot be evaluated as held-out on GSE116256.
- **V8-KI-TNK-002**（minor）：The frozen tissue_resident_like state recurred in GSE131907 and GSE188737 only within the primary run (GSE116256 did not contribute); the '3 GSE' claim adds GSE274229 from the classifier-extension run.
- **V8-KI-001**（moderate）：Review clusters are built on the 359-gene marker object that contains every scoring panel (partial circularity). Only Fibroblast/Endothelial have been sensitivity-tested; T/NK, Myeloid, B/Plasma not yet.
- **V8-KI-002**（moderate）：Identity is decided per cluster, so every cell in a cluster inherits the call; within-cluster heterogeneity was not recorded by the historical runners (lineage_review now writes cluster_core_dispersion).
- **V8-KI-003**（minor）：Historical runner scripts were not hashed into run provenance; lineage_review now records engine and spec SHA-256.

### Fibroblast
- **V8-KI-FB-001**（minor）：State recurrence in the frozen run counted all clusters, not only identity-positive clusters; the restriction was recomputed post hoc in EXECUTION_RECORD.md.
- **V8-KI-FB-002**（minor）：Clustering sensitivity: cell-level label agreement 92-97% and myCAF recurrence retained in all variants, but cluster partitions are unstable (within-sample ARI 0.22-0.47); excluding scoring panels moves 414 cells from Fibroblast_parent to contamination fallback.
- **V8-KI-001**（moderate）：Review clusters are built on the 359-gene marker object that contains every scoring panel (partial circularity). Only Fibroblast/Endothelial have been sensitivity-tested; T/NK, Myeloid, B/Plasma not yet.
- **V8-KI-002**（moderate）：Identity is decided per cluster, so every cell in a cluster inherits the call; within-cluster heterogeneity was not recorded by the historical runners (lineage_review now writes cluster_core_dispersion).
- **V8-KI-003**（minor）：Historical runner scripts were not hashed into run provenance; lineage_review now records engine and spec SHA-256.

### Endothelial
- **V8-KI-EN-001**（minor）：KDR is in both Endothelial_core and state_angiogenic_tip, so that state correlates with identity by construction (it fell back to parent).
- **V8-KI-EN-002**（minor）：Historical runner carried an unused audit matcher copied from Fibroblast (fibro/CAF/stromal terms); it had no effect (audit_unique_n=0) and is not reproduced in lineage_review.
- **V8-KI-EN-003**（minor）：Clustering sensitivity: cell-level label agreement 92.6-96.4% and arterial/venous stay cross-GSE recurrent in all variants, but support shrinks (venous 3 GSE/5 samples/1,718 cells -> 2 GSE/2 samples/841 cells at Leiden 0.3, exactly at the >=2/>=2 threshold); sample GSE256136|GSM8086247_N10 drops to 50.7% label agreement.
- **V8-KI-001**（moderate）：Review clusters are built on the 359-gene marker object that contains every scoring panel (partial circularity). Only Fibroblast/Endothelial have been sensitivity-tested; T/NK, Myeloid, B/Plasma not yet.
- **V8-KI-002**（moderate）：Identity is decided per cluster, so every cell in a cluster inherits the call; within-cluster heterogeneity was not recorded by the historical runners (lineage_review now writes cluster_core_dispersion).
- **V8-KI-003**（minor）：Historical runner scripts were not hashed into run provenance; lineage_review now records engine and spec SHA-256.

### B_Plasma
- **V8-KI-BP-001**（major）：45.4% of development cells (252,944) fall back to competing lineage. Cell-level: 29.6% carry no B/Plasma identity (likely V7 mislabel, e.g. Bcell_unresolved_CST3 / Bcell_coarse_provisional), 33.6% co-express identity and a foreign lineage (doublet/ambient), 28.0% are clean identity cells rejected with their cluster. Per-GSE fallback ranges 11%-100%.
- **V8-KI-BP-002**（moderate）：Held-out GSE predictions were computed and remain readable on disk; exclusion is enforced by code and policy, not access control.
- **V8-KI-001**（moderate）：Review clusters are built on the 359-gene marker object that contains every scoring panel (partial circularity). Only Fibroblast/Endothelial have been sensitivity-tested; T/NK, Myeloid, B/Plasma not yet.
- **V8-KI-002**（moderate）：Identity is decided per cluster, so every cell in a cluster inherits the call; within-cluster heterogeneity was not recorded by the historical runners (lineage_review now writes cluster_core_dispersion).
- **V8-KI-003**（minor）：Historical runner scripts were not hashed into run provenance; lineage_review now records engine and spec SHA-256.

## 证据等级定义

- `PILOT_FROZEN_HELDOUT_COMPROMISED`：决策使用了 held-out GSE，P1-2 对该谱系不能称为留出评估。
- `PILOT_FROZEN`：决策证据 ≤3 个 GSE。
- `DEVELOPMENT_FROZEN`：决策证据 >3 个非 held-out GSE。
- `HELDOUT_VALIDATED`：通过 P1-2 held-out benchmark（尚无）。
