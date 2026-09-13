# V8 Myeloid/DC 全量 review-only run（2026-09-13）

## 状态

`running`。Run ID：`20260913_084700_mydc_full_v1`。V7 保持冻结默认；
`writeback=false`、`independent_benchmark=false`。

## 全量范围

- V7 major Myeloid/Macrophage/Monocyte/Dendritic/Granulocyte：570,957 cells。
- Myeloid/DC-related uncertainty audit：371,230 unique cells，0 unmatched。
- 去重并集：882,415 cells、1,321 samples、54 GSE。

审计行仅在 original/marker/CellTypist/conflict/candidate/negative-marker 字段明确涉及
Myeloid、Macrophage、Monocyte、Dendritic、Granulocyte、neutrophil 或 cDC/pDC/DC3 时进入。

## 方法边界

- 359-gene 输入 X 为 full-library denominator 的 `log1p(CP10K)`；结构性缺失按
  `availability_by_gse` 排除；marker-only counts layer 禁止用于 QC/count inference。
- 每样本独立 PCA/neighbors/Leiden/marker review；少于 50 cells 显式 skipped。
- 固定阈值：top1 score >=0.5、top1-top2 margin >=0.1、至少 20% cells 检出 >=2 panel
  genes、panel coverage >=50% 且 >=2 genes、contamination <20%。
- Identity、连续 state、contamination、margin 与 rejection reason 分字段保存；macrophage
  C1QC/SPP1/FOLR2 programs 进入 state，不替代 cell type。

## 启动与资源

- `setsid + flock`；PID 559，PPID=1，PGID=SID=559。
- 1 worker、8 numerical threads、nice 10、128 GiB virtual-memory soft limit。
- 08:51 快照：25 success、30 skipped、1 running；Python RSS 约 2.99 GiB。
- 首次 launcher 在任务已启动后因旧版 `ps -o` 语法误报；未重复启动，补记 PID identity
  manifest 后修复 launcher。

## 验证

- smoke v3：2/2 samples success。
- unit tests：4 passed。
- launcher contract：PASS。
- cell type 禁用恶性/CNV/可评估性状态检查：0 illegal terms。

## 监控

```bash
bash ${PROJECT_ROOT}/Project_v3/v8_upgrade/fullscale/myeloid_dc_20260913/scripts/monitor.sh \
  ${PROJECT_ROOT}/Project_v3/v8_upgrade/fullscale/myeloid_dc_20260913/runs/20260913_084700_mydc_full_v1
```

只有全部样本 terminal、0 failed 且 merge summary 通过，才可把计算状态改为完成；即使完成也
仍是 review-only recurrence evidence，不是 accuracy 或 V8 biological freeze。
