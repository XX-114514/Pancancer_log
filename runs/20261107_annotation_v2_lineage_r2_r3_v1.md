# 20261107 注释 v2：第二轮审阅、第三轮重跑、低深度转移、稳定性检查、多种子共识（v1）

方法与问题见 [methods/tme-annotation-v2-lineage-rounds-20261107.md](../methods/tme-annotation-v2-lineage-rounds-20261107.md)。本文件只记录运行。

- 状态：**进行中**。第二、三轮聚类、检查和三个池的转移已完成；多种子共识第一步在运行；同步到旧服务器在进行。
- 计算环境：WEHI Milton，`${PC_ROOT}/envs/conda/scanpy`。
- 输入：TME v1 run v2（`${PC_ROOT}/tme_v1/run_v2`：`tme_v1.zarr` 和 shards）、`annotation_v2/doublets`、合并面板 v2。

## 作业

| 作业 | 内容 | 代码 commit | 开始–结束（UTC） | 结果 |
|---|---|---|---|---|
| 30676198 | 应用第一轮判定 | 84c00103 | 11-06 22:19–22:21 | 完成 |
| 30676201–30676209 | 第二轮，9 个池 | 84c00103 | 11-06 22:21–11-07 05:42 | 全部完成；簇数见方法文件 §2.3 |
| — | 应用第二轮判定（会话内运行 `apply_lineage_r2_review.py`） | worktree 0764824a + 未提交文件 | 11-07 约 07:45 | `cell_labels_r3.tsv.gz`（SHA256 74b0973c…） |
| 30677658 | 第三轮 Myeloid，656,808 个细胞 | 0764824a | 07:48–10:53 | 22 簇，30 个 PC |
| 30677659 | 第三轮 Endothelial，137,014 个细胞 | 0764824a | 07:48–09:58 | 19 簇，30 个 PC |
| 30677660 | LISI、Harmony 对照、高分辨率基因检出 | 0764824a + 未提交脚本 | 07:49–09:52 | `r3_checks/lisi_harmony.tsv` 等 |
| 30677663 | 标签转移，第一次 | 同上 | 07:51，已取消 | 读错标签列的 bug，结果作废（方法文件 §5 问题 14） |
| 30677673 | 标签转移：Stromal、B_Plasma、T_NK | 同上，已修复 | 07:52–08:05 | `r3_checks/transfer_*` |
| 30678135 | Leiden 种子和 Harmony 稳定性对照：Endothelial、Stromal | 同上 | 10:00–10:14 | `r3_checks/leiden_stability.tsv` |
| 30678224–30678227、30678229 | 多种子共识第一步：T_NK、B_Plasma、Stromal、Endothelial（第三轮）、Myeloid（第三轮） | 同上 | 11:23 起 | 运行中 |

上表中“未提交”的文件在本次提交中入库（见任务日志）。

## 输出（不入库；以 `${PC_ROOT}/tme_v1/annotation_v2` 为根）

`lineage_r2/`、`lineage_r3/{Myeloid,Endothelial,checks,transfer,consensus}/`、`lineage_r3/cell_labels_r3.tsv.gz`。汇总表复制在 [inventories/annotation_v2/r3_checks/](../inventories/annotation_v2/r3_checks/)，附 SHA256SUMS。

## 同步到旧服务器（用户要求，2026-11-07）

- **方式：** 走 R2 中转，复用 [迁移工具](../scripts/wehi_migration/README.md) 反向传输：WEHI 端用 `old_pack_upload.sh` 打包上传，旧服务器端用 `wehi_fetch_verify.sh` 下载校验。
- **约束：** 同一时间只有一个上传方，桶上限 9G。每块核对 SHA256 后才从桶里删除；解包后逐文件核对，并写收据。
- **目标：** 旧服务器的新目录 `from_WEHI_20261107/`，不覆盖已有数据。
- **不传：** 用户决定（2026-11-07）不传 `object.h5ad`；tme_v1 shards、integration_v2 的 decontx/harmonised_new、ingest 也不在本批。

| 项 | 内容 | 大小 | 状态（11-07 11:26 UTC） |
|---|---|---|---|
| w2o_records | WEHI 独有的 4 个运行、annotation_v2 的 major 和第一轮结果表、doublets、malignancy_v1、整合运行的报告和指纹、investigations | 13G，2,451 个文件 | 已验证，有收据 |
| w2o_annotation_r2 | 第二轮结果表和 `cell_labels_r2` | 107M | 已验证，有收据 |
| w2o_tme_v1_zarr | `tme_v1.zarr` | 64G，215,753 个文件 | 已验证（64/64 块，11:16 UTC） |
| w2o_analysis_ready_zarr | `pancancer_analysis_ready_v1.zarr` | 19G | 上传中 |
| w2o_integration_v2_zarr | `integration_v2.zarr` | 107G | 排队 |

旧服务器的接收速度在约 3.4–30 MB/s 之间波动。桶中原有的 `XX/` 前缀（846 个对象，414 MiB）已按用户指示删除。

第三轮结果和最终标签定稿后，另作一项同步。
