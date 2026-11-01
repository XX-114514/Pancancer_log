# Run：整合数据集 v2（2026-11-01/02）

说明见[docs/INTEGRATION_V2.md](../docs/INTEGRATION_V2.md)。环境：WEHI `${PC_ROOT}/envs/conda/scanpy`（scanpy 1.11.5、anndata 0.11.4、harmonypy 2.0.0）；decontX使用`${PC_ROOT}/envs/conda/scrna_r_full`中的celda 1.18.2。

| 步骤 | 作业 | 结果 |
|---|---|---|
| decontX，原51个来源 | 30625727（数组任务，另有测试作业30625589） | 51/51完成，没有单位失败 |
| v2下载/读取/QC | 30626034/30626678、30626731–30626735、30626948、30626949、30628397 | 14个来源；GSE244983因基因数不足排除 |
| v3下载/读取/QC | 30627097、30628425–30628427、30628438、30628722 | 7个来源；GSE117570（基因数不足）、GSE190888（2名患者）、GSE180286（深度0.27倍）排除 |
| 新来源统一字段 | 30628724及v3对应作业 | 21个来源 |
| decontX，新来源 | 30628745、30629143 | 21/21完成，没有单位失败 |
| 整合 | 30629654（bigmem lrg-n01，7 h，峰值66 GB） | 4,720,931个细胞；Harmony 5轮收敛 |
| 报告 + marker图 | 30634868 | 回读通过 |
| 指纹查重 | fp-v2 | 完整指纹跨来源相同为0 |

输出：`${PC_ROOT}/integration_v2/run_v2/`（`integration_v2.zarr`、`figures/`、`report/`、`final_report/`、`fingerprint/`）。
