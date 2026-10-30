# WEHI空间数据下载、Xenium入库与来源检查（2026-10-30）

第一次在WEHI数据服务器上正式运行。全部步骤只是数据准备和描述，不定义状态，也不计算效应。汇总见[WEHI数据汇总](../inventories/runtime/2026-10-30/wehi_spatial/DATA_SUMMARY.md)，
流程见[服务器工作流](../docs/WEHI_SERVER_WORKFLOW.md)，任务日志见[codex日志](../logs/codex/2026/2026-10-30_WEHI-server-setup-and-spatial-data.md)。

## 环境

- `${PC_SPATIAL_ENV}`：Python 3.11.18、scanpy 1.11.x、anndata、squidpy 1.8.2、spatialdata 0.7.3、spatialdata-io 0.7.1；版本清单见`scripts/wehi_xenium/envs/pc-spatial.versions.txt`。
- `${PC_REPO_ENV}`：Python 3.11.18，按[pc-repo.yml](../scripts/wehi_xenium/envs/pc-repo.yml)钉版本，只用于仓库检查。
- 调度器：Slurm；下载用datamover分区，计算用regular分区。作业号保存在服务器的运行目录中，这里只记录阶段和状态。

## 阶段

| 阶段 | 代码 | 输入 | 输出 | 状态 |
|---|---|---|---|---|
| Xenium 7来源下载与入库（14单元） | 移入仓库前的服务器代码 | `manifests/xenium_tables_v1.tsv`中除S18以外的行 | 14个counts h5ad＋QC；面板报告 | 2026-10-29完成；仓库版样本表生成脚本逐字节复现所用样本表 |
| S18 GSE300147下载 | 同上 | 34个文件（17切片的矩阵和cells） | `raw/S18_GSE300147_HNSCC` | 完成，0失败，字节数全部吻合 |
| S18入库 | 同上；`--array=15-31` | 样本表第15–31行 | 17个h5ad＋QC；面板报告更新为31单元 | 17/17完成；面板汇总完成 |
| 新来源下载（GEO） | `b5eae67` `jobs/download.sbatch` | `manifests/spatial_sc_sources_v1.tsv`中GEO部分，415个文件，约12.2 GB | `raw/{S20..S23,V01,R01..R03}` | 完成，27分20秒，415/415，0失败 |
| 新来源下载（Zenodo IMC） | `b5eae67` | README、作者RDS、Steinbock_Result.zip（16.2 GB） | `raw/I01_Zenodo18265339_IMC_NSCLC_LN` | **运行中**：README与RDS完成；Steinbock下载限速约1.1–1.5 MB/s |
| 来源检查 | `7891ca3` `inspect_sources.py` | 所有非Xenium来源 | `processed/summary/sources/{source_sections.tsv,source_gene_coverage.tsv,SOURCES.md}` | 完成，1分31秒；282个条目来自9个来源；IMC归档未完成，标记为不可读 |

## 失败与修正（保留）

1. **S18第一次入库提交错误。** GEO的文本以CRLF结尾，`\r`进入了GSM元数据表，导致样本表生成失败；但我在同一条命令里串联了作业提交，入库array按旧的14行样本表提交。
   在任何任务开始前已全部取消（sacct显示CANCELLED，`samples.tsv`没有改变）。去掉`\r`后重建样本表，确认第1–14行不变，再重新提交。
   修正后的`fetch_geo_samples.py`逐行去掉`\r`，并加了回归测试（`tests/test_wehi_xenium_metadata_v1.py`）。
2. **第一次来源检查结果不完整。** Python的tarfile无法读取GSE176021的归档（归档里有macOS的`._`条目；`gzip -t`和GNU tar都能正常读取）；
   CosMx表头的引号使对照探针和坐标列识别失败。改为用csv解析表头、用GNU tar流式读取后重跑，结果已与仓库此前记录的数量级一致。
3. 第一次在空间环境里运行仓库检查时有22个报错，原因是环境版本，不是代码问题。处理办法见工作流文档。

## 验证

- 下载：每个文件的字节数都与GEO filelist或Zenodo API一致；每个来源目录都生成了`SHA256SUMS`。
- 入库：31个单元的坐标全部完整（缺失0），S18的patient和run写入了obs和uns。
- 来源检查：GSE176021和GSE316195的barcode数与仓库此前记录的QC后数量级一致（902,463 vs 868k；83,651 vs 82,983）。
- 仓库检查：在`${PC_REPO_ENV}`和4 CPU节点上运行，620个测试有1个报错（`ATLAS_SYMPHONY_WHEEL`未设置，属已知环境缺口）；新增3个测试通过；仓库校验0错误。

## 未验证

- GSE328481的X是否为原始计数（目前只检查了前200行是否为整数值）。
- S16、S20、S22的core→患者映射。
- IMC的Steinbock归档内容（仍在下载）。
