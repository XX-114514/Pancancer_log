# WEHI数据服务器建立、空间数据下载与规划

## 用户要求

用户为本仓库配置了可写deploy key，要求：先处理好环境和云端工作流，然后下载需要的数据，汇总现有研究数据，最后规划后续任务。
此前同一会话里已经在WEHI上搭建了Xenium环境、入库7个来源（14单元），并按用户确认下载了GSE300147。

## 时间

- 任务开始：2026-10-30约01:35 UTC / 09:35 UTC+8（deploy key认证成功，克隆仓库）。
- 本日志取证：2026-10-30 02:33 UTC / 10:33 UTC+8 / 13:33 Australia/Melbourne。
- 更早阶段（Xenium 14单元入库、S18下载与入库）在仓库克隆之前完成，这里是补录，没有对应的更早提交。

## 开始前状态

- 仓库：`main@3dcece6`，工作树干净。WEHI上此前没有仓库检出，也没有Git凭据。
- 服务器：`~/spatial_xenium`下有一套不在仓库中的Xenium代码，Xenium 31单元已入库（含S18）。
- 读取：AGENTS、STATUS顶部、ROADMAP/TODO标题、[Cloud开发指南](../../../docs/CODEX_CLOUD_DEVELOPMENT.md)、[RUNTIME_SYNC](../../../docs/RUNTIME_SYNC.md)、
  [下一阶段方案](../../../methods/CM-project-next-stage-20261030.md)、[项目同步日志](2026-10-30_project-sync-and-next-stage-review.md)、目录约定。

## 计划

1. 验证deploy key，克隆仓库，跑基线检查。
2. 把WEHI代码移入仓库，去掉主机路径，建立仓库检查环境。
3. 按优先级下载空间来源和单细胞参考，做只读检查。
4. 写数据汇总、运行记录和规划，更新状态与索引，普通push并核对远端SHA。

## 执行与发现

**认证与Git。** 通过仓库级Host别名进行SSH认证，返回本仓库名。服务器git为2.47，直接使用ssh配置，不需要旧服务器那种包装器。
没有全局git身份；在仓库本地设置了与历史作者一致的身份。

**仓库检查环境（已观察）。** 系统python3为3.9，低于要求。用空间环境运行`cloud_check.sh`：618个测试，22个报错（zarr 3、numpy 2 StringDType、缺harmonypy/celltypist）。
只装标准库的3.11：77个导入报错，说明“仅依赖标准库”的说法已过时。按既往收据钉版本建了`${PC_REPO_ENV}`：在2 CPU上剩6个报错，其中
psutil缺失和harmonypy缺少`ncores`参数通过安装psutil、harmonypy 2.1.0解决；numba线程数报错只在2 CPU时出现。改用4 CPU批作业后，620个测试只剩1个报错：需要`ATLAS_SYMPHONY_WHEEL`，
即旧服务器的私有wheel。其余验证全部通过。

**代码移入仓库。** 新建`scripts/wehi_xenium/`：`config.sh`从未跟踪的本地文件读取路径；作业脚本从运行目录提交，并在第一行记录repo commit；
冻结面板直接读取仓库内的JSON（与服务器副本SHA相同，前缀`5c45c9e4`）。环境锁定文件去掉了本地prefix和conda构建路径。
仓库版样本表生成脚本重建的`samples.tsv`与入库所用的样本表逐字节相同。

**下载。** `make_geo_manifest.py`从GEO filelist生成清单，排除转录本级文件和图像；Zenodo IMC排除了原始TXT（7.9 GB）。GEO部分415个文件全部完成，0失败；
Zenodo的Steinbock归档仍在下载。

**数据检查的关键发现。**
- S16 Xenium、S20 CosMx和S22 MERFISH来自同一批TMA的连续切片（PMID 41006245），不是独立患者，但可用作平台校准。
- GSE328481是Stereo-XCR-seq：11位LUAD患者各1张切片，全转录组，带细胞边界和TLS分群。此前服务器笔记推测为MERFISH，是错误的。
- GSE273952（ccRCC Visium）和GSE316195（PDAC snRNA）是四队列中的两个独立来源；服务器笔记写成“配对”是错误的，已在汇总中纠正。
- 基因覆盖：CosMx缺LEF1、IRF1、TRAC、IKZF2；MERFISH缺6个干扰素基因中的5个和多数增殖基因；S18缺全部干扰素基因。

## 修改的文件

- `scripts/wehi_xenium/**`：WEHI代码、环境定义、清单（新增）。
- `tests/test_wehi_xenium_metadata_v1.py`：CRLF、patient/run解析和重复切片归并的回归测试（新增）。
- `docs/WEHI_SERVER_WORKFLOW.md`：服务器角色、环境、Git同步与运行约定（新增）。
- `methods/WEHI-spatial-line-plan-20261030.md`：WEHI空间线规划（新增）。
- `inventories/runtime/2026-10-30/wehi_spatial/DATA_SUMMARY.md`：来源级数据汇总（新增）。
- `runs/20261030_WEHI_spatial_ingest_and_download_v1.md`、`inventories/runs.tsv`、`inventories/datasets.tsv`、STATUS、TODO、ROADMAP、CHANGELOG、daily：状态与索引。

## 验证

见[运行记录](../../../runs/20261030_WEHI_spatial_ingest_and_download_v1.md#验证)。提交前运行`validate_repository.sh`并查看diff。

## 失败尝试

- S18第一次入库提交在样本表生成失败后仍被串联提交；在运行前已全部取消，没有产生输出。
- 第一次来源检查：tarfile读不了GSE176021，CosMx带引号的表头未正确解析；修正后重跑。
- 为等待作业，重复启动过等待命令，已停止重复的那个；不影响作业。

## 未验证

GSE328481计数是否为原始计数；三份同TMA数据的core→患者映射；IMC归档内容；新平台还没有入库成h5ad。

## 未解决

- S16用哪种分割作主分析（用户决定）。
- 旧服务器到WEHI的迁移通道：冻结T子集、私有登记、Symphony wheel。
- SCP2510需要用户登录下载。

## 下一步

按[WEHI空间线规划](../../../methods/WEHI-spatial-line-plan-20261030.md)：先做W1新平台读取器，再做W2/W3单位登记和测量mask；IMC下载完成后做W5输入复核。

## 建议的commit message

Record WEHI server setup, spatial downloads, data summary and plan
