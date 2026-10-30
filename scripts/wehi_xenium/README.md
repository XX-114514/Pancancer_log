# WEHI 空间与单细胞数据线（`scripts/wehi_xenium/`）

WEHI 数据服务器上的下载、Xenium 入库和面板覆盖代码。服务器从本仓库检出明确 commit 运行；
真实路径只写在服务器本地未跟踪的配置文件里，本目录不含主机路径、样本级表或数据。

## 角色与环境

| 名称 | 作用 | 定义 |
|---|---|---|
| `${PC_SPATIAL_ENV}` | 空间分析（scanpy 1.11、squidpy、spatialdata 0.7，zarr 3 / numpy 2） | [envs/pc-spatial.yml](envs/pc-spatial.yml)，版本见 `pc-spatial.versions.txt` |
| `${PC_REPO_ENV}` | 仓库检查 `scripts/cloud_check.sh`（按既往收据钉 numpy 1.26 / zarr 2.18 / anndata 0.11） | [envs/pc-repo.yml](envs/pc-repo.yml) |
| `${XENIUM_DATA_ROOT}` | `raw/`（只读下载 + 每来源 SHA256SUMS）、`processed/{h5ad,qc,summary}` | scratch，不备份，30 天未访问清理 |
| `${XENIUM_META_ROOT}` | `samples.tsv`、`geo/samples_GSE*.tsv`（GEO 公开样本信息） | home，持久 |

本地配置文件默认是 `~/.config/pancancer/wehi_xenium.local.sh`（可用 `XENIUM_LOCAL_CONFIG` 改），需定义
`XENIUM_DATA_ROOT`、`XENIUM_META_ROOT`、`PC_SPATIAL_ENV`、`PC_REPO_ENV`。[config.sh](config.sh) 读取它并导出派生路径；
冻结面板直接读取仓库内 `inventories/analysis/2026-10-18/priority_immune_CM_bulk_v1/frozen_plan_before_effects.json`。

## 运行

所有作业从日志目录提交，并导出 `XENIUM_CODE_ROOT`（sbatch 默认传递环境变量）：

```bash
export XENIUM_CODE_ROOT=${REPO}/scripts/wehi_xenium
source "$XENIUM_CODE_ROOT/config.sh"; cd "${LOG_DIR}"
sbatch "$XENIUM_CODE_ROOT/jobs/download.sbatch" "$XENIUM_CODE_ROOT/manifests/xenium_tables_v1.tsv"
python3 -I "$XENIUM_CODE_ROOT/fetch_geo_samples.py" GSE300147 "$META_ROOT/geo/samples_GSE300147.tsv"
python3 -I "$XENIUM_CODE_ROOT/make_samples.py" "$RAW_ROOT" "$SAMPLES_TSV"
n=$(($(wc -l < "$SAMPLES_TSV") - 1))
ingest=$(sbatch --parsable --array=1-$n "$XENIUM_CODE_ROOT/jobs/01_ingest.sbatch")
sbatch --dependency=afterok:$ingest "$XENIUM_CODE_ROOT/jobs/02_panel_coverage.sbatch"
```

仓库检查（测试中有多线程用例，请在 ≥4 CPU 的计算节点上运行）：

```bash
PATH=${PC_REPO_ENV}/bin:$PATH PYTHON=python bash scripts/cloud_check.sh
```

## 清单

- [manifests/xenium_tables_v1.tsv](manifests/xenium_tables_v1.tsv)：Xenium 处理后表格（矩阵、cells、analysis、面板、作者标签）。不含图像、转录本和 FASTQ。
- [manifests/spatial_sc_sources_v1.tsv](manifests/spatial_sc_sources_v1.tsv)：CosMx/MERFISH/空间 h5ad/Visium/IMC 和四队列 sc/snRNA。
  由 [make_geo_manifest.py](make_geo_manifest.py) 从 GEO filelist 生成；转录本级文件和 Zenodo 原始 IMC TXT 已排除。

## 口径

- `ingest_xenium.py`：X 为未修改的基因计数；`qc_pass` = 计数 ≥10，只标记不删除；patient/run 只在来源公开时填写。
- S07 两行是技术重复切片；S16 每张 TMA 有 unimodal 和 multimodal 两种分割，core→患者映射建立前不能作患者单位；
  S18（GSE300147）7 位患者有两次 run，统计单位是患者。
- `panel_coverage.py` 只做描述：面板覆盖、检出和 QC。不定义状态，不计算效应。
