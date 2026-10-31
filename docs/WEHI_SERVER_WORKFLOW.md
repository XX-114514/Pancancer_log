# WEHI数据服务器工作流

2026-10-30起，WEHI HPC作为新的数据服务器，承担空间与单细胞来源的下载、入库和面板资格审查。
旧数据服务器的私有登记、冻结图谱zarr和既往运行产物**不在**WEHI上。2026-10-30起WEHI可经SSH连到旧服务器，但实测只有约4–30 KB/s，只适合代码和小型登记；
大对象按[按需迁移规划](../methods/WEHI-data-migration-plan-20261030.md)经对象存储中转，公开数据在WEHI重新下载，并在日志中注明来源。本文件只写角色和流程，主机名、账户和路径保存在服务器本地配置。

## 角色

| 角色 | 内容 | 是否进入本仓库 |
|---|---|---|
| 代码 | 本仓库检出；WEHI专用代码在[scripts/wehi_xenium/](../scripts/wehi_xenium/README.md) | 是 |
| 本地配置 | `~/.config/pancancer/wehi_xenium.local.sh`：数据根、元数据根、两个环境前缀 | 否（权限600） |
| `${XENIUM_DATA_ROOT}` | 下载（只读）与处理输出，位于scratch：不备份，30天未访问会被清理 | 否 |
| `${XENIUM_META_ROOT}` | 样本表与GEO公开样本信息，位于home | 否 |
| 运行日志 | 每次运行单独目录，记录repo commit与Slurm作业号 | 否；摘要写入`runs/` |

## 环境

| 环境 | 用途 | 定义 |
|---|---|---|
| `${PC_SPATIAL_ENV}` | 空间分析：scanpy 1.11、squidpy、spatialdata 0.7（zarr 3 / numpy 2） | [pc-spatial.yml](../scripts/wehi_xenium/envs/pc-spatial.yml) |
| `${PC_REPO_ENV}` | 仓库检查：Python 3.11 + 按既往运行收据钉的numpy 1.26.4、pandas 2.2.3、zarr 2.18.3、anndata 0.11.4、scanpy 1.11.5、harmonypy 2.1.0 | [pc-repo.yml](../scripts/wehi_xenium/envs/pc-repo.yml) |

两套环境分开的原因（已观察）：在空间环境中运行`cloud_check.sh`，618个测试有22个报错，全部来自zarr 3的`LocalStore`接口、numpy 2的StringDType转换以及缺少harmonypy/celltypist；
只装标准库的Python 3.11则有77个测试在导入numpy/pandas时报错。因此[Cloud指南](CODEX_CLOUD_DEVELOPMENT.md)中“仅依赖标准库”的说法已不完整。
服务器系统自带的python3是3.9，低于检查脚本要求的3.10。

在`${PC_REPO_ENV}`和4 CPU计算节点上运行，620个测试只有1个报错：`test_single_source_symphony_replay_v1`需要`ATLAS_SYMPHONY_WHEEL`，它指向旧服务器上的私有wheel，WEHI上没有。
其余验证（运行快照、仓库校验、release校验、`git diff --check`）全部通过。多线程用例要求CPU≥4，因此要用`sbatch -c 4`运行，不要放在2 CPU的交互作业里。

## Git同步

- 认证：仓库级deploy key（可写），通过`~/.ssh/config`中的独立Host别名使用；私钥只在服务器上，不写进任何记录。
  服务器git 2.47可以直接使用ssh配置，不需要旧服务器那种`GIT_SSH`包装器。GitHub的22端口在这台服务器上可达。
- 提交身份沿用仓库历史作者，配置在仓库本地（`git config user.*`），不写入全局。
- 流程沿用[AGENTS](../AGENTS.md)和[RUNTIME_SYNC](RUNTIME_SYNC.md)：更新任务日志、daily、run和状态，运行验证并查看diff，显式stage，普通commit，普通push，再用`git ls-remote`独立核对远端main的SHA。
  禁止force push和改写历史。数据、h5ad、样本级临床字段和本地路径不进入仓库。

## 运行约定

1. 先检出并记录明确commit，再提交作业；作业日志第一行写`repo_commit`和`dirty`计数。
2. 从运行目录提交作业，并导出`XENIUM_CODE_ROOT`；`--array`必须在提交时显式给出，与样本表行数一致。
3. 下载按清单校验字节数，每个来源目录生成`SHA256SUMS`。
4. 样本表重建后，先确认已有行号没有变化，再只提交新增行。

## 迁移数据布局（2026-10-31）

旧服务器数据迁移后按旧布局放在`${PC_ROOT}`（本地配置`~/.config/pancancer/pancancer.local.sh`），运行配置中的旧绝对路径按根目录下的`PATH_MAP.tsv`做前缀替换即可对应。

| 位置 | 内容 |
|---|---|
| `Project_v5/` | 运行目录、代码、配置、冻结文件、来源记录 |
| `Project_v3/annotation_rebuild/` | 核心全基因zarr及两个注释目录 |
| `resources/` | CellTypist模型、ExperimentHub缓存、私有wheel |
| `envs/specs/` | 旧服务器环境的精确清单，用于在WEHI重建 |
| `views/` | 按分析线分组的符号链接，含每个系列的最新版本和关键输入入口 |
| `index/runs_by_line.tsv` | 运行与分析线对照 |
| `migration/` | 回执、逐文件SHA、移动记录与复核结果 |

目前位于scratch（不备份，30天未访问会被清理），持久空间待申请；旧服务器原件保留。详见[迁移记录](../runs/20261030_WEHI_old_server_migration_v1.md)。
