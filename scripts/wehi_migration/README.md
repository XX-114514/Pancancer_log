# 旧数据服务器 → WEHI 迁移（Cloudflare R2中转）

执行[按需迁移规划](../../methods/WEHI-data-migration-plan-20261030.md)的工具。旧服务器打包并上传到用户自有的R2桶，WEHI边下载边校验，解包后逐文件核对SHA256。
SSH直连只有约4–30 KB/s，所以M0也走同一通道。

## 文件

| 文件 | 运行位置 | 作用 |
|---|---|---|
| `config.sh` | 两端 | 读取本地配置`~/.config/pancancer/migration.local.sh`；清除旧服务器失效的代理变量；并行SHA256 |
| `setup_r2_remote.sh` | 两端，交互 | 写入rclone的R2 remote（密钥不回显、不出现在命令行），并做一次写/读/删探测 |
| `old_stage_m0.sh` | 旧服务器 | 用本地rsync收集M0：代码、配置、每个运行的config/provenance/scripts/tests/manifests和outputs顶层文本，外加指定文件（如私有wheel） |
| `old_pack_upload.sh` | 旧服务器 | 逐文件SHA256 → tar分块（默认1 GiB）→ 分块SHA256 → 先传清单，再在预算上限内分批并行传分块 → 写`_UPLOADED` |
| `wehi_fetch_verify.sh` | WEHI | 分块一到就下载并核对SHA256，核对通过即从桶里删除、留下`received/`标记；不一致的分块两边都删，由旧服务器重传 → 解包（不覆盖已有目标）→ 逐文件SHA256与旧服务器清单逐行一致 → `receipt.tsv` → 清空该项 |
| `old_m2_closure.py` | 旧服务器 | 从种子运行出发，递归解析config、provenance和outputs顶层的JSON中引用到的运行目录 |
| `old_stage_m2.sh` | 旧服务器 | 重新做依赖分析 → 去掉排除项 → 等引用这些运行的进程结束 → 用硬链接建M2目录（不额外占盘、不碰原件），符号链接只记录不打包 |
| `old_stage_list.sh` | 旧服务器 | 按清单（硬链接或复制）建立暂存目录，用于补充迁移和增量同步 |
| `old_export_env.sh` | 旧服务器 | 导出conda/pip/R精确清单，并复制非CRAN/Bioconductor的R包 |
| `wehi_organize.py` | WEHI | 把已核对的各项按旧服务器布局移入一个数据树，记录每次移动；默认只演练 |
| `wehi_verify_relocated.py` | WEHI | 按移动记录在新位置重新核对全部文件的SHA256 |
| `wehi_build_views.py` | WEHI | 按分析线建立运行索引和符号链接视图（含每个系列的最新版本） |
| `jobs/fetch.sbatch` | WEHI | 在datamover分区运行下载核对脚本 |
| `wehi_rebuild_env.sh`、`jobs/rebuild_env.sbatch` | WEHI | 按`envs/specs/<env>`重建旧环境：conda explicit（不求解）→ 补装`name==version`不符的PyPI精确版本（`--no-deps`）→ 复制非CRAN R库 → pip/R差异表 |
| `wehi_install_R_exact.R`、`jobs/build_R_subset_env.sbatch` | WEHI | conda-forge r-base 4.3.3 + 按旧记录版本从CRAN归档/Bioconductor 3.18装指定包及其依赖闭包（包名用冒号分隔：`sbatch --export`按逗号拆参数） |
| `wehi_oldpath_sandbox.sh`、`shims/quota` | WEHI | bwrap沙箱：按本地`PATH_MAP.tsv`把旧绝对路径只读绑定到WEHI位置、旧环境路径绑定到重建环境，只有指定输出目录可写；`quota -g <旧组>`改报当前用户在WEHI的真实配额 |
| `jobs/replay_python_stage.sbatch`、`jobs/replay_covarnet.sbatch` | WEHI | 在新目录中按原配置重放一个已接受阶段，再用`wehi_compare_outputs.py`逐文件与原输出比较 |
| `wehi_compare_outputs.py` | WEHI | TSV/TSV.GZ逐值、JSON逐键（列出不同的键）、其余按SHA256 |

注意：在一个Slurm作业内部（例如本身作为作业运行的会话）不要用`srun`，它会在当前作业里开step并等待资源；一律用`sbatch`。

两个脚本都可以重跑：已完成的哈希、打包步骤会保留，rclone跳过已传完的分块。有符号链接时打包会停止，先处理。

## 本地配置（不进仓库）

`~/.config/pancancer/migration.local.sh`，权限600：

```bash
MIG_STAGE=...          # 暂存目录：旧服务器放在空间充足的数据盘，WEHI放在scratch
MIG_BUCKET=...         # R2桶名
MIG_REMOTE=r2mig       # setup_r2_remote.sh 写入的remote名
MIG_RCLONE=...         # 旧服务器上rclone的路径；WEHI用 module load rclone
# 可选：MIG_CHUNK=1G MIG_TRANSFERS=8 MIG_PART_CONCURRENCY=4 MIG_PART_SIZE=16M MIG_JOBS=8
```

## 创建R2桶（用户在Cloudflare控制台操作）

1. R2 → Create bucket；位置提示选Asia-Pacific或Oceania。
2. R2 → Manage API tokens → Create API token：权限选**Object Read & Write**，范围**只限这个桶**，有效期设为迁移期间。
3. 记下Access Key ID、Secret Access Key和Account ID。不要发在聊天、日志或仓库里。
4. 在两台服务器上分别运行`bash setup_r2_remote.sh`，按提示输入。

## 预算控制（总预算不超过5美元）

- **存储**：`MIG_BUCKET_CAP`（默认`9G`）是整个桶的硬上限。旧服务器每批上传前都会查一次桶的已用量，超过上限就等待；WEHI核对一块就删一块。因此桶里同时存放的数据始终低于R2每月10 GB的免费存储额度。
- **操作次数**：140 GB按16 MiB分段约9千次分段上传，加上每块的创建、完成、标记、删除，以及两端每5分钟一次的列表，合计约2万次A类操作，远低于每月100万次的免费额度；下载（B类）只有几百次，也远低于免费额度。
- **流量**：R2下载不收费。
- 预计费用约为0美元。最坏情况下，假设整个140 GB在桶里放满一个月，存储费也只有约2美元。
- 不要调大`MIG_BUCKET_CAP`；如果这个Cloudflare账号的R2里还存着别的数据，要相应调小。
- R2没有消费硬上限，需要绑定付款方式。建议在Cloudflare的Billing → Notifications里设置用量提醒。
- 迁移完成后删除桶和令牌。

## 运行顺序

```bash
# 旧服务器（长任务用nohup）
bash old_stage_m0.sh "$OLD_V5" "$WHEEL"
# 一次只运行一个上传进程：上限检查按整个桶计，并行上传可能短暂超出上限
nohup bash -c 'bash old_pack_upload.sh m0_records "$MIG_STAGE/m0_tree" && bash old_pack_upload.sh immune_authority_zarr "$IMMUNE_ZARR"' > upload.out 2>&1 &

# WEHI，在运行目录提交
export MIG_CODE_ROOT=<仓库>/scripts/wehi_migration
sbatch --export=ALL,MIG_ITEM=m0_records,MIG_DEST="$DEST" "$MIG_CODE_ROOT/jobs/fetch.sbatch"
```

M2由`old_stage_m2.sh`收集为一项；排除清单只放在旧服务器本地，内容为可重新下载的输入，以及作为M1/M3单独传输的数据存储。M3的两个原生存储各为一项。
每项完成后把`receipt.tsv`的计数与SHA摘要写入`runs/`记录；目录、样本表和数据对象不进仓库。

## 已验证

在WEHI上用rclone本地后端代替R2做了端到端测试：301个文件（含带空格的路径）、91个分块，上传、下载、解包以及逐文件SHA全部一致；人为改坏一个分块后，脚本在解包前报错退出，目标目录保持为空。同时运行两端，上限设为3 MB、分块1 MB时，桶的峰值约2.1 MB；人为改坏桶内一个分块后，WEHI删掉它，旧服务器自动重传，最终逐文件一致。真实R2吞吐要等桶建好后实测。
