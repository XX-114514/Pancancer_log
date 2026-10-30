# 旧数据服务器到WEHI的按需数据迁移规划（2026-10-30）

日期：2026-10-30，Australia/Melbourne。本文件是规划，未传输任何数据对象。它承接[下一阶段研究方案](CM-project-next-stage-20261030.md)、[WEHI空间线规划](WEHI-spatial-line-plan-20261030.md)第五节“迁移通道”，以及[服务器工作流](../docs/WEHI_SERVER_WORKFLOW.md)。
原则：**按下一阶段工作包的实际依赖迁移，不做整盘复制**；公开原始数据在WEHI重新下载；旧服务器保持权威副本，直到WEHI副本逐文件校验通过。

路径记号：`${OLD_V5}`为旧服务器`Project_v5`，`${OLD_V3}`为旧服务器`Project_v3`；`${WEHI_STAGE}`为WEHI暂存（scratch），`${WEHI_PERSIST}`为WEHI持久项目空间（待申请）。实际路径只写在服务器本地配置。

## 一 已观察的约束

### 网络（2026-10-30实测）

| 路径 | 实测 | 结论 |
|---|---|---|
| WEHI经SSH从旧服务器拉取 | 单流约4 KB/s；6流合计约30 KB/s；每次连接约30 s；RTT约354 ms且丢包 | 只适合代码、配置和小型登记 |
| 旧服务器HTTPS上传到境外（Cloudflare） | 单流约68 KB/s；4流合计约440 KB/s，随并行数上升 | 大对象需经对象存储中转，并行上传 |
| 旧服务器HTTPS下载 | 约1.5–2 MB/s | — |
| WEHI公网下载 | 约150 MB/s | 公开数据在WEHI重新下载远快于迁移 |

旧服务器本地代理当时无响应，Google不可直连，因此已有的Google Drive rclone配置不能作为当前通道。16流以上的合计速度尚未实测。

### WEHI存储

- home配额约30 GB（已用约7.6 GB），放不下数据对象。
- scratch配额充足，但不备份，30天未访问会被清理；只作暂存。
- 当前账号没有可写的WEHI项目目录。冻结对象需要持久项目空间，**需用户向实验室或管理员申请**。

### 旧服务器数据量

`${OLD_V5}`约350 GB，其中`runs/`约317 GB、383个运行目录。当前TME主线的输入链条（`runs/20261026_joined_native_atlas_delivery_v1`及其上游配置）最终只依赖三个原生存储：

| 存储 | 来源数 | 行数 | 大小 | 组成 |
|---|---|---|---|---|
| `${OLD_V3}`核心全基因zarr（2026-09-23） | 54 | 4,676,787 | 约60 GB，约8.3万文件 | layers 29.8、X 20.7、obsp 8.2、obsm 2.3 GB |
| `${OLD_V5}`公开预处理合并zarr（2026-10-12 v2） | 9 | 3,180,973 | 约41 GB | X 23.0、layers 18.3 GB |
| `${OLD_V5}`免疫治疗权威zarr（2026-10-14 v1） | 9 | 2,670,057 | 约28 GB，24,431文件 | X 14.9、layers 10.3、obsm 1.8、obsp 1.2 GB |

2026-10-23以后的运行中，小于5 GB的共332个目录、合计约32 GB；当前链条里的工作层多为MB级（例如工作注释v5交付约6 MB、身份目录约110 MB、T/NK谱系Harmony输出约0.4 GB、scVI图约0.6 GB）。

## 二 按工作包确定需求

| 下一阶段需要 | 依赖的旧服务器对象 | 迁移层级 |
|---|---|---|
| 三个接口：source_role_registry、identity_state_crosswalk、patient_view_manifest | 私有来源登记、作者样本桥接、当前工作注释及crosswalk | M0、M2 |
| WEHI W4标签迁移：用冻结T子集替换临时参考 | 免疫权威zarr | M1 |
| CM02支线、IMMUNE_ATLAS交付门槛复核 | 免疫权威zarr及其公开manifest | M1 |
| 仓库检查唯一缺口`ATLAS_SYMPHONY_WHEEL` | 私有wheel | M0 |
| 来源×注释四格诊断、1–2个癌种状态测量试点、同癌种复现 | 当前TME工作层＋三个原生存储 | M2、M3 |
| 空间W1–W3、Xenium/CosMx/MERFISH等 | 无（已在WEHI重新下载） | 不迁移 |

只要主线全来源计算还在旧服务器上运行，M3就不是阻塞项；M0–M2让WEHI先承担空间和图谱参考工作。

## 三 迁移层级

### M0 代码、配置与私有登记（小，SSH直接传）

- `${OLD_V5}`下的`scripts/`、`config/`、`frozen/`（约5 MB）、`sources/`（约29 MB）。
- 每个运行目录**只取**`config/`、`provenance/`、`scripts/`和`outputs/`里的JSON清单，不取数组、checkpoint和图件，用于在WEHI上追溯来源和复核SHA。
- 私有来源登记、作者样本/患者桥接表、`ATLAS_SYMPHONY_WHEEL`。
- 方法：`rsync -a --partial --compress`，设`--max-size`上限并显式指定include列表；传前先在旧服务器统计实际大小。预计小于300 MB；按约30 KB/s计约3小时内完成，可后台断点续传。
- 存放：`${WEHI_STAGE}/old_server/M0`，私有登记权限设为700，不进入仓库。

### M1 冻结免疫权威zarr（约28 GB，优先级最高的大对象）

- 整个store逐字节迁移，不重新打包内容；使用现有公开manifest的逐文件SHA验收。
- 先按目录打成约1 GiB的tar分块（不再压缩，chunk已压缩），每块记录SHA256；24,431个小文件不逐个上传。
- 预计：按0.5–1.5 MB/s计约6–17小时。

### M2 当前TME工作层（预计5–15 GB，按依赖闭包选取）

- 起点是当前指针及以下终点：工作注释v5交付、joined原生图谱交付、T/NK父群模型范围桥、已接受的谱系HVG/Harmony与GPU scVI输出（含保存的模型和latent）、同图Leiden/UMAP，以及10-30的CM02身份解析和ABC试点终态。
- 用脚本从这些终点的`config/spec*.json`和`provenance/`递归解析引用路径，形成闭包清单；**不按日期整批复制**。审计目录只取其验收JSON，preview与checkpoint不取。
- 旧服务器上仍在运行的链条（官方Wilcoxon、T/NK图链、ABC指标等）等其终态被接受后再取，避免复制中间状态。

### M3 原生底层存储（约101 GB，主线计算迁移时才需要）

- 核心全基因zarr（约60 GB）和公开预处理合并zarr（约41 GB）。
- 迁移前先核实哪个层是原生计数（部分来源`count_origin`为`X`），不能假设X可以重算。
- 核心store的`obsp`（约8 GB，旧图）若在M2闭包中没有消费者，可以列为可选不迁移，并在清单中明确标注。不缺省删除任何层。
- 预计：按0.5–1.5 MB/s计约1–2.5天。
- 全来源重拟合需要约400 GB可用内存，WEHI上对应bigmem分区；GPU scVI可用gpuq的A100/A30。

### 不迁移

| 对象 | 约大小 | 处理 |
|---|---|---|
| `runs/20261010_clinical_scrna_expansion_v1` | 170 GB | 公开临床队列；需要时在WEHI按accession重新下载，与旧服务器SHA清单核对。用户私下提供的表格（如TableS1）走M0 |
| `runs/20261011_GSE243013_preflight_v1` | 25 GB | 历史预检；需要时重新下载 |
| `runs/20261029_CM02_original9_UMAP_compare_v1` | 9 GB | 历史对照图；可在WEHI重算 |
| `runs/20261026_full_archive_projection_v1` | 5.5 GB | 不迁移，除非出现在M2闭包中 |
| `data_expansion_20260909`、`spatial_expansion_*`、`xenium_cosmx_20260904`、`gpu_transfer_*`、`review_bundle_*` | 约40 GB | 空间来源已在WEHI重新下载；打包件为历史交接 |
| 2026-10-23以前的其他运行目录、审计中间文件、checkpoint、preview | — | 只保留M0中的清单 |
| conda环境 | — | 按仓库yml在WEHI重建 |

迁移总量预计约135–145 GB（M0–M3），约为`Project_v5`的40%；不做M3时约35–45 GB。

## 四 传输通道

1. **SSH直传**：用于M0及任何小于约100 MB的补传。
2. **对象存储中转**：用于M1–M3。旧服务器用rclone开8–16路并行、分块上传到用户自有的桶；WEHI下载、校验、解包后删除桶内对象。候选是Cloudflare R2（旧服务器已确认能连到Cloudflare，下载不收流量费）。需要用户创建桶，并在两端本地rclone配置中写入只限该桶的密钥；密钥不进仓库、日志或对话。
3. 若用户恢复旧服务器代理，可复测Google Drive通道再决定。
4. 不使用的方式：SSH直传大对象（28 GB按30 KB/s计约11天）。

## 五 校验与验收

- 每个分块：上传前、下载后SHA256一致；解包后逐文件SHA与旧服务器源清单一致。
- 免疫权威zarr：重放`immune_authority_release_public.json`、`immune_final_interface_audit_v1_public.json`对应的文件与计数检查（9来源、400文库、2,670,057细胞、49,563特征）。
- 原生存储：在WEHI上对原生档案的65个分区重算`native_rows_sha256`、`native_ids_sha256`和`native_units_sha256`，与旧服务器配置一致。
- M2：闭包清单中每个文件SHA一致；工作注释在WEHI上的实际读取结果与旧服务器读取结果一致。
- 只有验收通过后，WEHI副本才能作为分析输入；旧服务器副本保留，直到用户确认。
- 每层迁移写入`runs/`记录：清单SHA、分块数、字节数、失败重传与验收结果。

## 六 顺序

| 步骤 | 内容 | 前提 |
|---|---|---|
| 1 | M0：代码、清单、私有登记、wheel | 无；可立即开始 |
| 2 | 申请WEHI持久项目空间；创建中转桶并配置两端rclone | 用户 |
| 3 | 16流上传实测，确认真实吞吐后定分块大小与并行数 | 步骤2 |
| 4 | M1免疫权威zarr，验收后WEHI的W4改用冻结T子集 | 步骤3 |
| 5 | 编写依赖闭包脚本，迁移M2 | 旧服务器当前链条终态被接受 |
| 6 | M3原生存储，验收后再决定主线全来源计算是否转到WEHI | 步骤5及用户决定 |

## 七 需要用户决定或提供

- WEHI持久项目空间（路径和配额）。
- 中转桶及其仅限该桶的密钥（或恢复旧服务器代理）。
- 主线全来源计算是否迁移到WEHI；若不迁移，M3可以暂缓。
- 核心store的`obsp`旧图是否可以不迁移。

## 八 停止条件

- 任一分块或文件SHA不一致：重传该分块，不用部分数据做分析。
- 闭包解析发现引用的旧服务器路径已不存在：记录缺口，不用相近版本替代。
- 持久空间未到位前，WEHI上的副本只作暂存，不作为冻结输入。
