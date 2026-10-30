# 2026-10-30 旧数据服务器到WEHI的按需迁移规划

时间：约06:00–07:30 UTC（17:00–18:30 Australia/Melbourne，14:00–15:30 UTC+8）。

## 用户要求

“我已经可以在这个服务器连接到我[旧数据服务器]了，现在可以根据情况研究如何迁移数据分析。”
随后：“根据pancaner项目最新的情况规划迁移数据，不要全部迁移，根据需要。”

## 开始前状态

- 仓库HEAD `64f0056`，工作区干净。
- WEHI已承担空间来源下载与入库；冻结免疫图谱、私有登记和旧运行产物只在旧服务器上（[工作流](../../../docs/WEHI_SERVER_WORKFLOW.md)）。
- 旧服务器上主线仍在运行：官方Wilcoxon、T/NK谱系Harmony/scVI图链、ABC指标等。

## 计划

测量两台服务器之间的网络通道；按下一阶段工作包追溯实际依赖；给出分层迁移清单、通道、校验与顺序。不传输数据对象。

## 检查过的文件

README、STATUS、AGENTS、`methods/CM-project-next-stage-20261030.md`、`methods/WEHI-spatial-line-plan-20261030.md`、
`docs/WEHI_SERVER_WORKFLOW.md`、`docs/IMMUNE_ATLAS_V1_HANDOFF.md`、`docs/DATA_POLICY.md`、
`inventories/analysis/2026-10-23/full_context_annotation_v1/`中10-29至10-30的技术收据。
旧服务器（只读）：T/NK谱系Harmony运行的`config/spec_v1.json`、joined原生图谱交付的配置、原生档案范围`native_archive_collection_v1.json`及其`frozen_context`。

## 执行（只读）

- SSH吞吐：不同负载大小及6路并行的`head -c | wc -c`计时；SSH上传方向计时。
- 旧服务器出站：带/不带本地代理的HTTPS下载、单路与4路并行上传到Cloudflare测速端点。
- 大小：在旧服务器后台对`Project_v5`顶层、`runs/*`、三个原生存储及其组成部分执行`du --apparent-size`，结果写在旧服务器用户缓存目录。
- WEHI：`quota`、可写项目目录检查。

## 关键发现（已观察）

- SSH单流约4 KB/s，6流合计约30 KB/s；旧服务器HTTPS上传单流约68 KB/s，4流合计约440 KB/s；WEHI公网下载约150 MB/s。旧服务器本地代理无响应。
- 当前TME主线（joined原生图谱→原生档案范围）最终依赖三个原生存储：核心全基因zarr（约60 GB）、公开预处理合并zarr（约41 GB）、免疫权威zarr（约28 GB）。
- `runs/`约317 GB，其中临床scRNA扩展170 GB、公开预处理43 GB、免疫图谱29 GB、GSE243013预检25 GB；当前链条的工作层多为MB到百MB级。
- WEHI home配额约30 GB；scratch会清理；当前账号没有可写的项目目录。

## 修改的文件

- 新增`methods/WEHI-data-migration-plan-20261030.md`：分层迁移规划（M0–M3与不迁移清单）。
- `docs/WEHI_SERVER_WORKFLOW.md`：通道现状改为“SSH可用但只适合小文件”，链接规划。
- `STATUS.md`、`logs/daily/2026/2026-10-30_WEHI-spatial-server.md`：追加节点。

## 验证

`git diff --check`；检查新增文本中没有主机名、账号、绝对路径和凭据。

## 失败尝试

- 首次用`rsync`和`cat`整文件测速，在超时内没有完成任何可计量传输；改为小负载计时。
- 8路并行上传测试超时未返回；改为4路×1 MB。
- 第一次`du * runs/*`因du去重未列出`runs/*`，改为单独扫描。

## 未验证

- 16路以上并行上传的真实合计速度。
- M0实际大小（预计小于300 MB）和M2闭包实际大小（预计5–15 GB）。
- 核心store中哪个层是原生计数，`obsp`是否有当前消费者。

## 未解决的问题

见规划第七节：持久项目空间、中转桶、主线计算是否迁移、`obsp`是否迁移。

## 建议的下一步

立即执行M0（SSH，小）；用户提供持久空间和中转桶后实测16流吞吐，再迁移M1。

## 建议的commit message

Plan needs-based data migration from the old server to WEHI

## 续：确定全部迁移，采用Cloudflare R2（约07:30–07:45 UTC）

用户：“确定全部100多GB，考虑cloudflareR2方案”。

- 旧服务器16路×1 MB并行上传到Cloudflare，合计约1.1 MB/s（14 s内16 MB，含连接建立）。两端都能解析并连接R2端点；用假账号ID时TLS被拒，属于预期，真实连通要等有账号后验证。
- WEHI有`module load rclone`（1.67）；旧服务器有rclone 1.75；旧服务器数据盘剩余约3.7 TB，可作暂存。
- 私有wheel位于旧服务器某个运行目录（2026-10-23参考运行时检查）；旧仓库检出没有未跟踪的私有登记，私有登记应在运行的outputs中，由M0收集。
- 新增`scripts/wehi_migration/`：配置、R2 remote设置（密钥不回显、不进命令行）、M0本地收集、打包上传、下载校验、datamover作业。
- 验证：所有脚本通过`bash -n`；用rclone本地后端做端到端测试（301个文件、91个分块，路径含空格），逐文件SHA一致，`diff -r`一致；改坏一个分块后，在解包前失败退出；M0的include规则在合成目录上排除了npz、子目录和checkpoint。
- 未验证：真实R2吞吐与费用；M0实际大小（旧服务器后台统计时SSH超时，未读到结果）。
- 下一步：用户创建桶和只限该桶的令牌，在两端运行`setup_r2_remote.sh`；先传M0和M1，并实测吞吐。

## 续：总预算不超过5美元（约09:20 UTC）

用户要求整体预算不超过5美元，并请我运行`setup_r2_remote.sh`。我没有运行：桶和令牌还没建好，密钥应由用户在自己的终端输入，而且本会话没有交互终端。

- 改为预算硬控制：旧服务器每批上传前查整个桶的已用量，`MIG_BUCKET_CAP`默认9G，低于R2的10 GB免费存储额度；WEHI核对一块就删一块，并写`received/`标记；核对失败的分块两边都删，由旧服务器自动重传。
- 估算：A类操作约2万次（免费额度每月100万），B类几百次，下载不收流量费；预计约0美元，最坏情况低于2美元。
- 验证：两端同时运行，上限3 MB、分块1 MB时桶峰值约2.1 MB；人为改坏桶内一个分块后自动重传，最终逐文件一致。
- 上传一次只运行一个进程，避免两个进程同时判断余量。

## 续：完整启动迁移（约09:50–10:20 UTC）

用户已在旧服务器配置`r2mig`和WEHI专用令牌`r2wehi`，并要求完整启动迁移。

- `r2wehi`配置段经SSH从旧服务器写入WEHI的rclone配置（权限600，内容未显示）；WEHI上写、读、删探测通过；随后从旧服务器配置中删除该段，两台服务器各自只保留自己的令牌。
- M0：13,265个文件、1.8 GB、2个分块；约8分钟完成上传、下载、逐文件核对（VERIFIED），桶内该项已清空。清空时出现读取版本设置的403警告，是令牌权限范围所致，不影响结果。
- M1（免疫权威zarr，28个分块）上传中，第一批8个分块实测约6.8 MiB/s，远高于curl测速。
- M2依赖分析：从10-29、10-30的161个运行出发，闭包包含291个运行（约315 GB），其中多数是上游原始下载。排除临床scRNA扩展、GSE243013预检、CM02历史UMAP对照，以及已作为M1/M3传输的两个数据存储；另加Project_v3的celltypist标签和表达注释下游结果。估计约35 GB，高于原先估计的5–15 GB。Project_v3旧CoVarNet运行（185 GB）和核心zarr的上游输入（22 GB）不迁移。
- 新增`old_m2_closure.py`和`old_stage_m2.sh`，在合成目录上测试通过（递归依赖、排除、符号链接记录、原件未改动）。
- 旧服务器上传顺序：M0 → M1 → M3核心 → M3公开 → M2（M2打包前重新做依赖分析，并等引用这些运行的进程结束）。每一项只在前一项上传完成后才开始，同时只有一个上传进程。
- WEHI：每项一个datamover作业；受用户CPU配额限制，最多同时运行两个，其余排队。
- 预算：桶硬上限9 GB；按GB-月计费，进行中的分段上传只是短暂占用，预计仍在免费额度内。

## 续：完成状态检查（2026-10-30 22:40 UTC）

5项全部VERIFIED，桶为空，WEHI上的datamover作业全部结束。详见[run](../../../runs/20261030_WEHI_old_server_migration_v1.md)。原定12:18 UTC的检查提醒没有留下仓库记录，本次补写运行记录。

## 续：按三条在研分析线核对补充迁移（2026-10-31）

用户要求根据三条在研分析（新TME的CM发现与复现性、免疫图谱CM结合患者信息、状态NMF）读取日志和服务器文件，找出仍需迁移的关键文件。

- 读取：`methods/TME-annotation-CM-state-goal-20261030.md`，10-28至10-30的CM02与状态日志，迁移后新增的提交。
- 旧服务器只读审计：快照后新运行594个（约1.5 GB）；已迁移运行在快照后没有文件被更新；快照前未迁移的运行151个，其中148个小运行约9.6 GB；对已迁移和新运行的JSON做引用解析，发现直接引用临床扩展运行中的9个来源对象（约20 GB）和临床侧表，以及R/Python环境、CellTypist模型、ExperimentHub缓存。
- 结论写入规划第九节（M4）。尚未传输。

## 续：启动M4，并整理本机文件清单（2026-10-30 23:05 UTC）

- 新增`old_stage_list.sh`（按清单用硬链接或复制建立暂存目录，可用于以后的增量同步）和`old_export_env.sh`（conda/pip/R精确清单，并复制非CRAN/Bioconductor的R包）。用合成目录测试了收集清单脚本。
- 旧服务器按顺序运行：m4a_inputs（9个来源对象、临床扩展的非原始部分、快照后594个新运行（复制）、CellTypist模型、ExperimentHub缓存；99,132个文件，24 GB）→ m4b_runs（148个运行，硬链接）→ m4c_env（scrna_r、scanpy、singler_bioc、scvi）。启动时没有进程在使用这594个新运行。
- WEHI：三个datamover作业，目标目录在scratch的`m4/`下。
- 本机文件：根据仓库日志整理出只在用户本机的资料类别（CM02私有交付与临床来源桥、用户上传的作者资料原件、待补的Table S1、黑色素瘤IMC公开分段下载、审阅图册/工作簿/研究稿）。本会话无法访问用户本机，清单需用户核对。
