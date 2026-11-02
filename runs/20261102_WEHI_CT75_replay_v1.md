# CT75可视化与协变审查：WEHI复现

目标：在WEHI上按原配置重放[2026-11-01 CT75检查点](../logs/codex/2026/2026-11-01_TME-CT75-visual-review-public-checkpoint.md)涉及的已接受阶段，逐文件与旧服务器原输出比较。2026-11-02 01:30–02:20 UTC执行（09:30–10:20 UTC+8）。本记录只含技术状态，不含来源标识、表型名称、科学效应或矩阵。

## 增量同步（M6）

- 旧服务器用`old_m2_closure.py`从5个CT75种子运行解析依赖闭包：236个运行，其中WEHI缺的只有28个2026-11-01运行（另有一个按迁移计划不迁移的预检运行和两个无效引用）。10-26至10-31的上游全部已在WEHI。
- `old_stage_list.sh`（copy）→ `old_pack_upload.sh`：1,234个文件、123 MiB、1个分块，01:37:40 UTC上传完成；WEHI作业30636926在01:42:20 UTC `VERIFIED`并从桶中清除。桶容量远低于9G上限。
- 28个运行移入`${PC_ROOT}/Project_v5/runs/`，按旧服务器清单在新位置复核1,232个文件SHA256，全部一致。
- 读取R2版本状态的403（AccessDenied）为已知的无害探测，不影响校验。

## 重放

工具：[replay_python_stage.sbatch](../scripts/wehi_migration/jobs/replay_python_stage.sbatch)、[replay_covarnet.sbatch](../scripts/wehi_migration/jobs/replay_covarnet.sbatch)（repo_commit `a8465a8`，dirty 0），CPU作业均加`--constraint=AVX512`。Python阶段用重建的`scanpy`环境；CoVarNet用`scrna_r_covarnet2`。每个阶段只读原上游输入（隔离重放）。

| 阶段（原运行） | 作业 | 结果 |
|---|---|---|
| 免疫CM面板准备 | 30637028 | 56个数据文件逐字节一致；26个JSON只差输出路径/路径相关SHA |
| 主方案CoVarNet拟合（K由秩调查选出） | 30637035 | 输入与作者代码SHA通过；秩调查、全部网络边/节点表、权重等24个文件逐字节一致；W/A最大差≈1e-16，H≈1e-14；RDS二进制不同（未判定，其导出表一致）；`numeric_readback_acceptance`由Python包装层生成，本作业不产生 |
| 免疫CM来源/背景复现 | 30637029 | 2个数据文件一致；接受JSON只差路径 |
| CT75全局边来源/背景 | 30637031 | 3个一致；接受JSON只差路径 |
| CT75固定K4来源/背景 | 30637032 | 3个一致；接受JSON只差路径 |
| CT75共享cluster条件丰度诊断 | 30637034 | 3个一致；接受JSON只差完成时间 |
| CT75 UMAP与CM图（渲染） | 30637030 | 全部PNG与表格逐字节一致；9个SVG只差matplotlib随机marker id和日期（规范化后一致）；摘要只差时间。sbatch状态FAILED仅因比较脚本对缺失文件返回非零 |
| CT75原生表达/程序/配体受体兼容性 | 30637033 | 33 min（旧服务器约18 min）；11个文件（含两个NPZ、全部LR表和4张图）逐字节一致；摘要只差时间 |

## 未能复现的部分

- UMAP运行目录下的`secondary_source_review_v1`、`capture_review_v1`、`same_cell_review_v2`，以及图包zip、带注释的LR top20表、参考概率诊断表和导出清单：生成脚本不在运行目录中，旧服务器上按文件名搜索也未找到（搜索在500 s超时内无命中，结论不完全）。推测由本机Codex工作区经SSH执行。需要用户提供这些脚本才能重放。
- `visual_addendum_v1`与共享cluster运行的输出逐字节相同，后者已重放一致。

## 结论

可检查的所有已接受阶段在WEHI上数值复现：数据表和PNG逐字节一致，NMF因子在1e-14内一致，CM成员、网络边和来源/背景统计完全相同。差异只有路径、时间戳和SVG随机id。复现输出在`${PC_ROOT}/Project_v5/runs/20261102_WEHI_CT75_replay_v1/`（不进仓库）。
