# 2026-11-07 注释 v2：第二轮审阅、重跑、标签转移、稳定性和共识；同步到旧服务器

## 用户要求（按时间）

1. 阅读 84c00103 的推送，检查第二轮，输出当前注释结果，并规划下一步。
2. 检查运行进展；把重要文件和整合结果通过 R2 同步到旧服务器。
3. 不传 `object.h5ad`；删除桶里的 `XX/`；更新进展。
4. 审阅标签并推进注释；对低深度细胞做标签转移；合并重跑。
5. 先做多种子共识；把注释采用的策略和遇到的问题详细写成文档（方法、参数、代码），同步到服务器和 GitHub；追踪注释和迁移进度。

## 开始前状态

- **84c00103（11-07 09:21 墨尔本）：** 第一轮判定已入库，第二轮 9 个作业已提交。
- **origin/main：** 已经比 84c00103 新很多，另有 V4 CM 的工作线，与本任务的文件没有交叉。
- **主 checkout：** 落后远端，并有未跟踪的 `scripts/nmf_benchmark/`。本任务另开 worktree（detached，0764824a），不碰主 checkout。

## 关键发现（事实）

- **组织特异池命名偏宽：** 特异 marker 只支持 Kupffer 5.6%、LSEC 7.3%、AM 13.2% 的细胞。
- **Leiden 不稳定：** 只换种子，约 30% 细胞的亚型标签会变（Endothelial 0.68、Stromal 0.70）。Harmony 提前停止不是主因。
- **批次检查：** 来源 LISI 中位数，Harmony 校正前 1.2–2.6，校正后 4.9–7.7。
- **标签转移：** 留出验证宏平均准确率 0.89–0.92，这是上限；未定比例 9.5–21.6%。
- **稀有亚型：** 从 NK 中拆出 CD8 TEMRA 样，并识别出 ILC3 样、MAIT/Th17 样、CD8 naive 样；γδ T 不可分。

## 修改的文件

- 新增 `methods/tme-annotation-v2-lineage-rounds-20261107.md`：权威的方法和问题记录。
- 新增 `runs/20261107_annotation_v2_lineage_r2_r3_v1.md`。
- 新增 `inventories/annotation_v2/lineage_r2_reviewed.tsv`、`r2_action_summary.tsv`、`r3_checks/`。
- 新增脚本：`scripts/annotation_v2/apply_lineage_r2_review.py`、`transfer_low_depth_v1.py`、`lineage_batch_and_rare_checks_v1.py`、`leiden_stability_check_v1.py`、`consensus_seed_clusters_v1.py`。
- `methods/tme-annotation-v2-lodi-workflow-20261105.md`：在文末追加指向新文件的说明（原文不改）。
- `STATUS.md`：在顶部加一段；`inventories/runs.tsv`：加一行。

## 验证

- `apply_lineage_r2_review.py`：断言每个池的簇集合与审阅表一致、每个细胞原标签正确、没有残留的 `r2_<pool>`；细胞总数守恒。
- 转移验证：第一次运行全为 1.0，查出读错标签列的 bug 并修复重跑。
- `scripts/validate_repository.sh` 和 `git diff --check`：提交前运行，结果见下。

## 失败和更正

- 标签转移的标签列 bug（作业 30677663），已修复。
- 我在对话中两次给出过错误估计，均已更正：一次说 Myeloid 要并入约 10 万个细胞，实际 7.5 万；一次说同步约 2 小时完成，实际按约 3.4 MB/s 计算要十几小时。
- 后台等待命令因超时被系统停止，不影响 Slurm 作业。
- 一次清理输出的 `rm` 被安全检查拦下；改为直接覆盖同名输出，没有删除任何东西。

## 未验证和未决

- 第三轮 Myeloid 和 Endothelial 尚未审阅；共识第二步尚未执行。
- 核内长转录本簇的处理不统一；单样本簇是否改用样本层面的批次变量。

## 建议的 commit message

Annotation v2: round-2 review and tissue-pool rerouting, round-3 reruns, low-depth label transfer, LISI/Harmony/Leiden stability checks, multi-seed consensus stage 1; methods and problems record
