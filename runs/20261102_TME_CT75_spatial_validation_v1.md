# CT75 CM独立空间验证（20261102_TME_CT75_spatial_validation_v1）

方法：[TME-CT75-spatial-CM-validation-20261102.md](../methods/TME-CT75-spatial-CM-validation-20261102.md)；结果：[报告](../docs/reports/20261102_CT75_spatial_validation/README.md)；汇总：[inventories](../inventories/analysis/2026-11-02/CT75_spatial_validation_v1/)。
运行目录`${PC_ROOT}/Project_v5/runs/20261102_TME_CT75_spatial_validation_v1/`（逐细胞、逐spot文件和模型留在服务器）。作业时间为集群本地时间（Australia/Melbourne，UTC+11），2026-11-02 22:29至11-03 01:26，即11-02 11:29–14:26 UTC。

## 环境与代码

- `scripts/spatial_cm_validation/`。参考抽样和注释在旧路径沙箱中运行，使用重建的`scanpy`环境（CellTypist 1.7.1、scikit-learn 1.7.2，与图谱标签所用相同）。
- 新建`${PC_ROOT}/envs/conda/c2l`：cell2location 0.1.5、scvi-tools 1.4.2、torch 2.15.0+cu132、pyro-ppl 1.9.0（新版pyro去掉了`deep_getattr`，cell2location 0.1.5无法导入，因此固定版本）。配方见`jobs/build_c2l_env.sbatch`，`pip_freeze.txt`留在环境目录。
- GPU：A100-PCIE-40GB与A30，驱动580.159.03。
- 前几个作业提交时仓库为`3ed60f1`（dirty），所用脚本与随后提交的`b3820a9`逐字节一致（`git diff`为空）。之后作业分别记录`581407d`、`8049d39`、`6e170ed`。

## 阶段

| 阶段 | 作业 | 结果 |
|---|---|---|
| 参考v1（只含免疫区室） | 30646404 | 完成；因Visium需要肿瘤与基质的参考谱而**弃用**，目录改名为`reference_immune_only_superseded` |
| 参考v2（全部区室） | 30646436（bigmem） | 186个节点，参考90,410个细胞，留出14,458个，12,783个基因 |
| c2l环境 | 30646405 | 建成；之后为可导入降级pyro到1.9.0 |
| Visium全部切片合并拟合 | 30646439 | 参考回归完成（参考谱已保存）；空间模型全批量**CUDA显存不足（40 GB）而失败** |
| Visium逐切片拟合 | 30646515（array 0–18） | 19/19完成，每张切片约35–45分钟 |
| Xenium注释（23张切片、4种面板） | 30646440 | 完成；留出验证见报告 |
| 全基因可分性诊断 | 30646452 | 区室×cluster准确率0.89；75个节点中33个F1≥0.5 |
| Xenium邻域v1 | 30646495/6 | 完成；因加入谱系匹配零分布而被v2**取代**，原目录保留为`superseded_no_lineage_null_*` |
| Xenium邻域v2（主要分析与探索性） | 30646524/5 | 完成 |
| Visium分数 | 30646527 | 完成 |
| Visium参考谱控制（事后） | 交互运行，约1分钟 | 完成 |
| 汇总报告 | `report_spatial_cm_v1.py` | 写入inventories |

## 失败与修正（保留）

1. 参考v1只在免疫分析域内抽样，没有非免疫细胞；改为加入同一批单位中的非免疫细胞，免疫部分抽样规则不变。
2. 普通分区的用户内存配额被另一个会话的作业占满，参考作业改到bigmem分区。
3. Visium合并全批量显存不足，改为逐切片全批量拟合。
4. 看到探索性Xenium结果后加入谱系匹配零分布；看到Visium边相关异常高之后加入参考谱相似度控制。两项都已在方法文件中注明为事后加入。

## 验证

- 打分函数的合成测试：共变成员分数0.98，零分布中位0.48；双变量Moran能识别共享的空间模式。
- 注释：区室×cluster模型与官方模型在留出图谱细胞上重放；官方模型SHA与图谱的`reference_model_SHA256`一致。
- 未做：每个模型的重复种子稳定性；Visium的CellCharter生态位（原文的正交验证，本轮未做）。
