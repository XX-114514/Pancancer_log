# WEHI Xenium逐切片统计、测量mask与单位结构（2026-10-30）

对已入库的31个Xenium单元做描述性统计。不定义细胞状态，不计算效应；Leiden群只是编号，不是身份。
对应[WEHI空间线规划](../methods/WEHI-spatial-line-plan-20261030.md)的W3（Xenium部分）和W2（单位结构）。任务日志见[codex日志](../logs/codex/2026/2026-10-30_WEHI-xenium-section-stats.md)。

## 代码与环境

- 代码：commit `70f2ba6`，`scripts/wehi_xenium/section_stats.py`（每切片，array）与`cross_section_stats.py`（汇总）；作业`jobs/03_section_stats.sbatch`、`jobs/04_cross_section_stats.sbatch`。
- 环境：`${PC_SPATIAL_ENV}`（scanpy 1.11.5、squidpy 1.8.2、scipy 1.17.1、pynndescent 0.7.0、leidenalg 0.12.0）。
- 正式运行：31/31切片与汇总作业COMPLETED，作业日志首行均为`repo_commit 70f2ba622d44 dirty 0`。每切片1分15秒–10分09秒。

## 方法（运行前写定）

| 步骤 | 规则 |
|---|---|
| 细胞 | 只用`qc_pass`（基因计数≥10）细胞 |
| 背景 | 从源h5重读每个阴性对照特征（Negative Control Probe + Codeword）在qc_pass细胞中的总计数；基因总计数 > 本切片最大阴性特征总计数 → `above_background` |
| 特征类别 | S18面板含`<run>_<CDR3>_TRA/TRB`克隆型探针（69–77个/切片），单列为`tcr_clonotype_probe`，不进PCA，不计入普通基因数 |
| 标准化 | 每细胞缩放到本切片中位库大小，log1p；原始计数存`layers['counts']` |
| 降维聚类 | 靶向面板用全部普通基因；C05 Prime 5K用2000个`seurat` HVG（`seurat_v3`需scikit-misc，锁定环境中没有）；scale(max 10)，30 PC，kNN k=15（pynndescent），Leiden res 0.5（igraph，seed 0），UMAP |
| marker | 每个Leiden群Wilcoxon，前10 |
| 一致性 | Leiden vs 10x graphclust、作者标签：ARI、NMI |
| 空间自相关 | squidpy 6-NN质心图，Moran's I，解析P，切片内BH（只对有定义P的基因） |
| 测量mask | 每切片×基因列表：全部基因above_background=measurable，部分=partial，无=unmeasurable；来源级取最差切片；缺失基因不删 |
| 单位结构 | S07 Rep1 vs Rep2；S16同TMA unimodal vs multimodal；S18同患者run1 vs run2：伪bulk log1p CPM的Pearson/Spearman、检出率Pearson、背景判定一致率；S18另算所有不同患者切片对 |

## 结果（来源级汇总）

- 普通基因above_background：C05 4994/5001；C11、S10–S12、S16全部或差1–2个；S07 232/313与288/313；S18 257–364/399–407。
- Leiden群9–20个；与10x graphclust ARI 0.28–0.67；S07与作者supervised标签ARI 0.54/0.49。
- 单位一致性：S07技术重复Pearson 0.992；S16两种分割0.982–0.993；S18同患者两次run中位0.987（0.971–0.992，n=7），不同患者切片对中位0.833（0.519–0.940，n=129）。
- 测量mask（来源级，最差切片）：没有任何来源对全部冻结列表measurable。interferon在C11、S07、S16、S18为unmeasurable（6基因均不在面板）；
  S18 cytotoxic为unmeasurable，只由一张切片（一位患者的run1，23,177细胞）造成，其余16张为partial（缺CTSW、CD8B）。
  在面板上但处于背景的基因只出现在S07和S18：GZMK、PDCD1最常见，S18该切片另有GZMB、PRF1、NKG7、CD8A、GZMA、CCR7。

## 失败与修正（保留）

1. 第一次冒烟测试失败：scanpy 1.11.5对pynndescent返回的`csr_array`调用`getnnz()`，scipy 1.17已移除。改为薄封装，把结果转为`csr_matrix`，算法与参数不变。
2. C05第一次失败：`seurat_v3`需要scikit-misc。改用`seurat`，不改锁定环境。
3. v1全量运行（未提交代码）发现：S18两张切片各有1个零计数TCR探针，Moran P为NaN，squidpy自带BH使整列NaN；改为只对有定义P的基因做BH。
   同时把TCR克隆型探针移出PCA并单独计数。v1输出已被v2（`70f2ba6`）覆盖，v1只作调试记录。
4. S07 Rep1作者unsupervised标签（GSM7780153 Fig 3）：barcode范围0–170,429，超过cells表的167,780个细胞；与supervised标签一致性为随机水平（ARI≈0.01，+1平移也一样），
   所以不是Xenium cell_id。该ARI不可解释，在重新映射之前这组标签不可用。另：作者supervised表把barcode 100000写成`1e+05`，每个文件有1个细胞没有标签。入库未重跑。

## 结论边界

- `above_background`是可检测性判断，同时受探针性能和该切片中表达细胞多少的影响；一张T细胞稀少的切片就会让来源级mask变成unmeasurable。它不是面板失效的证据。
- 最大阴性特征阈值较严；只用于资格判断，不用于定义状态。
- 输出（每切片`stats/`、`*.norm.h5ad`、`summary/STATS_REPORT.md`及TSV）留在服务器scratch，不进仓库。
