# CT75 CM独立空间验证

## 请求

用户要求用CT75结果在独立空间转录组上验证：Visium按CoVarNet原文方法；其他平台用细胞邻域方法；可以只验证CM中的部分细胞群；注释沿用同一方案；结果同步到GitHub。

## 起点

- CT75的28个11-01运行已在当天同步到WEHI（见[CT75复现](2026-11-02_WEHI-CT75-replay.md)）。
- WEHI已有Xenium入库数据（31张切片）和GSE273952 Visium原始文件；CosMx、MERFISH和Stereo-seq尚未入库。
- 另一个会话同时在跑恶性评分作业，并向同一仓库提交；它的未提交文件没有动。

## 过程

1. 从PMC12240829全文读取原文的空间方法：cell2location（N=5，α=20，参考每群≤1,000个细胞）、W加权CM活动、Spearman共定位分数和双变量Moran聚集分数。
2. 写方法文件，与脚本一起先提交（`b3820a9`），之后才看结果。
3. 参考抽样、Xenium注释与留出验证、cell2location、细胞邻域、打分、报告，详见[run记录](../../../runs/20261102_TME_CT75_spatial_validation_v1.md)。
4. 事后加入两项控制并在方法文件中注明：谱系匹配零分布（`8049d39`）和Visium参考谱相似度控制。

## 关键发现

- 靶向面板上，CM成员节点几乎都无法按同一方案准确注释；按事先规则，Xenium不可计分。全基因下可以分开（F1 0.61–0.83）。
- Visium按原文分数看，CM01很强；但这主要是cell2location把同一信号分给几个参考谱几乎相同的节点所致。控制之后，两个CM都没有得到支持。
- N_238（单核，C00）~N_44（调节T相关，C06）在Visium（16/19名患者，p=0.004）和Xenium探索性分析（8/8单位）中都为正。

## 修改的文件

- `methods/TME-CT75-spatial-CM-validation-20261102.md`（新建，含事后补充一节）
- `scripts/spatial_cm_validation/`（新建：脚本、作业文件、环境配方）
- `inventories/analysis/2026-11-02/CT75_spatial_validation_v1/`（只含汇总表和图，没有细胞或spot行）
- `docs/reports/20261102_CT75_spatial_validation/README.md`
- `runs/20261102_TME_CT75_spatial_validation_v1.md`、本日志

## 失败与未验证

失败与修正记录在run文件中。未验证：重复种子稳定性、CellCharter生态位、CosMx/MERFISH/Stereo-seq平台。

## 下一步

见报告第5节。

建议提交：spatial: CT75 CM validation on Visium (CoVarNet method) and Xenium neighbourhoods
