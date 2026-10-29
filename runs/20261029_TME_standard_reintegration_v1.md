# 当前TME标准HVG/PCA/Harmony重整合

用户要求以已有流程和有据方法重新整合当前TME，并核实历史doublet/环境RNA。按[标准方法](../methods/TME-standard-reintegration-20261029.md)在独立新运行目录执行；原5M/免疫图谱、原始计数、旧身份及日期CM保持。具体输入source/精确细胞轴与已采用分析层在私有冻结合同，不上传个体或矩阵。

主阶段：标准Scanpy seurat_v3按来源选择HVG，沿用旧6500候选→4000非技术基因，完整已采用RNA归一化后取HVG，ARPACK中心化PCA50，直接完整当前TME的Harmony（来源batch，最多20轮，固定种子和8线程）。原始PCA保留，不复用旧全档案参考查询坐标。整数counts用于HVG，PCA使用已采用ambient扣除的分析计数，二者不混称；原QC/去污染不重做。

执行03:38:10 UTC /11:38:10 UTC+8 /14:38:10 Melbourne真实启动；实际标准LOESS/HVG接口测试通过，完整范围和共同唯一实测特征核对后开始native integer输入。取证03:44:09 UTC进程仍活跃、有来源输出；完整标准HVG/PCA/Harmony和整合生物学审阅尚待。审核/同图近邻Leiden及UMAP链03:43 UTC启动等待，有限24小时监督，实际退出码与进程startticks保留。观察超时不构成科学失败，不据此重跑。

新群体图只显示数字Leiden；同一connectivities用于Leiden和UMAP，原生marker与正负标志审阅再命名。不得以Harmony收敛或图片分离代替身份/状态保留评价。详见[同任务日志](../logs/codex/2026/2026-10-29_TME-standard-reintegration.md)。
