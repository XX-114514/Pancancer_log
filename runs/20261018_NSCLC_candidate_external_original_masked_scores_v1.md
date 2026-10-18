# 全原候选外部固定去LR指数

## 启动前定义

在真实单细胞panel提取/独立审核仍运行时，准备这一独立输入依赖，不读取其未接受结果。输入为已接受common-native counts/TMM及原外部全因子assay d745f978，原去LR全宇宙52dfdcf1和同一私有原候选轴b87fa57f。方法见[合同](../methods/NSCLC-candidate-external-original-masked-scores-v1.md)。

源码measure_NSCLC_candidate_external_masked_scores_v1.py语法已解析；原定义的排除集、A/RMS、source gates及同slot共享关系冻结。所有有限完整指数先复现旧assay，再用直接保留、减去排除贡献及增广最小二乘核对每原候选/来源的去重叠指数。16MiB预算、500GB+1GiB预留，既有scanpy环境及PANCANCER_PREFLIGHT_ROOT执行；不新做TMM、拟合或关联。

真实部署/终态尚待，不能称独立验证或科学效应。新指数只是相容冻结池上的外部固定基因assay，不是原肺癌完整作者细胞池的精确latent复制。真实RNA接受、来源模型/原患者分组/全控制及跨来源解释仍待。科学候选轴和指标仅server/授权对话。

## 实际完成与双公式内核核对

源码6cab75c5434fa33bd495e62e254adf7a1e7b0073453b589acb3e84e5a7b954a7于2026-10-18 09:23:49.659857 UTC /17:23:49.659857 UTC+8启动，09:24:07.445174 UTC /17:24:07.445174 UTC+8完成，终态c04e2169d130d73e5ad22917cf1205c9e6f040b9f1bbeaed07e20f44995d6669。原完整assay最大误差6.76e-14，直接保留与减贡献7.47e-14，与增广最小二乘8.18e-14；同slot共享响应、全部测量可用性和两个原source gates分别保存。没有新的独立raw/TMM或独立生物学验证声明。

原科学定义、排除集、权重和源归一化保持。下一新同文库case输入按d4141c94冻结合同组装，不直接继承程序—程序问题的发送者全programme gate；没有新外部效应/P。
