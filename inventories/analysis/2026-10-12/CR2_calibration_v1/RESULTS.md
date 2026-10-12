# CR2公式与合成校准已完成 真实筛查保持探索性

6个合成设计的β、CR2 SE、系数特异Satterthwaite自由度、p及95%区间全部与真实CRAN clubSandwich0.7.0对照一致，最大差异约1e-11。[独立对照](independent_genuine_R_comparison_v1.json)。单次/重复/不等样本/高维controls均覆盖，无患者表达参与。

30组各5000次零假设校准，固定seed2026101201。多数场景大小接近或低于名义5%；singleton30的暴露相关异方差出现7.42%，不删除这一不利结果。10相关TF的Bonferroni union在当前簇截距模拟中为0.14%—3.94%，不能由此证明真实任意误差/任意TF数可控。[完整校准](null_calibration_summary.tsv)。

因此真实CR2 p/BH/BY只做自动探索排序，不能命名已验证FDR或受体激活。原全部59131 LR情境/所有TF保留，min df4预先固定；数值独立audit完成后，入选关系继续块重采样和来源一致性/独立来源核查。原CM旧结果保持，不因为换方法得到显著晋升机制。

Python合成worker10:39:25.552075—10:39:26.993502 UTC exit0。观察脚本把末尾newline序列写成literal，执行JSON文件读法需strip该尾部后另建正确取证，原文件不覆盖；scientific inventory和校准结果本身JSON正确。真实R对照10:42:53取证通过，R4.5.2/clubSandwich0.7.0/sandwich3.1.1；隔离库在任务目录。
