# Run 小分析块CR2和LR联合TF检验校准

2026-10-12T18:39:09.600994+08:00 /UTC 2026-10-12T10:39:09.600994+00:00预声明尚未运行。[合同](../inventories/analysis/2026-10-12/CR2_calibration_plan_v1.json)。全候选设计已独立验收，真实RNA效果未计算。

先使用无患者/表达数据的6个合成设计（单次/重复/不等取样及高维controls），分别5000次iid、簇截距、依赖暴露的簇方差和t5误差。联合TF模拟保留响应间相关；主LR全TF联合p采用Bonferroni union，主BH家族保持全部59131LR情境单位，设计及推断失败内部p1/公开NA。不是基于真实p挑TF。

Python CR2采用单位working covariance和系数特异Satterthwaite自由度；min df4在真实效果之前固定。全部公式结果与CRAN clubSandwich0.7.0真实参考包对照，隔离安装在任务目录，不改变服务器或现有R环境。报告所有size及95%二项区间；小范围模拟不能证明真实RNA满足误差/身份独立条件。CR2不能消除来源/分选/状态混杂，不能把方法改变直接当生物验证。

[维护者说明](https://jepusto.github.io/clubSandwich/articles/panel-data-CRVE.html)支持有限簇时使用修正标准误/自由度，仍是依赖working model的近似；此处具体门槛、联合检验及真实分析合同为本项目设计。
