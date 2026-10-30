# 2026-10-30 WEHI Xenium逐切片统计处理

时间：开始约04:00 UTC（15:00 Australia/Melbourne，12:00 UTC+8）；正式运行完成约04:30 UTC。

## 用户要求

“先对所有的xenium数据继续统计处理，需要详细告诉我你是如何处理的。”

## 开始前状态

- 31个Xenium单元已入库（counts h5ad＋qc.json），面板覆盖报告已完成；无状态或效应分析。
- 仓库HEAD `b009a15`，工作区干净。下一步按[规划](../../../methods/WEHI-spatial-line-plan-20261030.md)为W2/W3。

## 计划

在不定义状态、不计算效应的前提下：每切片背景校正后的可测性、标准化与数字分群、空间自相关；跨切片测量mask与重复/分割/run一致性。规则在看结果前写进代码。

## 检查过的文件

README、STATUS、TODO、AGENTS、`methods/WEHI-spatial-line-plan-20261030.md`、`docs/WEHI_SERVER_WORKFLOW.md`、
`inventories/runtime/2026-10-30/wehi_spatial/DATA_SUMMARY.md`、`scripts/wehi_xenium/*`、冻结计划JSON、S07作者标签原始CSV。

## 执行

- 新增`section_stats.py`、`cross_section_stats.py`及两个sbatch；commit `70f2ba6`。
- 冒烟测试（1切片）→ v1全量（未提交代码，调试）→ v2正式全量（`70f2ba6`，31/31＋汇总COMPLETED，dirty 0）。
- 详细方法、结果与失败记录见[run](../../../runs/20261030_WEHI_xenium_section_stats_v1.md)。

## 关键发现（已观察）

- S18面板有69–77个TCR克隆型探针，多数处于背景；普通基因399–407个。
- 同患者两次run的伪bulk相关（中位0.987）明显高于不同患者（中位0.833），支持S18按患者作为统计单位、run作为重复。
- S16两种分割伪bulk相关0.98–0.99，背景判定几乎一致；主分割选择对伪bulk层面影响小（细胞层面未比较）。
- S07 Rep1作者unsupervised标签无法映射到Xenium cell_id。

## 未验证

- 细胞层面S16两种分割的差异；Leiden群的生物学含义（未做身份）；S07 unsupervised标签的正确映射方式。

## 下一步建议

1. 修正入库读取作者表时的科学计数barcode（`1e+05`），并为S07 unsupervised标签寻找作者映射或标为不可用。
2. W4粗身份迁移后，用T细胞内而非全切片的检出重算mask的敏感性版本。
3. S16主分割方案仍需用户决定。

## 建议commit message

Record WEHI Xenium section statistics run, measurement mask and unit concordance
