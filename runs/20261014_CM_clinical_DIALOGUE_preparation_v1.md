# CM临床患者机会核验与DIALOGUE独立实现

2026-10-14继续active目标；上一goal turn为progress，336真实模型及受体检查改变下一动作。完整目标仍为具有生物或临床价值的互作，不以环境安装/输入准备替代。

## 临床入口

只读现有同SHA210行clinical sidecar及3874旧细胞群计数，所有3来源QC总量与当前2670057对象冻结一致。按原治疗时点/组织/臂、来源患者及双方捕获计算16机会，无新p。[真实机会与限制](../inventories/analysis/2026-10-14/clinical_CM_patient_opportunity_v1/RESULTS.md)。精确部位明细移到私有工作目录保存，public只粗汇总。代码路径通过`--sidecar --counts --output`提供，Python3.12/pandas；失败首次父目录缺失未写产物，修正后成功；独立16逐库复算通过。

## DIALOGUE方法准备

官方固定commit `9c146ccf28d7706aaa60d00947a9126b4e75fd69`，R4.3.3 scrna_r复用，仅新私有库。官方CRAN版本：PMA1.2-4、unikn1.0.0、ppcor1.1、UpSetR1.4.1、lmerTest3.1-3；最新版lmerTest3.2-1会引入额外reformulas，保留下载但不使用。源归档及SHA留私有run/public provenance。没有新表达矩阵下载。

v1 UTC00:38:32—00:38:58失败：PMA配置的conda编译器不在PATH。v2 UTC00:41:22—00:41:49：PATH修正后所有官方包安装成功，合成输入缺少PC列名触发索引维度错误。v3 UTC00:42:54—00:43:02：列名修正，PMD进入HLM，但所有混合模型值NA，原实现捕获错误后最终空表失败。纯合成显式模型取证确认`chm_factor_ldetL2`未由Matrix提供：共享lme4 1.1-34与Matrix1.6-5的二进制接口不匹配。不是本研究数据的阴性结果。

v4在新`library_v4`针对当前Matrix重编译**同版本**lme4，保留原共享R库和所有失败。合成用24模拟sample、两类各720模拟cell、2共同潜变量/100模拟gene；真实cell0。实际终态/验收之后补录，尚未完成真实DIALOGUE拟合或临床机制。

## 可比单位与后续

已有独立状态任务正在按5来源240533原细胞准备完整RNA/原数字身份/患者块；读取其已提交合同而不另建重复原生输入。DIALOGUE来源内优先CRC未分选T—巨噬，患者整块发现/留出需先冻结；后续分清连续状态、身份混合、共同炎症和空间。现有T—巨噬不能替代其他谱系/全9CM；原空间阴性及所有模型敏感记录保持。

当前任务日志：[持续记录](../logs/codex/2026/2026-10-13_finite-CM-communication-and-pivot.md)。
