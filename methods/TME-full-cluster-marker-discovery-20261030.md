# 当前TME同图数字群的全量marker证据

接续已完成并独立数值验收的新Harmony近邻图、主Leiden数字群和全量UMAP。原输入计数、作者标签、既有工作身份和日期CM不改；本阶段产出marker证据与旧标签对应，不自动宣布最终身份。

## 输入和标准检验

固定当前全部拟合细胞与主分辨率群成员；使用各来源均唯一实测的全部共同基因，范围不限于整合HVG，不抽样细胞。非共同基因保留在来源内全原生表达统计中，不给未测来源补生物学零。marker成员通过来源、原生行及群成员SHA冻结。

表达为历史已采用的分析计数：存在采用差值时恢复counts减ambient_removed_counts，其余继承发布计数。完整RNA总量CP10K、log1p后用于标准[Scanpy rank_genes_groups](https://scanpy.readthedocs.io/en/latest/api/generated/scanpy.tl.rank_genes_groups.html)，method=wilcoxon、tie_correct=True、reference=rest、use_raw=False，按群对全部共同基因BH校正。Harmony坐标不当作表达矩阵，不对浮点去污染估计取整。

这是用于身份发现的细胞级marker筛选，不是患者级DEG或独立临床效应。来源和癌种不均衡可能影响排序，因此同时输出来源内群计数、全部原生唯一实测基因的CP10K/logCP10K均值所需充分统计量及检测数；未观察到群和基因未测不同于真实零。后续按来源内效应、阳性/排除标志和多基因一致性审阅，不单凭极小P值命名。

## 计算与证据

核实已安装Scanpy的实际分块排名实现、真实标准Wilcoxon已知反向marker算例、当前RAM及物理/组配额。在RAM中保留共同基因稀疏分析表达，不写出重复全量表达矩阵；原生全基因统计按来源写侧表，保留原生基因轴、测量范围和确切成员。复用既有有界CSR读数，对完整原始RNA及采用后RNA逐块复算并匹配已验输入。

实际阶段先输出来源内完整表达统计、旧工作标签到新数字群的计数crosswalk，再做全量标准Wilcoxon。保存全基因检验结果及每群候选，不把第一名自动当最终显示marker。独立读回与表达/来源复核通过后再审大类、谱系内细分及一个展示marker加数字的命名。拟合来源全部参与，本轮为传导式发现，不冒称来源留出验证或完成A/B/C全部比较。
