# GSE299886原坐标恢复与计数行保真

## 来源与问题

[官方说明](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE299886)为MERFISH多平台TMA研究，涉及肺腺癌和间皮瘤。3份assay不是3名独立患者；core/FOV/切片/患者需分别登记。该批现有基因轴完整覆盖原固定CM02两臂及pan-T，但不因此通过表达、状态或几何验收。

原转换H5AD的obsm/spatial、center_x/y与FOV均缺失，额外材料以前的基因机会/技术校准不能替代坐标验收。完整RNA和位置在旧转换过程未正确挂接；本阶段不覆盖或删除旧矩阵、标签、归档和历史结果。

## 恢复规则

1. 由已绑定SHA的原归档member map选对应RNA/count与cell_metadata文件，不按表型效果选文件。
2. 原count CSV的cell与metadata EntityID必须用精确字符串一对一且全覆盖；不经二进制float、不模糊匹配、不填空ID。
3. 每个原RNA行ID必须等于现有barcode同一行；每个保留feature的全部原值须逐项复现现有X与counts两份存储，不仅比较总量或抽样。
4. 全部匹配后，以同一原ID取回center_x/y，保存全行原ID/metadata行号、坐标、FOV已知标记、体积及RNA向量的小型旁表。读回核原ID/坐标/RNA一致，输入SHA与归档stat前后不变。
5. 缺失FOV、未知core/患者和绝对单位继续未知；旧um_assumed不当来源证明。坐标恢复不代表正式细分状态、临床或空间互作验证。

原矩阵保留500个feature，原CSV另含50个阴性对照探针。元数据transcript_count与保留基因RNA合计的差异需另核完整对照计数，不能直接把差异当错配或擅自将其作为相同RNA量。本阶段全量复算范围为全部保留feature，不冒称已复算被去掉的对照探针。

所有ID、坐标和数组仅服务器私有。公共Git只代码、合同、技术数值范围/轨迹与SHA。实际事件和限制见[run](../runs/20261018_GSE299886_original_coordinates_v1.md)。
