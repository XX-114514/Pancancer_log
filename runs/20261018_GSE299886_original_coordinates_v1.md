# 既有MERFISH：原坐标恢复

前序[固定面板与新空间来源资格](20261018_GSE328481_spatial_bridge_qualification_v1.md)。现有14assay/4唯一probe轴的全部98panel/28axis元数据预检10:40:27完成67155cf4；进入原CM02两臂和pan-T完整的全部3个既有输入，不挑表达/空间效果。

## 结构失败与完整缺失量化

结构v1源码942e0c8ca68cbd0894db0a19d2c71c94edfb236edce14f42c2289194d5fded01于2026-10-18 10:45:47.937623 UTC /18:45:47.937623 UTC+8启动，因实际坐标不满足有限性断言停止，未给几何通过或做空间检验。保留原源码/日志/plan。v2改为完整保存所有无效坐标数量而不掩盖失败，源码02297fa545183fa6b98bb48cdb50cbb1b430774d859a102639f8c05afc2c6e0a于10:47:56.548773 UTC /18:47:56.548773 UTC+8启动，10:47:58.899608 UTC /18:47:58.899608 UTC+8完成f8cb853f9defe10d9b6dccc68f69c7d4ab9101df4198efd577f078fa314546ec。全3输入165205/205768/153981行，旧spatial/center/FOV全缺失，不能把schema完成当输入合格。

## 原归档ID诊断

源码938fa08f8aeb96e4988092dddb302e107aabccbe39e782943f2b61f3512e8f6d于10:51:08.035701 UTC /18:51:08.035701 UTC+8启动，10:51:23.177209 UTC /18:51:23.177209 UTC+8完成；全部原count cell与metadata EntityID精确唯一全覆盖，原坐标全有限；现有barcode/cell_id逐行等于原count ID。精确及decimal格式诊断两口径不是6份独立输入；实际仅3份，恢复只用无需归一化的EXACT口径。

该事实支持修复旧元数据连接，而不是认为原组织没有空间位置或用缺失坐标作阴性。此阶段不修改旧H5AD，不输出实际ID/患者行。

## 实际旁表恢复与全部RNA核查

[方法](../methods/GSE299886-original-coordinate-sidecars-v1.md)。源码ac9981d77777d5e732d541b5242630d5a4c4500d96f40df38ff41d45b4564456于10:54:34.125170 UTC /18:54:34.125170 UTC+8启动，10:55:11.785058 UTC /18:55:11.785058 UTC+8完成a1731ea5731ccc62fe8154cfb6b2d0b87d02cde12d235ed0461bad54a247e528。

全524954原RNA行/500保留feature，共262477000原值分别与现有X/counts完整逐项相等；全部原barcode行序、metadata唯一双射、全坐标有限及旁表读回一致，原矩阵SHA与归档stat保持。3小旁表约15.9MB，原数据/模型/注释无修改，未产生新空间P或细胞状态。

FOV在部分输入仍缺失，TMA core/患者关系未审；source center_x/y绝对物理单位未独立核实，不提升旧um_assumed。原550列中50阴性对照未保留于现有500列，metadata transcript_count与500列合计不同，需下一对照QC核对后解释，不能自动当RNA行错配或忽略。

下一对固定面板RNA捕获、阴性探针、原几何分区和适当技术单位完成资格；再冻结有限空间问题和null。原CM02临床代理与精确CM/互作的边界保持，完整goal仍在推进。全部科学/ID/数组只server/授权对话，公共仅方法/技术范围和SHA。
