# 与并发 TME 注释 v2 的兼容版本

本轮发布期间，远端新增了TME注释v2（commit `3400d49`）。该版本有179模块，包含本轮原26组中若干尚未拆分的上皮组合。直接把两份JSON拼接会产生5个新增模块同名冲突，并保留胃主/壁、前列腺基底/club、肾脏不同节段的混合组合。

后续使用TME v2工作流时，采用本目录的[185模块兼容registry](tme_v2_registry_reviewed.json)，而不是把137模块的V7扩展版覆盖到整个TME v2资源上。两份原始文件均保留：

| 资源 | 模块数 | 作用 |
| --- | --- | --- |
| V7上皮扩展 `registry_merged.json` | 137 | 112旧V7模块基础上的25项扩展，适用于本轮V7对照审阅 |
| 原TME v2 `../tme_annotation_panel_v2/registry.json` | 179 | 并发任务的独立历史资源，含Lodi作者代码与DEG参照 |
| `tme_v2_registry_reviewed.json` | 185 | 将TME v2中19个项目整理上皮模块替换为本次25项，并更新近端小管；其余159个模块内容不变 |

41个Lodi作者代码/DEG参照模块逐字段保留；7个TME v1非上皮新增模块保留。引用ID重叠的Kumar2023和Lake2023为同篇论文不同URL/范围描述，在兼容资源中保留原条目并嵌入本轮审阅范围，不覆盖作者原始材料。基因符号别名仍沿用TME v2，未移除其已有映射。

19→25的逐项替换、输入/输出SHA在[兼容manifest](tme_v2_overlay_manifest.json)；生成代码为[compose_tme_v2.py](../../../scripts/v7_epithelial_extension_v1/compose_tme_v2.py)。源版本范围改变或发生未知ID/引用冲突时拒绝合并，需要先审阅。11项检查包括Lodi模块与非上皮模块保留、混合身份退出和两份registry的字节重放。

这只是资源兼容，不代表已切换正在运行的生产任务或改变细胞标签。生产工作流的registry路径应显式记录为本兼容文件；逐基因核实状态、42条原记录未纳入以及缺失类型边界仍以本轮完整表为准。
