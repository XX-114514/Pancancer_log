# V7 组织上皮文献扩展构建与反例验证

状态：资源构建及9项反例/重放检查通过，未进行全V7细胞注释。2026-11-05 07:04:08 UTC / 15:04:08 UTC+8取证；任务开始精确时刻未登记。本次为本地资源构建，无服务器科学作业、无新细胞数。

输入为已发布的112模块registry及100行项目marker TSV；其中26行组织上皮、121条原组合×基因记录。源TSV SHA256为`9311ff5fc7d3ae254759b2a2a9b91915513084fbf83cf948519a2b0260fde88d`；旧registry SHA256为`83f7f01a5072a325eb65dffaf65565e000df43677a8085f4fb6172079f203268`。完整输入/产物SHA见[manifest](../resources/annotation/v7_epithelial_extension_v1/manifest.json)。

输出32个组织上皮/相关分泌候选模块（新增25、复用7），合并137模块、473个不同正向基因。当前32模块119条基因记录中85条本轮核实、34条继承旧trace；42条原基因记录未纳入相应模块，仍完整记录。新副本仅修订旧近端小管模块，原始112模块文件及所有细胞标签保持不变。

命令：`python scripts/v7_epithelial_extension_v1/build.py --out resources/annotation/v7_epithelial_extension_v1`，然后运行`python scripts/v7_epithelial_extension_v1/test_extension.py`。重建测试使用临时新目录逐字节对比全部生成文件与manifest。输出目录存在时构建器拒绝覆盖；首次未发布草稿仅为纠正论文图号重新生成，未重写历史资源。

9项检查涵盖：PAX8不足以补齐TG/TPO；集合管不补齐足细胞；四种胰岛激素不作同细胞AND；缺器官背景仍保留候选；SMOC2不补齐ISC；全零测量不等于阳性；混合原组合拆分；源基因全数记账和旧资源保存；字节一致重建。首次测试有2条文件句柄ResourceWarning，已用with关闭，修正后9项通过且无该警告。肾脏图号复核纠正为Lake2023 Fig.3d p.589（健康/改变状态Visium表达），不再误写Fig.1。

本轮未执行实际RNA矩阵重聚类、DEG、CellTypist重跑、来源留出或最终命名，反例通过不代表生物学分类准确率。甲状腺参考的人体部分是组织/类器官比较，不将同篇小鼠单细胞图当人体单细胞验证。完整论文补充表未全部逐行读取。

下一步：在既有V7谱系审阅流程中显式选用新registry，按来源检查实际测量、同细胞共表达、DEG与参考映射；补尚待核实的基因及原26项未覆盖的上皮类别。详见[方法](../methods/V7-epithelial-extension-20261105.md)和[任务日志](../logs/codex/2026/2026-11-05_v7-epithelial-extension.md)。

执行环境：构建Python 3.9（仅标准库）；验证Python 3.12.14、pandas 2.2.3、numpy 2.3.5。合并registry SHA256：`b86d81388094b163029fdd4473441542f9fb164dbbfbde246aa0e720d26b9dc3`。
