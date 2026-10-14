# 免疫专用单细胞档案论文图组 v1

**状态：已完成45组可审阅论文图稿（6主图、39补充图）。** 输入为已验收immune_validation_authority_v1.zarr：2,670,057细胞、9来源、400文库、270描述性节点。计数、UMAP、身份和模型保持不变。

## 图集组织

主图依次为队列与UMAP、组成与主要谱系marker、T细胞细分、T连续状态、CNV证据、whole-TME CM。补充图覆盖18主要谱系、15CM范围、队列、PCA/Harmony、深度/校正来源、捕获、技术与未确定身份。

- [全部图注](../inventories/runtime/2026-10-14/immune_publication_figures_v1/FIGURE_LEGENDS.md)
- [设计与解释方法](../methods/immune-atlas-publication-figures-v1.md)
- [最终验收](../inventories/runtime/2026-10-14/immune_publication_figures_v1/visual_review_final_public.json)
- [持续任务日志](../logs/codex/2026/2026-10-14_immune-publication-visualization.md)
- [6主图预览](../inventories/runtime/2026-10-14/immune_publication_figures_v1/previews/)

## 实际执行与验收

v1生成43图后视觉复核发现布局问题；v2/v3被文字边界/空标签数组门槛中止，全部保留。v4实际2026-10-14 06:37:22 UTC exit0；context_v1 06:35:50 exit0。最终v5精确复用42图，仅重绘T状态标题；context_v2重绘2背景图并换行捕获长标签。最终135图文件取回SHA一致，45页图册可读。

269个命名marker的原生全量均值和阳性比例对照通过（最大均值误差0），全270节点计数和49,563基因轴匹配。输出183mm宽、高<=170mm、5–7pt文字和8pt面板号；PDF/SVG保留可编辑文字，PNG为600dpi。全45图联系表和关键图放大视觉复核完成。详细展示抽样在图注写明，全量统计不抽样。

## 产物位置与科学边界

服务器：`${PROJECT_ROOT}/Project_v5/runs/20261014_immune_publication_figures_v1/outputs/{figures_v5,extra_context_v2}`。本地交付包含HTML总览、45页PDF图册、各图PDF/SVG/PNG、全部图注、匿名汇总表和ZIP。完整版不进入Git；公开仓库只收代码、图注、匿名manifest和小预览，不含患者/细胞标识。

细胞身份仍为候选/描述性节点；CNV仍为单方法证据；CM仍为组成共变。富集采样不能解释为组织自然比例，新的临床检验未执行。正式投稿前仍需完成最终生物学注释审查、患者层验证与统计论证。
