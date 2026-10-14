# 免疫专用单细胞档案论文图组

## 用户要求与初始状态

充分可视化、按顶刊发表逻辑绘制UMAP等图。2026-10-14T05:56:26.488061+00:00 / 2026-10-14T13:56:26.488061+08:00开始；正式v1对象2670057×49563、9来源/400文库、270描述群、单方法CNV及15范围CM已实际验收。已有图仅100k来源平衡QA，不作为全量论文图。并发未提交内容保留。

## 计划与规范

使用server scientific-visualization skill及既有单细胞/CoVarNet合同，核对Nature官方图稿规范。6主图：范围/UMAP、多基因证据、T描述分群、T状态、CNV、CM；补充全部18谱系/15CM、批次捕获及可用质量证据。全局UMAP全部细胞；细图seed20261014固定展示抽样，全量比例和marker统计使用所有细胞。旧UMAP坐标保持，T图是atlas embedding局部展示，不伪称新拟合。

183mm双栏、<=170mm高、5—7pt正文、8pt面板号、RGB、PDF42/可编辑SVG文字、600dpi PNG。密集点云可栅格化，轴/标签保留矢量。仅描述可用数据，不补零缺测、不画虚构显著性，不重训练/改变注释与CM。

## 已检查

正式archive obs/var/uns接口、270×49563全量CP10K/阳性/可测计数、节点与基因表；CM W/H/global/each网络表；已有plot_immune_atlas.py；正式结果/交接/源资格、AGENTS及运行/方法/daily。必要图的基因计数使用native chosen表达充分统计，不重载大RNA矩阵。

## 边界和待执行

18major含未决与边界、fine候选未晋升；单方法CNV与未评估分别画；CM是source-adjusted组成共变，独立临床/旧9CM复现未冻结。成品缺标准mito/nGene全量字段，不能用0填补；可画真实RNA分母和纠正/排除来源。所有图逐页视检、PNG尺寸/矢量文字/数据映射验收后才交付。实现与真实绘制待，失败历史另版保留。

建议commit：Build complete publication figure suite for audited immune atlas。
