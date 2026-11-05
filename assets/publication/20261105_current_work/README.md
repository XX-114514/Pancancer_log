# 当前研究总览｜图件与证据

配套[当前研究总览](../../../docs/homepage/07_current_work.md)。页面图号按科学论证排列；文件名保留绘制批次名称，因此页面图7对应`fig6_spatial_membership`。新增六幅主文研究图及一幅补图均提供PNG和PDF；本目录不含逐细胞表、坐标、患者信息或表达矩阵。

| 页面图号 | 科学问题与范围 | PNG / PDF | 直接来源 |
|---|---|---|---|
| 图1 | integration_v2的癌种/来源覆盖及当前TME分析范围 | [PNG](fig1_data_scope.png) / [PDF](fig1_data_scope.pdf) | [癌种](../../../inventories/integration_v2/composition_by_cancer.tsv)、[来源批次](../../../inventories/integration_v2/composition_by_batch.tsv)、[版本表](../../../inventories/server_analysis/2026-11-05/version_crosswalk.tsv) |
| 图2 | 相同20,000细胞上的六种既有嵌入：来源与marker粗身份 | [PNG](fig2_batch_umap.png) / [PDF](fig2_batch_umap.pdf) | 服务器既有TME v1 `batch_benchmark_v1`；[输入hash与参数](umap_receipt.json) |
| 图3 | 混合与结构保留的量化权衡 | [PNG](fig3_batch_metrics.png) / [PDF](fig3_batch_metrics.pdf) | [graph / scIB原表](../../../inventories/server_analysis/2026-11-05/tables/batch_benchmark/) |
| 图4 | GSE138709候选注释的路由、冲突与缺少参考 | [PNG](fig4_annotation_routes.png) / [PDF](fig4_annotation_routes.pdf) | [v2试点收据](../../../inventories/server_analysis/2026-11-05/annotation_pilots.json) |
| 图5 | V7 strict rank与NMF尺度不变性 | [PNG](fig5_cm_identifiability.png) / [PDF](fig5_cm_identifiability.pdf) | [rank](../../../inventories/server_analysis/2026-11-05/tables/v7_cm/primary_strict/rank_survey_v1.tsv)、[尺度审计](../../../inventories/server_analysis/2026-11-05/nmf_scale_sensitivity.tsv) |
| 图6，复用 | CT75参考谱控制后的空间边 | [原PNG](../../../inventories/analysis/2026-11-02/CT75_spatial_validation_v1/figures/05_visium_signature_adjusted_edges.png) | [原报告](../../../docs/reports/20261102_CT75_spatial_validation/README.md) |
| 图7 | V7全部九模块、两seed、两成员定义的空间方向 | [PNG](fig6_spatial_membership.png) / [PDF](fig6_spatial_membership.pdf) | [四份50µm表](../../../inventories/server_analysis/2026-11-05/tables/v7_spatial/) |
| 图8，复用 | 2026-10-14免疫CM02 v1固定RNA代理的历史临床关联 | [原PNG](../20261026/cm02_clinical_associations.png) | [固定代理、模型与比较家族](../../../docs/PROJECT_PUBLICATION_EVIDENCE_20261026.md) |

## UMAP读法与复现

图2来自300,000细胞、71来源标签的旧TME基准；随机种子20261105均匀无放回抽取相同20,000行。对六份已完成嵌入分别拟合展示UMAP：15近邻、min_dist=0.5、欧氏距离、300 epochs、单线程。相同点序、来源颜色及粗标签颜色跨面板固定；图不显示细胞ID。[完整71来源图例PNG](figS1_source_legend.png) / [PDF](figS1_source_legend.pdf)。

这里新增的是二维展示计算，不是新的Harmony/scVI/Scanorama训练，也不是CM重拟合。UMAP各面板坐标任意，不能比较簇面积或跨图距离；20,000细胞的诱导邻域不能替代原30万/10万评价。粗marker标签不是独立真值。BBKNN只有图表示，本次只展示其原图指标，没有以其他嵌入替代其UMAP。图3中的Scanorama scIB缺值与图方法未计算嵌入分数分别标记，不补0。

```bash
# 聚合图：只读取仓库汇总表，依赖版本见figure_manifest.json
python scripts/render_current_work_figures.py

# 研究环境：只读取已有obs与嵌入，不读取表达矩阵；输出目录必须为新版本
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
python scripts/render_current_work_umap.py \
  --benchmark-root "${PROJECT_ROOT}/tme_v1/batch_benchmark_v1" \
  --output "${VISUAL_ROOT}/figures" --n 20000
```

[聚合图脚本](../../../scripts/render_current_work_figures.py)及[输入/输出hash与本地环境](figure_manifest.json)；[UMAP脚本](../../../scripts/render_current_work_umap.py)及[嵌入hash、抽样轴hash与服务器环境](umap_receipt.json)。PNG用于GitHub阅读，PDF保留矢量文字/图形，密集散点栅格化。图中的汇总条、点和线来自现有统计表，未新增科学检验；没有误差估计的汇总不自行添加误差条。

## 原归档补图

[来源UMAP](atlas_umap_source.png)和[16基因marker](atlas_umap_markers.png)来自`integration_v2/run_v2/figures`，原样下载并按[原图收据](atlas_receipt.json)核对hash。4,720,931细胞的数据层和正在验收的新TME run_v2不是同一对象。来源原图分辨率较低；marker原图各基因按99.5分位设置色阶、没有独立色标，原脚本将缺测绘为0，因此仅保留作历史定位与来源追踪，不作为正式测量mask或身份验证图。主文已明确这些局限。

[作图执行与质量检查](../../../runs/20261105_current_work_visualization_v1.md)。科学状态仍按各分析自己的取证时刻维护，不随图件生成刷新。

## 发布前新增：V7完整来源留出与材料分层

`fa2d26a`经`4285a4b`并入main后，补图S2增加[PNG](figS2_v7_source_material.png) / [PDF](figS2_v7_source_material.pdf)。数据直接取[53来源留出表](../../../docs/reports/20261103_V7_CM_blisa/tables/v7_primary_K9_LOSO_r.tsv)和[材料活动份额](../../../docs/reports/20261103_V7_CM_blisa/tables/CM_activity_share_by_sample_class.tsv)，不重新拟合或推断因果。CM04原报告写32/53达到r≥0.8；公开三位小数表重数为33/53（其中一值恰为0.800），保留精度/口径差异，待未舍入结果确认。图内按公开表重数并在图注说明。共享节点Pearson匹配与九算法审计的行归一余弦是不同指标，不能混同。
