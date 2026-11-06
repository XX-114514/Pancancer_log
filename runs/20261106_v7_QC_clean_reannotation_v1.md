# 20261106_v7_QC_clean_reannotation_v1

状态：新QC视图/输入独立审计PASS，T6与Mast clean refit已实际启动，末端marker/身份待完成。事件与取证时间见[结果登记](../inventories/analysis/2026-11-06/v7_qc_clean_reannotation_v1/README.md)。

授权：用户2026-11-06明确同意清理。输入为冻结V7（4,676,787细胞）及20261105 Mast/T6全基因raw counts和文献复核证据；输出为`${PROJECT_ROOT}/runs/20261106_v7_QC_clean_reannotation_v1`，不覆盖旧annotation_rebuild对象。

事件UTC / UTC+8：A00:26:31/08:26:31完成；B00:30:16—00:31:35/08:30:16—08:31:35，terminal0；C00:39:54/08:39:54完成；D两谱系00:42:18/08:42:18启动，各16线程。主机/PID/完整日志保留私有。R源body SHA256为3c8dadc8819623e07e736fe137ff09f9bb5a044021a6744ec3c18b7a963338fc。

软件：scanpy环境Python/NumPy1.26.4/pandas2.3.3/anndata0.11.4/Scrublet0.2.3；scrna_r环境Seurat5.3.1/Harmony1.2.4。full feature轴95,038；清理前后counts逐元素一致。实际代码body hash见登记；公开Python只增加PROJECT_ROOT环境引导，私有政策中的逐单元清单不公开。

逻辑命令（需既有服务器输入/私有政策；脚本拒绝覆盖新run，不可盲目重跑）：

```bash
python scripts/v7_qc_clean_reannotation_20261106/profile_cleaning_evidence.py
python scripts/v7_qc_clean_reannotation_20261106/supplement_unscored_doublets.py
# 先完成直方图和政策裁决，再生成新视图
python scripts/v7_qc_clean_reannotation_20261106/apply_cleaning_views.py
python scripts/v7_qc_clean_reannotation_20261106/audit_clean_views.py
Rscript scripts/v7_qc_clean_reannotation_20261106/clean_lineage_refit_v1.R "${RUN_ROOT}/stage_c_clean_views/T6/refit_config.json"
Rscript scripts/v7_qc_clean_reannotation_20261106/clean_lineage_refit_v1.R "${RUN_ROOT}/stage_c_clean_views/Mast/refit_config.json"
```

真实运行由detached supervisor记录退出和complete；不以PID消失、文件存在或CPU进度代替终验。aggregate_clean_markers.py在REVIEW_LINEAGE指定谱系且真实terminal0后运行。R模型不使用旧簇、CellTypist/scATOMIC标签作特征。

已知纠正：C的初版旧簇聚合表holdout计数用了不一致的字符串；独立审计新表直接从掩码重算，掩码和模型输入本身正确。警告单元的1个候选reason13属于数值待审，不能沿初版泛化显示字符串误解为锚点不足。其余失败/警告和终验在后续同任务记录追加。

## 后续实际完成（v2优先）

Mast额外NK组件核对36个均校准有效：01:02:17 UTC / 09:02:17 UTC+8实际追加排除；独立raw重算和精确count子集通过，当前全局2003排除／236暂存／4674548保留。Mast-v2于01:03:15—01:09:07 UTC / 09:03:15—09:09:07 UTC+8，terminal0；T6于00:42:18—01:09:01 UTC / 08:42:18—09:09:01 UTC+8，terminal0。R4.3.3、Seurat5.3.1、Harmony1.2.4；Python3.10.19及版本见结果软件收据。

两个模型各95038基因完整RNA轴；source-res.4为12 T6／18 Mast诊断簇，不等于细分数。T6共有11091基因DE，Mast12929。T6扩展118基因／1416条、Mast110基因／1980条marker均独立重算，CP10K误差<1.1e-12。当前无正向DE分区T6=0/7/10/11，Mast=1/5/7；T6初始化Quick-TRANSfer10条、tiny2 DE跳过。工作大类按多marker簇及来源人工复查，单细胞阳性规则作为独立证据，未要求每个细胞足量TCR。侧表01:19:39 UTC / 09:19:39 UTC+8写出并读回一致，fine final接受0。

详细结果、原样TSV、工作标签、掩码合同、警告、转移表与10张图见[当前登记](../inventories/analysis/2026-11-06/v7_qc_clean_reannotation_v1/clean_refit_review_v2/README.md)。完整对象留私有，旧V7 SHA在最终取证再次匹配。最新Mast模型是stage_i_Mast_NK_clean_refit_v2，T6仍为stage_d_T6_clean_refit；v1各对象保留。
