# 冻结V7的doublet过滤溯源方法

本轮核验的对象为54来源、4676787细胞的冻结V7及其同轴全基因core。区分：原工具执行、原工具预测、实际剔除、后续继承和本轮新增复查。原标签/矩阵不改，不重新拟合Scrublet，不用95分位重写历史判定。

## 两个历史处理分支

旧队列：`scripts/run_phase06_sample_qc_worker.py`的`Scrublet.scrub_doublets()`按原样本调用，expected rate 0.04；在应用基础QC排除前计算分数/预测，随后把doublet原因与低genes/counts/高mt一起合并到filter_reason。`filter_reason == pass`形成真正写出的filtered对象。每样本最多移除8%；<500细胞不运行模型，分数NaN、False是占位值。33队列记录中当前V7只采用31队列，必须按实际来源范围取数。

新增队列：`01_new_cohort_prep/scripts/03_per_sample_qc_doublet.py`先基础QC，再对通过细胞按样本调用`scanpy.pp.scrublet(expected_doublet_rate=0.04, random_state=0)`；写出`~doublet_removed`的filtered H5AD。配置中存在若干内部参数字段，但该调用没有逐项显式传入，不宣称与旧直接调用逐参数完全相同。<500或异常时返回NaN与False占位；本次549样本单元的实际清单是536 scored/13 skipped_small_sample/0失败/0 capped。

两批的8%上限本次均未触发。expected rate 0.04不是强制删4%的配额；本轮95分位也不是原工具概率或统一阈值。

## 状态字段不等于预测字段

`02_merge_new_cohort/scripts/01_merge_into_2M.py:96–102`将新增分数非缺失映射为`scrublet_success_new_cohort_20260818`，否则not_run，并继承predicted/removed字段。后续core只保留了方法/分数/执行状态，没有原预测与移除列。因此不能从core执行状态推断工具阳性/阴性；本次回到54个原QC对象恢复这两个字段。

<500分支的代码实际未打分，纠正文档“只打分不去除”的措辞。not_run与False占位不能视作阴性；历史调用成功也不保证生物学纯度。

## 实证核验

1. 精确核对冻结V7与core全细胞index，核对冻结checkpoint SHA。
2. 每个来源按当前source_obs_name对齐原QC对象index，恢复原predicted_doublet、doublet_removed、filter_reason及score。检查54个对象均有原预测/移除字段，全4676787行完整匹配，分数含NaN逐项相同。
3. 独立实现不复用第一次的indexer：检查来源键唯一和排序后完整多重集相同，逐项重查flags/分数；从旧2148018行cell_qc原账本重算doublet原因、与其他QC重叠及最终pass数；新增sample_qc_thresholds与原移除汇总交叉核对。
4. 全量按cohort/来源/父类汇总not_run与原预测；按精确index追溯T6、K-09全部及245候选、Mast。嵌套范围不相加。

[完整结果](../inventories/analysis/2026-11-05/v7_doublet_provenance_audit_v1/README.md)与[run](../runs/20261105_v7_doublet_provenance_audit_v1.md)。只发布来源/父类汇总与代码，原细胞条码、样本/患者行、矩阵和完整stdout留私有。PASS只接受过滤历史与成员对应，不能当作doublet识别准确率或最终注释验收。
