# NSCLC固定签名：完整病例级临床家族

## 目标、范围与当前状态

承接[全基因测量与投影](20261016_NSCLC_frozen_programme_measurement_v1.md)，全部相容谱系程序及作者亚群的完整临床家族已执行、独立验收。科学值及全部反例留服务器/授权对话；此记录不发布患者行、学习签名或效应表。

这是治疗后来源内固定签名与作者病理疗效的关联，尚不是细胞对联合验证、直接LR机制、治疗前预测或正式细胞状态。不能将一个上下文MCP1的临床线索移用为同上下文MCP2或某条LR的证明。

## 实际时间线（UTC /UTC+8）

- 09:48:53 /17:48:53：任何临床效应前，写入完整212槽合同及全部病例/亚群/设计；09:48:55预检完成。
- 09:50:05 /17:50:05：唯一全家族正式启动；09:51:08计算终止。无按新结果改字段、亚群、权重或控制。
- 09:53:51 /17:53:51：独立审计启动；09:53:56全部病例/暴露/设计、164可估计HC3、完整BH、23,594项留一及首/中/末三例完整wild/bootstrap通过。
- 09:55后当前任务读取真实审计及临床汇总，取证晚于实际终态。本节点未以运行完成代替互作目标完成。

## 参数和真实验收范围

[完整预冻结方法](../methods/NSCLC-frozen-signature-clinical-family-v1.md)：原作者结局/暴露、病例单位、固定12控制、50病例/每组10/细胞20/设计门槛、20,000次零假设wild、2,000次病例bootstrap、完整212 BH及未调整Mann–Whitney。旧残差/原固定权重分数调整后推断不变性通过。原/残差分数不是两次独立验证。

独立每个实际模型以Statsmodels HC3核系数/SE/t，所有留一以删除观测的线性恒等式核查；全部212计划及排除集不丢失。三例按可估计顺序首/中/末预固定，全部20,000 wild用QR实现重放，全部2,000病例系数bootstrap用SVD秩及QR求解复核，包含区间有效次数不足场景；未冒称全部164模型的每次重采样都独立重放。

## 来源、环境与证据

现有scanpy Python/NumPy/SciPy/Statsmodels，BLAS线程1，源码由各真实launch_receipt SHA固定。临床表逐病例不变性、原全基因64单元独立审计、全部投影分数/训练外推审计先完成。源不是当前免疫图谱的自动细身份映射；过去已有其他临床问题探索，不称项目级完全独立验证。

本轮结束实际环境读取：Python3.10.19、NumPy1.26.4、pandas2.3.3、SciPy1.11.4、Statsmodels0.14.6、Zarr2.18.3。冻结RDS导出继承原R4.3.3，仅导出既有参数，没有重新训练DIALOGUE。

`${PREFLIGHT_ROOT}/NSCLC_frozen_signature_clinical_preflight_20261016_v1/plan_before_any_clinical_effect.json` SHA `60a0146fd102fd1057005ee27012fd6be95bb0d6ae3e746dac8df4e536754a4c`；设计表SHA `4a07a989f0ba4ecb7a807455e720fbb53c535c92576ce3cede62ee6f8c698236`。

`${PREFLIGHT_ROOT}/NSCLC_frozen_signature_clinical_results_20261016_v1/private_all_planned_signature_clinical_results.tsv` SHA `6ab08b1f21a6c4188e52c1f1487adcfbc0ec8526ee1753e1369429b604c767da`；同目录`private_independent_audit.json`状态PASS_ALL_PLANNED_CLINICAL_SCOPE_EXPOSURE_DESIGNS_HC3_BH_LOO_AND_THREE_FULL_WILD_BOOTSTRAP_REPLAYS，审计源码SHA `de1e6aeb4974f5ca267d5760e9742e50706608fb808ef95cd366b82a66481046`。

## 下一问题与停止规则

本完整临床家族停止追加阈值。对所有保留线索先比较固定基因/分数冗余及功能，避免把相关记录当独立发现；再预冻结相容两端的来源内协调及接口问题，保留已经观测的原LR/TF或临床不支持项，不能换MCP编号背书。临床关联仍可能受未知处理/中心及未审病例重复影响，不作为临床预测或因果耐药结论。下一新问题须另立完整家族及对照，不把事后来源选择称独立验证。
