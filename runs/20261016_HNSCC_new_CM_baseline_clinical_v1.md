# 新免疫图谱CM基线临床关联：完整30槽验收

## 用户目标、起点与计划

继续完整细胞互作生物学/临床goal。上一轮真实完成独立乳腺IFN家族、数值审计及全部新CM临床范围审查，已改变下一动作，记progress。当前先核本地/服务器权威状态，保留另一对话DIALOGUE/空间及未跟踪进展review文件。

本轮落实上一预声明四范围/30组件：精确原缩放、固定模型投影真正基线、完整临床家族、独立输入/数值审计和覆盖解释。全程不改原2026年10月CM（v20261010）172/9CM、原HNSCC程序或权威Zarr，不复用混合时点H当基线，不新增表达下载。

## 事件时间线（2026-10-16 UTC；括号UTC+8）

- 首次只读原函数定位的正则转义失败、一次本机glob未匹配；改用普通字符串定位与rg --files，无科学输入或参数变化。
- 原R导出v1在定位freq_normalize AST时遇到命名空间调用的非symbol头，长度逻辑条件失败；未导出模型或算临床。v1错误与目录保留，v2只增加head类型检查，原代码/模型参数保持。
- 04:20:07（12:20:07）：原归一化及NMFns W/S/H/fitted四范围逐项核对exit0。平滑S确实属于原模型，已在首次效果前明确以A=W S作主要固定投影，原W-only仅描述。未重学W/S/缩放。
- 04:20:33（12:20:33）：全部30组件真正基线投影完成；从原模型F重演H、NNLS KKT、原作者患者/取样/治疗/HPV对齐通过，没有主要协变量缺失。唯一患者和范围内患者行分开，范围重叠不能加总为患者。
- 04:22:12（12:22:12）：独立BVLS从原频率重算所有基线投影、原minmax与临床键PASS，效果尚未计算。
- 04:23:07（12:23:07）：完整30患者级临床检验及120描述/LOO真实完成；全家族不按效果删槽，不变二分类/临床切点，不补有利P。
- 04:25:33（12:25:33）：全部30活动/捕获、重复秩型和队列相对投影残差诊断完成，无新P/BH或事后排除。相对残差不当临床R2，少活动患者与不适配不当全癌种生物学否定。
- 04:25:38（12:25:38）：全部30系数/t/20000置换/2000患者bootstrap/区间及统一BH、全部120描述与全LOO独立manual-rank/full-QR复算PASS。
- 04:27:07（12:27:07）：完整私有证据报告exit0，科学结果/患者行/模型矩阵只保留server及获准对话。
- 04:27:28（12:27:28）前重新核既有DIALOGUE：旧chain v1实际04:14:54 failed，输入接受清单/模型尚未生成；其审计为检出比例1%边界失败，另一个对话已有v3输入及chain v2恢复。当前两个新预期child真实live，输入13/14上下文，尚未官方拟合验收。没有重启、派发或更改其他支线。

## 结论与下一阶段

本轮完整临床家族和独立审计已完成，未得到可晋级的新临床CM；所有反例、缺区间、年龄/性别及原W敏感性保留。原巨噬表达程序的单端线索不改写，不能由它宣称两端通讯或模块耐药。

捕获和投影诊断支持来源/队列适配限制需单独解释；不因结果不显著就移除少数活动组件重算FDR，不换rank/归一化追阳性。新图谱CM与原2026年10月CM（v20261010）有明确不同版本/端点，不能把本轮阴性扩展到原9CM或全部癌种。

下一沿已有来源内表达DIALOGUE链、固定患者留出/功能与接收响应推进；实际输入接受、官方基因支持、训练分数重建和留出要分开验收。当前恢复链不等于已发现MCP。合格空间仍受真实测量/患者/端点限制。详见[方法](../methods/HNSCC-frozen-new-CM-clinical-validation-v1.md)。完整goal active，未验证直接细胞互作或临床预测。

## 复现与完整证据

服务器scanpy Python与实际原NMF0.28/R环境；PANCANCER_PROJECT_ROOT/PANCANCER_PREFLIGHT_ROOT提供路径。export_HNSCC_CM_frozen_projection_v1/v2.py、prepare_HNSCC_new_CM_baseline_v1.py、audit_HNSCC_new_CM_baseline_inputs_v1.py、test_HNSCC_new_CM_clinical_v1.py、audit_HNSCC_new_CM_clinical_v1.py、diagnose_HNSCC_new_CM_capture_v1.py、report_HNSCC_new_CM_clinical_private_v1.py均有新版本目录、代码/输入SHA和终态；v1失败明确保留，其余真实exit0。

私有逻辑目录：${PREFLIGHT_ROOT}/HNSCC_new_CM_frozen_model_export_20261016_v1/v2、HNSCC_new_CM_baseline_preparation_20261016_v1、HNSCC_new_CM_clinical_results_20261016_v1、HNSCC_new_CM_capture_diagnostics_20261016_v1、HNSCC_new_CM_clinical_private_evidence_20261016_v1。主要结果SHA a54a83d0e31ef6192d55460b04ed84d1a282e5ff2ec1f0673321af9bd649afb2；基线分数SHA193ac25609c8b5e2bea0bbc5d32aeac2b04e28c565909dbe2df2c7207e09e870；临床键SHA52790e6975069caa5b02a7f7fffc58e03ffa83d54ba022e29924eae4eefa881e；完整报告SHAb657f3688a955b126905fcb0b667230ab1ef07d55586445f98d9ce50932ce59c。

轻量compile、仓库隐私/链接/TSV和diff在发布前检查；只暂存本任务及共享文件本段，保留并发历史。普通commit/push后独立HTTPS核云端完整SHA，实际receipt单列。不公开科学统计、学习权重或患者/细胞矩阵。建议commit：analysis: validate fixed baseline CM clinical family and audit source applicability。
