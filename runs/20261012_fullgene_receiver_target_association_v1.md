# Run：先验靶表达与发送配体的样本级探索关联

- Run ID：20261012_fullgene_receiver_target_association_v1；parent已审v17/fullgene+620preflight。
- 预声明UTC 2026-10-12T03:40:41.817236+00:00 /UTC+8 2026-10-12T11:40:41.817236+08:00；PREDECLARED_NOT_STARTED。
- scanpy Python3.10.19；620cases×M0/M1/M2=1860记录，3已用于发现source，非独立验证。

## 输入与方法

[完整合同/代码与输入SHA](../inventories/analysis/2026-10-12/fullgene_receiver_target_association_plan_v1.json)固定score/协变量/失败规则/2000resampling和每variant620diagnosticBH。Weightedunsigned prior-target expression不是NicheNetactivity或受体特异通路。原node分开、缺失不当阴性，控制共线不删项。

## 执行结构

提供PANCANCER_PREFLIGHT_ROOT后：

```bash
python fullgene_receiver_target_association_v1.py
python audit_fullgene_receiver_target_association_v1.py
```

真实数值与输入构建审计分别验收；选后diagnosticFDR不称独立确认。私有design/samplekeys留服务器，只有源/情境效应和SHA入仓库。实际首尾与exit待派发后记录；无状态/临床机制发布。

## 实际完成与审计：2026-10-12T12:07:15.946820+08:00

[结果与限制](../inventories/analysis/2026-10-12/target_association_v1/RESULTS.md)：1860记录/566可估计独立数值复算、325控制构建通过；生产与修复后两审计exit0，v1失败保留。M1有241可估计/65诊断q阳性/10稳定区间，仅NODE0002同时q阳性且区间不跨零；M2不稳定不能当阴性。当前COMPLETE_EXPLORATORY_AUDITED_NOT_AXIS_VALIDATED，非独立验证。
