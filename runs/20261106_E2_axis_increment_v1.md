# 20261106_E2_axis_increment_v1

V5/K9 历史分支的 bulk 轴代理增量校准：全轴零分布与混杂校正。[冻结方案](../methods/E2-axis-increment-frozen-plan-20261106.md) · [SHA 收据](../inventories/analysis/2026-11-06/E2_axis_increment_v1/) · [脚本](../scripts/e2_axis_increment_v1/)

状态：方案和输入已在计算新效应之前冻结，尚未产生结果。运行目录 `${PROJECT_V5_ROOT}/e2_axis_increment_20261106`；患者级输入留在研究环境。这是对已知结果的校准，不是盲法确证。

## 结果（2026-11-06，脚本 v1.0.2）

复现门 24/24 通过（缩放为 global_sd）。SPP1 判定 B（N2a 0.009，N2b 0.29，N3L 0.047，N3R 0.030，轴语义 Holm 0.093；S4 HR 1.18，保留 80%）。LGALS1 判定 C（N2a 0.18，轴语义 Holm 0.76）。特异性未通过：103 条非主轴中 54.4% q<0.05。[汇总表](../inventories/analysis/2026-11-06/E2_axis_increment_v1/results_v1.0.2/)

结论边界：现有 bulk 证据不支持“轴承载预后信号”的表述；SPP1 预后关联本身稳健。这是对已知结果的校准，观察性，不涉及通讯是否发生。

偏差：22 例年龄缺失未在逐癌种和嵌套模型中过滤。v1.0.1 丢弃整个癌种，方向不对；v1.0.2 逐行剔除。修复在看到零分布之后，只影响 robust 步骤，判定未变，门槛未动。S7 不可估计（V7 bulk 打分未构建）。
