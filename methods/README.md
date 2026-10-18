# 方法文档

解释边界与后续资格：[原生程序生物学及空间来源](native-programme-biological-and-spatial-qualification-20261018.md)，公开研究设计、原身份/程序的区别及未下载空间资源门槛。

当前结果入口：[原生全家族解释](NSCLC-native-interface-result-review-v1.md)，依赖完整数值审核，保留所有因子、敏感性及依赖记录，不新增推断。

当前新增入口：[原生接口外部测量机会](NSCLC-native-interface-external-opportunities-v1.md)，仅基因可测性和冻结捕获元数据，不是效果检验或完整验证就绪。

`methods/` 只保存当前认可、仍应执行的方法逻辑。历史方案、被放弃方案及选择理由进入 `decisions/`；某次运行的具体参数进入 `runs/`。

每个方法文件建议包含：

- 目的和适用范围；
- 输入与前置验证；
- 处理步骤；
- 参数与默认值来源；
- 输出；
- 验证；
- 已知限制；
- 相关代码与 run；
- 版本和变更依据。

当前入口：

- [10月14日之后的任务进展与下一阶段](project-progress-and-next-phase-20261016.md)：墨尔本时间线、2026年10月CM程序与空间闭环、临床及外部功能验证依赖。
- [10月12日之后的进展与下一阶段](project-progress-and-next-phase-20261014.md)：按角色和事件时间对照三条主线、依赖与停止条件。
- [CM定义与免疫来源恢复证据合同](CM-definition-and-recovery-v1.md)：固定成员/权重、情境重塑与免疫骨架；区分投影、生态和功能机制。
- [Project_v3 统一工作流](project-v3-unified-workflow.md)

方法发生实质变化时，应同时创建 Codex log；如果变化涉及方案取舍，再新增 ADR。
