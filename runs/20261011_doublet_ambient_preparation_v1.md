# 双细胞、环境RNA与整合输入准备 v1

状态：方法及环境预检。输入为`${COLLECTION_RUN_ROOT}/objects/qc_<dataset>.zarr`，4队列993,312细胞149文库；新输出写入独立`${PREPROCESS_RUN_ROOT}`。原始计数与旧QC对象保持原样。

输入：CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea（GSE241783）、GSE286827、GSE292046、GSE294109。GSE286827应用已验收patientfix_v2覆盖；稳定ID不变。

双细胞方法及污染校正参数在环境验证后冻结。每文库需保存scores、flags、移除原因、污染率、校正前后计数及marker检查。缺空液滴或单一谱系限制明确登记。

验收：原始counts与上游完全一致，corrected非负有限且不凭空增加表达；新对象ID唯一且为上游ID子集，所有新过滤可追溯；逐文库成功状态与来源AC状态分开。整合准备保留原始基因覆盖，技术／组织／富集分支，禁止缺失基因被解释为真实零表达。

完整执行、失败和证据见[任务日志](../logs/codex/2026/2026-10-11_preprocess-doublet-ambient.md)。本次尚未进行跨研究整合。
