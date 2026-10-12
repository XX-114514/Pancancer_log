# 双细胞、环境RNA与整合输入准备 v1

状态：方法及环境预检。输入为`${COLLECTION_RUN_ROOT}/objects/qc_<dataset>.zarr`，4队列993,312细胞149文库；新输出写入独立`${PREPROCESS_RUN_ROOT}`。原始计数与旧QC对象保持原样。

输入：CXG_5b9b89d2_2e3e_4ad1_8965_666f64540fea（GSE241783）、GSE286827、GSE292046、GSE294109。GSE286827应用已验收patientfix_v2覆盖；稳定ID不变。

双细胞方法及污染校正参数在环境验证后冻结。每文库需保存scores、flags、移除原因、污染率、校正前后计数及marker检查。缺空液滴或单一谱系限制明确登记。

验收：原始counts与上游完全一致，corrected非负有限且不凭空增加表达；新对象ID唯一且为上游ID子集，所有新过滤可追溯；逐文库成功状态与来源AC状态分开。整合准备保留原始基因覆盖，技术／组织／富集分支，禁止缺失基因被解释为真实零表达。

完整执行、失败和证据见[任务日志](../logs/codex/2026/2026-10-11_preprocess-doublet-ambient.md)。本次尚未进行跨研究整合。

## 21:46 UTC+8 方法冻结与范围更新

初始两个完整文库试运行均经独立全细胞验收。正式运行另立`${PREPROCESS_RUN_ROOT}`=prepare_v2；范围9分支792文库，输入3338473细胞（4QC+5预QC），不能作为最终通过数。参数及限制见[方法](../methods/doublet-ambient-preintegration-v1.md)。运行命令：

```bash
python scripts/run_ambient_preparation.py --run-root "${PREPROCESS_RUN_ROOT}" --collection "${COLLECTION_RUN_ROOT}" --rscript "${SCRNA_RSCRIPT}" --code-git-sha "${CODE_GIT_SHA}" --parallel-datasets 4 --library-workers 4
```

实际启动、终止及独立验收结果将追加；现阶段不声称全量完成。

## 2026-10-12T02:12:48.762758+00:00 恢复记录

prepare_v2因组quota errno122终止，保留部分输出、无全量验收。data4组上限现2048G/已用约1495G；恢复输出为另卷私有`${PREPROCESS_RUN_ROOT}`=prepare_v3，输入原位只读。预算782893303100 bytes+500GiB，小于输出卷约1.598TB实际余量。按原9分支792文库恢复；新增必需参数`--storage-plan "${STORAGE_PLAN}"`。不恢复新下载。

## 2026-10-12T10:38:45.707906+08:00 prepare_v3实际运行

启动02:15:17 UTC/10:15:17 UTC+8；代码cdcea29，4来源×4文库，500GiB预留，存储计划已核实。10:27快照44文库/0失败，275653保留细胞、额外9608；完整独立count/ID/标准化审计将在全队列终态执行。final source AC仍0，未开始全局整合。

```bash
python scripts/run_ambient_preparation.py --run-root "${PREPROCESS_RUN_ROOT}" --collection "${COLLECTION_RUN_ROOT}" --rscript "${SCRNA_RSCRIPT}" --code-git-sha cdcea29b45e989cb8d03235b9a92521611b6dd17 --parallel-datasets 4 --library-workers 4 --storage-plan "${STORAGE_PLAN}"
```

[可读快照](../inventories/runtime/2026-10-12/AMBIENT_PROGRESS_v3_checkpoint2.md)不等于实时状态或最终验收。

## 2026-10-12T10:54:17+08:00 第一批完整来源验收／审计I/O优化

- 实测10:53:50（UTC02:53:50）182/792文库，0失败，944249保留细胞、额外38637个预测双细胞，274673严格技术候选。完整运行仍活跃。
- 3个完整来源独立验收：CAR-T 396593、GSE292046 52525、GSE294109 130243，合计579361细胞/92文库；原始审计02:46:16 UTC和缓存审计02:54:17 UTC均PASS；所有来源/文库科学检查、计数、候选124513及基因覆盖诊断完全一致。报告及SHA匹配证据见runtime/2026-10-12。
- 审计缓存只改读取方式：每次最多一源36GiB稀疏参考，较大来源保持原backed读取，校验数学/门槛不变；原始报告保留。生产代码、双细胞和DecontX参数未改。全量阶段审计脚本升级尚待显式版本登记。
- 数量为技术处理与验收，不等于来源别名排重已完成；final AC仍0，未启动跨来源整合、新下载仍暂停。CAR-T产品及富集snRNA用途限制保留。
