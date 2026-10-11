# Run：新增四份QC全量技术验收

- ID：20261011_expansion_QC_acceptance_v1，范围为Project_v5公共新增分支，不改变Project_v3身份/空间结果。
- 输入：四份pending QC共993,312细胞/149文库，原始H5AD、过滤账本、阈值/registry、基线数据库和来源SHA。
- 方法：5000行/块逐一核对全部X/counts/原始计数，重新计算QC指标及文库阈值；验证稳定ID、来源顺序、患者/临床registry映射；不重跑Scrublet，不修改生物过滤。
- 环境：服务器scanpy Python3.10；numpy/pandas/scipy/h5py/zarr/anndata；OMP/BLAS各1线程。
- 启动：13:11:23 UTC /21:11:23 UTC+8；证据：[执行计划](../inventories/runtime/2026-10-11/qc_acceptance_execution_v1.json)。
- 代码：[accept_qc_shards.py](../scripts/accept_qc_shards.py)，实际SHA记录于执行证据；运行代码commit待本节点提交。
- 命令：`python -u scripts/accept_qc_shards.py --run "${RUN_ROOT}" --output "${NEW_ACCEPTANCE_AUDIT}"`。
- 输出：${RUN_ROOT}/artifacts/independent_qc_acceptance_20261011_v1.json，进度日志只留服务器；[任务日志](../logs/codex/2026/2026-10-11_QC-acceptance.md)。

本节点为running，不声明技术或来源final acceptance。正式来源验收涉及研究别名、文库重叠、模态与临床用途，不能仅由计算通过自动授予。GSE292046为两患者/8文库的富集snRNA；CAR-T源别名GSE241783且为产品队列。完整结果与适用边界终态后追加。

## 終態：21:15—21:18 UTC+8

全量審計13:15:27.440923 UTC結束，四分片4×36條件PASS/993,312細胞，原始/保留計數一致。13:18:17.691405 UTC元數據v2完成、13:18:52.908477 UTC獨立驗收PASS；7,463細胞覆蓋映射，原矩陣/registry保留。13:18:52.909137 UTC登記technical_QC_ACCEPTED_source_final_pending。[報告](../inventories/runtime/2026-10-11/QC_ACCEPTANCE.md)含數量、用途、分布和全部證據SHA。

啟動節點代碼commit1417788bc00bbf0cfe7c704eb4dfba65dcace772已普通push且獨立HTTPS核驗；數值審計實際腳本SHA與此commit相同。v2修復代碼本終態節點提交，真实已執行且覆核。正式来源新增0，不能由技術验收替代別名/基線文庫審查。
